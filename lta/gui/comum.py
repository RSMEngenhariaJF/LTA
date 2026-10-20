"""Componentes compartilhados da interface: estado, formulários, gráficos."""
import functools
import tkinter as tk
import traceback
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Polygon, Rectangle

from ..core.microstrip import Substrato
from ..exportacao import simuladores

RAIZ = Path(__file__).resolve().parents[2]   # pasta linhas_periodicas/
PASTA_SAIDAS = RAIZ / "saidas"


def relativo(caminho):
    """Mostra caminhos relativos à raiz do projeto, quando possível."""
    try:
        return str(Path(caminho).resolve().relative_to(RAIZ))
    except ValueError:
        return str(caminho)


class Estado:
    """Dados compartilhados entre as abas."""

    def __init__(self):
        self.substrato = Substrato()
        self.layouts = {}       # nome -> Layout
        self.resultados = {}    # nome -> (f_hz, S(Nf,2,2), z0)
        self._ouvintes = []
        self.status = None      # função para escrever na barra de status

    def ouvir(self, func):
        self._ouvintes.append(func)

    def notificar(self, evento):
        for f in self._ouvintes:
            f(evento)

    def publicar_layout(self, layout):
        self.layouts[layout.nome] = layout
        self.notificar("layouts")

    def publicar_resultado(self, nome, f, s, z0):
        self.resultados[nome] = (f, s, z0)
        self.notificar("resultados")

    def definir_substrato(self, sub):
        self.substrato = sub
        self.notificar("substrato")

    def msg(self, texto):
        if self.status:
            self.status(texto)


def protegido(metodo):
    """Mostra erros de cálculo numa caixa de diálogo em vez de travar a GUI."""
    @functools.wraps(metodo)
    def envolto(*a, **k):
        try:
            return metodo(*a, **k)
        except ValueError as e:
            messagebox.showerror("Erro", str(e))
        except Exception as e:  # noqa: BLE001
            traceback.print_exc()
            messagebox.showerror("Erro inesperado", f"{type(e).__name__}: {e}")
    return envolto


class Formulario(ttk.Frame):
    """Grade de campos. campos: [(chave, rótulo, padrão, unidade)].

    Se o padrão for tupla/lista, o campo vira combobox com essas opções
    (o primeiro item é o valor inicial). Strings viram campo de texto.
    """

    def __init__(self, mestre, campos, titulo=None):
        super().__init__(mestre)
        self.vars, self.tipos, self.widgets = {}, {}, {}
        alvo = ttk.LabelFrame(self, text=titulo, padding=6) if titulo else self
        if titulo:
            alvo.pack(fill="x")
        for i, (chave, rotulo, padrao, unidade) in enumerate(campos):
            # rótulos longos quebram em várias linhas para não empurrar os campos para fora do painel
            ttk.Label(alvo, text=rotulo, wraplength=150).grid(row=i, column=0, sticky="w", padx=(0, 4), pady=1)
            if isinstance(padrao, (tuple, list)):
                var = tk.StringVar(value=padrao[0])
                w = ttk.Combobox(alvo, textvariable=var, values=list(padrao), state="readonly", width=16)
                self.tipos[chave] = str
            else:
                var = tk.StringVar(value=str(padrao))
                w = ttk.Entry(alvo, textvariable=var, width=10)
                self.tipos[chave] = type(padrao)
            w.grid(row=i, column=1, sticky="ew", pady=1)
            ttk.Label(alvo, text=unidade, foreground="#555").grid(row=i, column=2, sticky="w", padx=(4, 0))
            self.vars[chave] = var
            self.widgets[chave] = w
        alvo.columnconfigure(1, weight=1)

    def get(self):
        out = {}
        for k, v in self.vars.items():
            tipo = self.tipos[k]
            txt = v.get().strip().replace(",", ".") if tipo in (int, float) else v.get()
            try:
                out[k] = int(float(txt)) if tipo is int else (float(txt) if tipo is float else txt)
            except ValueError:
                raise ValueError(f"Valor inválido no campo '{k}': {v.get()!r}") from None
        return out

    def set(self, valores):
        """Preenche campos sem disparar os observadores (alteração programática)."""
        self._mudo = True
        try:
            for k, v in valores.items():
                if k in self.vars:
                    self.vars[k].set(str(v))
        finally:
            self._mudo = False

    def ao_mudar(self, chave, func):
        self.vars[chave].trace_add("write", lambda *_: getattr(self, "_mudo", False) or func())

    def observar(self, func):
        """Chama func quando o usuário altera qualquer campo."""
        for chave in self.vars:
            self.ao_mudar(chave, func)


class OpcoesModelo(ttk.LabelFrame):
    """Caixas de seleção do modelo físico da microstrip: perdas e dispersão."""

    def __init__(self, mestre):
        super().__init__(mestre, text="Modelo físico da linha", padding=6)
        self.perdas = tk.BooleanVar(value=False)
        self.dispersao = tk.BooleanVar(value=False)
        ttk.Checkbutton(self, text="Perdas: condutor (σ, t, rugosidade) e dielétrico (tan δ)",
                        variable=self.perdas).pack(anchor="w")
        ttk.Checkbutton(self, text="Dispersão de εeff e Z0 (Kirschning-Jansen)",
                        variable=self.dispersao).pack(anchor="w")

    def get(self):
        return self.perdas.get(), self.dispersao.get()

    def set(self, perdas=False, dispersao=False):
        self.perdas.set(perdas)
        self.dispersao.set(dispersao)

    def descrever(self, sub):
        p, d = self.get()
        if not (p or d):
            return "Modelo: quase-estático, sem perdas."
        partes = []
        if p:
            partes.append(f"perdas (σ = {sub.sigma / 1e6:g} MS/m, t = {sub.t_um:g} µm, "
                          f"rug. = {sub.rug_um:g} µm, tan δ = {sub.tand:g})")
        if d:
            partes.append("dispersão Kirschning-Jansen")
        return "Modelo: " + " + ".join(partes) + "."


class Grafico(ttk.Frame):
    """Figura matplotlib com barra de ferramentas (zoom, salvar imagem...)."""

    def __init__(self, mestre):
        super().__init__(mestre)
        self.fig = Figure(figsize=(7, 5), dpi=100, layout="constrained")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        barra = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
        barra.update()
        barra.pack(side="bottom", fill="x")
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def eixos(self, linhas, colunas=1, **kw):
        self.fig.clear()
        ax = self.fig.subplots(linhas, colunas, squeeze=False, **kw)
        return ax.ravel()

    def desenhar(self):
        self.canvas.draw_idle()


class AreaTexto(ttk.Frame):
    def __init__(self, mestre, altura=10):
        super().__init__(mestre)
        self.txt = tk.Text(self, height=altura, width=42, wrap="word", font=("Consolas", 9))
        sb = ttk.Scrollbar(self, command=self.txt.yview)
        self.txt.configure(yscrollcommand=sb.set, state="disabled")
        self.txt.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")

    def escrever(self, texto):
        self.txt.configure(state="normal")
        self.txt.delete("1.0", "end")
        self.txt.insert("end", texto)
        self.txt.configure(state="disabled")


class AbaBase(ttk.Frame):
    """Aba com painel de controle à esquerda e gráficos à direita.

    Para criar uma ferramenta nova: herdar desta classe, definir TITULO e montar
    os controles em self.painel; depois incluí-la num grupo de GRUPOS em gui/app.py.

    PREVIA: classe de pré-visualização (ex.: PreviaEstrutura) exibida na
    sub-aba "Estrutura"; a aba implementa atualizar_previa() e chama
    self.observar_previa(form) para atualizar enquanto o usuário edita.
    RESULTADOS: se False, não cria a sub-aba de gráficos de resultado.
    """
    TITULO = "Aba"
    PREVIA = None
    RESULTADOS = True

    def __init__(self, mestre, estado: Estado):
        super().__init__(mestre)
        self.estado = estado
        esquerda = ttk.Frame(self, padding=6)
        esquerda.pack(side="left", fill="y")
        canvas = tk.Canvas(esquerda, width=330, highlightthickness=0)
        sb = ttk.Scrollbar(esquerda, orient="vertical", command=canvas.yview)
        self.painel = ttk.Frame(canvas)
        self.painel.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=self.painel, anchor="nw", width=330)
        canvas.configure(yscrollcommand=sb.set)
        canvas.pack(side="left", fill="y")
        sb.pack(side="right", fill="y")
        self.previa = None
        self._agendado = None
        if self.PREVIA is None:
            self.grafico = Grafico(self)
            self.grafico.pack(side="right", fill="both", expand=True)
            return
        self.direita = ttk.Notebook(self)
        self.direita.pack(side="right", fill="both", expand=True)
        self.previa = self.PREVIA(self.direita)
        self.direita.add(self.previa, text="Estrutura (pré-visualização)")
        if self.RESULTADOS:
            self.grafico = Grafico(self.direita)
            self.direita.add(self.grafico, text="Resultados")

    def ver_estrutura(self):
        if self.previa is not None:
            self.direita.select(self.previa)

    def ver_resultados(self):
        if self.previa is not None and self.RESULTADOS:
            self.direita.select(self.grafico)

    def observar_previa(self, *forms):
        for f in forms:
            f.observar(self._agendar_previa)

    def _agendar_previa(self):
        """Espera o usuário parar de digitar (~0,4 s) antes de redesenhar."""
        if self._agendado:
            self.after_cancel(self._agendado)
        self._agendado = self.after(400, self._previa_por_edicao)

    def _previa_por_edicao(self):
        self._agendado = None
        self.ver_estrutura()
        self.atualizar_previa()

    def atualizar_previa(self):
        """Sobrescrever: gera a geometria e chama self.previa.mostrar(...)."""

    def botoes(self, *itens):
        f = ttk.Frame(self.painel)
        f.pack(fill="x", pady=4)
        for i, (texto, cmd) in enumerate(itens):
            ttk.Button(f, text=texto, command=cmd).grid(row=i // 2, column=i % 2, sticky="ew", padx=2, pady=2)
        f.columnconfigure((0, 1), weight=1)
        return f

    def salvar_s2p(self, nome, f, s, z0, comentario=""):
        PASTA_SAIDAS.mkdir(exist_ok=True)
        caminho = filedialog.asksaveasfilename(
            initialdir=PASTA_SAIDAS, initialfile=f"{nome}.s2p", defaultextension=".s2p",
            filetypes=[("Touchstone 2 portas", "*.s2p")])
        if caminho:
            simuladores.touchstone(Path(caminho), f, s, z0, comentario)
            self.estado.msg(f"Touchstone salvo em {relativo(caminho)}")


def desenhar_layout(ax, layout):
    """Vista superior: cobre do topo (laranja), furos no terra (tracejado), placa."""
    x0, y0, x1, y1 = layout.placa
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc="#dfe8d8", ec="#556b2f", lw=1))
    for cx, cy, r in layout.furos_terra:
        ax.add_patch(Circle((cx, cy), r, fc="white", ec="#1f4e9a", ls="--", lw=1))
    for p in layout.topo:
        ax.add_patch(Polygon(p, closed=True, fc="#d98c1f", ec="#8a5208", lw=0.6, alpha=0.9))
    for porta in layout.portas:
        ax.annotate(f"P{layout.portas.index(porta) + 1}", (porta.x, porta.y),
                    xytext=(-14 if porta.lado == "xmin" else 4, 6), textcoords="offset points", fontsize=8)
    ax.set_xlim(x0 - 1, x1 + 1)
    ax.set_ylim(y0 - 1, y1 + 1)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title(f"Layout: {layout.nome}  ({layout.comprimento:.2f} x {layout.largura:.2f} mm)", fontsize=9)


def desenhar_smith(ax):
    """Carta de Smith simplificada, com os círculos r = 1 e g = 1 destacados."""
    t = np.linspace(0, 2 * np.pi, 361)
    ax.plot(np.cos(t), np.sin(t), color="k", lw=0.8)
    for r in (0.2, 0.5, 2, 5):
        c, raio = r / (1 + r), 1 / (1 + r)
        ax.plot(c + raio * np.cos(t), raio * np.sin(t), color="#ccc", lw=0.5)
    for x in (0.2, 0.5, 1, 2, 5):
        for sx in (1, -1):
            z = 1j * sx * x + np.logspace(-3, 3, 300)
            g = (z - 1) / (z + 1)
            ax.plot(g.real, g.imag, color="#ccc", lw=0.5)
    ax.plot(0.5 + 0.5 * np.cos(t), 0.5 * np.sin(t), color="C2", lw=1, ls="--", label="r = 1")
    ax.plot(-0.5 + 0.5 * np.cos(t), 0.5 * np.sin(t), color="C4", lw=1, ls="--", label="g = 1")
    ax.axhline(0, color="#ccc", lw=0.5)
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])


def salvar_csv_s(aba, nome, f, s):
    """Pergunta o arquivo e grava S de N portas em CSV (dB e graus; separador ';')."""
    PASTA_SAIDAS.mkdir(exist_ok=True)
    caminho = filedialog.asksaveasfilename(initialdir=PASTA_SAIDAS, initialfile=f"{nome}_S.csv",
                                           defaultextension=".csv", filetypes=[("CSV", "*.csv")])
    if not caminho:
        return
    n = s.shape[1]
    cab = ["f_GHz"] + [f"S{i + 1}{j + 1}_{q}" for i in range(n) for j in range(n) for q in ("dB", "graus")]
    linhas = [";".join(cab)]
    for fi, si in zip(f, s):
        vals = [f"{fi / 1e9:.9g}"]
        for i in range(n):
            for j in range(n):
                vals += [f"{20 * np.log10(max(abs(si[i, j]), 1e-15)):.6f}", f"{np.degrees(np.angle(si[i, j])):.4f}"]
        linhas.append(";".join(vals))
    Path(caminho).write_text("\n".join(linhas) + "\n", encoding="ascii")
    aba.estado.msg(f"S salvo em {relativo(caminho)}")

