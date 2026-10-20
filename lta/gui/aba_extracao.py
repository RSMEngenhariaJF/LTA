"""Aba Extração de parâmetros: modelo de circuito a partir de S medido ou simulado (Apêndice G)."""
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, ttk

import numpy as np

from ..core import abcd
from ..core import extracao as ex
from ..exportacao.simuladores import ler_touchstone
from .comum import PASTA_SAIDAS, AbaBase, AreaTexto, Formulario, desenhar_smith, protegido, relativo

ARQUIVO = "Arquivo .s2p…"
UNIDADES = {"L": ("nH", 1e9), "C": ("pF", 1e12)}

# Exemplos com os valores extraídos publicados no livro (dados sintéticos do próprio modelo)
EXEMPLOS = {
    "Exemplo: Tabela 3.3, CSRR + gap": ("CSRR + gap (microstrip)",
                                         dict(L=5.78e-9, Cg=0.31e-12, C=6.43e-12, Lc=1.49e-9, Cc=3.10e-12),
                                         (0.5e9, 5e9)),
    "Exemplo: Tabela 3.3, CSRR sem gap": ("CSRR (microstrip)",
                                           dict(L=5.78e-9, C=0.95e-12, Lc=1.65e-9, Cc=2.81e-12), (0.5e9, 4e9)),
    "Exemplo: Tabela 3.1, SRR + fita (CPW)": ("SRR + fita (CPW)",
                                               dict(C=2.07e-12, Lp=1.38e-9, L=4.91e-9, Ls=1.68e-9, Cs=3.90e-12),
                                               (0.5e9, 5e9)),
}


def _unid(chave):
    return UNIDADES[chave[0]]


class AbaExtracao(AbaBase):
    TITULO = "Extração de parâmetros (Ap. G)"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        caixa = ttk.LabelFrame(self.painel, text="Dados (uma célula simétrica)", padding=6)
        caixa.pack(fill="x")
        self.fonte = tk.StringVar(value=next(iter(EXEMPLOS)))
        self.combo = ttk.Combobox(caixa, textvariable=self.fonte, state="readonly", width=34)
        self.combo.pack(fill="x")
        self.combo.bind("<<ComboboxSelected>>", lambda e: self._fonte_mudou())
        self.lbl_dados = ttk.Label(caixa, text="", foreground="#555", wraplength=300)
        self.lbl_dados.pack(anchor="w", pady=(2, 0))
        self.form = Formulario(self.painel, [
            ("modelo", "Modelo", tuple(ex.MODELOS), ""),
            ("fmin", "Usar de", 0.0, "GHz"),
            ("fmax", "até (0 = tudo)", 0.0, "GHz"),
        ], titulo="Modelo de circuito")
        self.form.pack(fill="x", pady=(8, 0))
        self.man = Formulario(self.painel, [(k, f"{k}: {d}", 0.0, "GHz") for k, d in ex.FREQUENCIAS.items()],
                              titulo="Frequências manuais (0 = automático)")
        self.man.pack(fill="x", pady=(8, 0))
        self.botoes(("Extrair", self.extrair), ("Exportar S2P do modelo", self.exportar_modelo))
        self.saida = AreaTexto(self.painel, 18)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Fluxo do livro (Apêndice G):\n"
            "1. Simule no EM (sem perdas) ou meça UMA célula.\n"
            "2. Exporte .s2p e carregue aqui (ou use o resultado\n   de outra aba).\n"
            "3. Escolha o modelo e clique em Extrair.\n\n"
            "As frequências características saem da carta de\n"
            "Smith (interseções com r = 1 e g = 1), do zero de\n"
            "transmissão e de βl = ±π/2 (de cos βl = (A+D)/2).\n"
            "Com vários candidatos, fica a combinação de menor\n"
            "erro contra os dados.")
        self._dados = None
        self._ultimo = None
        estado.ouvir(lambda ev: ev == "resultados" and self._atualizar_fontes())
        self._atualizar_fontes()
        self._fonte_mudou()

    # ------------------------------------------------------------------ dados
    def _atualizar_fontes(self):
        res = [f"Resultado da aba: {n}" for n in self.estado.resultados]
        self.combo.configure(values=list(EXEMPLOS) + [ARQUIVO] + res)

    def _fonte_mudou(self):
        fonte = self.fonte.get()
        try:
            if fonte in EXEMPLOS:
                modelo, p, (fa, fb) = EXEMPLOS[fonte]
                f = np.linspace(fa, fb, 3001)
                self._dados = (f, ex.s_modelo(modelo, p, f), 50.0, fonte)
                self.form.set({"modelo": modelo})
            elif fonte == ARQUIVO:
                caminho = filedialog.askopenfilename(initialdir=PASTA_SAIDAS, filetypes=[("Touchstone", "*.s2p")])
                if not caminho:
                    return
                f, s, z0 = ler_touchstone(caminho)
                self._dados = (f, s, z0, relativo(caminho))
                self.fonte.set(f"Arquivo: {Path(caminho).name}")
            elif fonte.startswith("Resultado da aba: "):
                f, s, z0 = self.estado.resultados[fonte.split(": ", 1)[1]]
                self._dados = (f, s, z0, fonte)
        except ValueError as e:
            self._dados = None
            self.lbl_dados.configure(text=f"⚠ {e}", foreground="#b0122b")
            return
        if self._dados:
            f, s, z0, nome = self._dados
            self.lbl_dados.configure(text=f"{nome}: {len(f)} pontos, {f[0] / 1e9:.3f} a {f[-1] / 1e9:.3f} GHz, "
                                          f"Z0 = {z0:g} Ω", foreground="#555")

    # ------------------------------------------------------------------ extração
    @protegido
    def extrair(self):
        if not self._dados:
            raise ValueError("Escolha uma fonte de dados.")
        f, s, z0, nome = self._dados
        v = self.form.get()
        sel = (f >= v["fmin"] * 1e9) & ((f <= v["fmax"] * 1e9) if v["fmax"] > 0 else True)
        if sel.sum() < 20:
            raise ValueError("Poucos pontos na faixa escolhida.")
        f, s = f[sel], s[sel]
        manual = {k: x * 1e9 for k, x in self.man.get().items() if x > 0}
        p, cond, err = ex.extrair(v["modelo"], f, s, z0, manual)
        sm = ex.s_modelo(v["modelo"], p, f, z0)
        self._ultimo = (f, sm, z0, v["modelo"], p)
        self.estado.publicar_resultado("modelo_extraido", f, sm, z0)
        self._graficos(f, s, sm, z0, cond)
        mod = ex.MODELOS[v["modelo"]]
        linhas = [f"{mod.nome} ({mod.figura})", mod.descricao, "", "Elementos extraídos:"]
        for k in mod.parametros:
            u, esc = _unid(k)
            linhas.append(f"  {k:3s} = {p[k] * esc:.4g} {u}")
        linhas += ["", "Condições usadas:"] + [f"  {k}: {x:.4f}" for k, x in cond.items()]
        linhas += ["", f"Erro rms |S_modelo - S_dados| = {err:.4f}",
                   "(< 0,02: ótimo; > 0,1: modelo inadequado ou faixa",
                   " de frequência grande demais para o modelo)", "",
                   f"Dados: {nome}"]
        self.saida.escrever("\n".join(linhas))

    def _graficos(self, f, s, sm, z0, cond):
        fg = f / 1e9
        ax = self.grafico.eixos(2, 2)
        desenhar_smith(ax[0])
        ax[0].plot(s[:, 0, 0].real, s[:, 0, 0].imag, color="C0", lw=1.5, label="Dados")
        ax[0].plot(sm[:, 0, 0].real, sm[:, 0, 0].imag, color="C1", ls=":", lw=1.5, label="Modelo")
        ax[0].set_title("S11 na carta de Smith", fontsize=9)
        ax[0].legend(fontsize=7, loc="lower left")
        for i, (nome, cor) in enumerate((("S21", "C0"), ("S11", "C3"))):
            idx = (1, 0) if nome == "S21" else (0, 0)
            ax[1].plot(fg, abcd.db(s[:, idx[0], idx[1]]), color=cor, label=f"|{nome}| dados")
            ax[1].plot(fg, abcd.db(sm[:, idx[0], idx[1]]), color=cor, ls=":", label=f"|{nome}| modelo")
        ax[1].set_ylim(max(-60, ax[1].get_ylim()[0]), 2)
        ax[1].set_ylabel("|S| [dB]")
        ax[1].set_title("Módulo", fontsize=9)
        ax[2].plot(fg, np.degrees(np.angle(s[:, 1, 0])), label="dados")
        ax[2].plot(fg, np.degrees(np.angle(sm[:, 1, 0])), ls=":", label="modelo")
        ax[2].set_ylabel("∠S21 [°]")
        ax[2].set_title("Fase de transmissão", fontsize=9)
        for dados, rot, ls in ((s, "dados", "-"), (sm, "modelo", ":")):
            m = abcd.s_para_abcd(dados, z0)
            c = (0.5 * (m[:, 0, 0] + m[:, 1, 1])).real
            ax[3].plot(np.where(np.abs(c) <= 1, np.degrees(np.arccos(np.clip(c, -1, 1))), np.nan), fg, ls,
                       label=rot)
        ax[3].set_xlabel("|βl| [°]")
        ax[3].set_ylabel("f [GHz]")
        ax[3].set_title("Dispersão da célula, cos βl = (A+D)/2", fontsize=9)
        marcas = [("π/2" if k.startswith("f(βl") else k.split(" ")[0], x) for k, x in cond.items() if "[GHz]" in k]
        for a in ax[1:3]:
            for rot, x in marcas:
                a.axvline(x, color="gray", ls=":", lw=0.7)
                a.annotate(rot, (x, a.get_ylim()[1]), fontsize=7, ha="center", va="bottom", color="gray",
                           annotation_clip=False)
            a.set_xlabel("f [GHz]")
        for a in ax[1:]:
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()

    def exportar_modelo(self):
        if not self._ultimo:
            self.extrair()
        if self._ultimo:
            f, sm, z0, modelo, p = self._ultimo
            txt = " ".join(f"{k}={p[k] * _unid(k)[1]:.4g}{_unid(k)[0]}" for k in ex.MODELOS[modelo].parametros)
            self.salvar_s2p("modelo_extraido", f, sm, z0, f"Modelo extraido {modelo}: {txt}")
