"""Janela principal: navegação por menu, uma ferramenta visível por vez.

Para adicionar ferramentas, inclua a classe (subclasse de AbaBase) no grupo
desejado de GRUPOS; para um capítulo novo, acrescente um grupo. O menu
"Ferramentas", os atalhos Ctrl+número e o título são montados a partir daqui.
"""
import os
import tkinter as tk
from pathlib import Path
from tkinter import messagebox, ttk

from .aba_aplicacoes_ebg import AbaAplicacoesEBG
from .aba_carga import AbaCarga
from .aba_cl_loaded import AbaCLLoaded
from .aba_celula import AbaCelula
from .aba_crlh import AbaCRLH
from .aba_dispositivos import AbaDispositivos
from .aba_csrr import AbaCSRR
from .aba_ebg import AbaEBG
from .aba_exportar import AbaExportar
from .aba_extracao import AbaExtracao
from .aba_floquet import AbaFloquet
from .aba_multibanda import AbaMultibanda
from .aba_substrato import AbaSubstrato
from .comum import PASTA_SAIDAS, RAIZ, Estado

# (nome do grupo no menu, [ferramentas])
GRUPOS = [
    ("Cap. 2 – Linhas periódicas", [AbaCelula, AbaCarga, AbaEBG, AbaFloquet, AbaAplicacoesEBG]),
    ("Cap. 3 – Metamateriais", [AbaCRLH, AbaCLLoaded, AbaCSRR, AbaExtracao]),
    ("Cap. 4 – Aplicações", [AbaMultibanda, AbaDispositivos]),
]
# ferramentas gerais, fora dos capítulos (acessadas pelo menu Arquivo/Exportar e pela barra do substrato)
GERAIS = [AbaSubstrato, AbaExportar]

# lista plana mantida para scripts e testes (ordem de criação)
ABAS = [AbaSubstrato] + [c for _, cs in GRUPOS for c in cs] + [AbaExportar]

VERSAO = "1.5"


class Aplicacao:
    def __init__(self, raiz):
        self.raiz = raiz
        raiz.title("LTA – Linhas de Transmissão Artificiais")
        try:
            ttk.Style().theme_use("vista")
        except tk.TclError:
            pass
        self.estado = Estado()
        self.status = tk.StringVar(value="Pronto.")
        self.estado.status = self.status.set

        self._barra_substrato()
        self.titulo = tk.StringVar()
        ttk.Label(raiz, textvariable=self.titulo, font=("Segoe UI", 11, "bold"),
                  padding=(8, 4), foreground="#1F3A5F").pack(fill="x")
        ttk.Label(raiz, textvariable=self.status, relief="sunken", anchor="w",
                  padding=(6, 2)).pack(fill="x", side="bottom")
        self.area = ttk.Frame(raiz)
        self.area.pack(fill="both", expand=True)
        self.area.rowconfigure(0, weight=1)
        self.area.columnconfigure(0, weight=1)

        # cria todas as ferramentas uma vez (cada uma guarda seu estado) e empilha
        self.ferramentas = {}
        self.caminho = {}
        for cls in ABAS:
            f = cls(self.area, self.estado)
            f.grid(row=0, column=0, sticky="nsew")
            self.ferramentas[cls] = f
        for grupo, classes in GRUPOS:
            for cls in classes:
                self.caminho[cls] = f"{grupo} › {cls.TITULO}"
        self.caminho[AbaSubstrato] = "Substrato"
        self.caminho[AbaExportar] = "Exportar"

        self.atual = tk.StringVar()
        self._menu()
        self.mostrar(GRUPOS[0][1][0])

    # ------------------------------------------------------------------ partes
    def _barra_substrato(self):
        barra = ttk.Frame(self.raiz, padding=(8, 4))
        barra.pack(fill="x")
        self.txt_sub = tk.StringVar()
        ttk.Label(barra, text="Substrato:", font=("Segoe UI", 9, "bold")).pack(side="left")
        ttk.Label(barra, textvariable=self.txt_sub).pack(side="left", padx=6)
        ttk.Button(barra, text="Alterar…", command=lambda: self.mostrar(AbaSubstrato)).pack(side="left")
        ttk.Separator(self.raiz).pack(fill="x")

        def atualizar(ev=None):
            s = self.estado.substrato
            self.txt_sub.set(f"{s.nome}   εr = {s.er:g}   h = {s.h_mm:g} mm   t = {s.t_um:g} µm   "
                             f"tan δ = {s.tand:g}")
        self.estado.ouvir(lambda ev: ev == "substrato" and atualizar())
        atualizar()

    def _menu(self):
        barra = tk.Menu(self.raiz)
        self.raiz.config(menu=barra)

        arq = tk.Menu(barra, tearoff=False)
        arq.add_command(label="Substrato…", command=lambda: self.mostrar(AbaSubstrato))
        arq.add_separator()
        arq.add_command(label="Abrir pasta de saídas", command=self._abrir_saidas)
        arq.add_separator()
        arq.add_command(label="Sair", command=self.sair, accelerator="Alt+F4")
        barra.add_cascade(label="Arquivo", menu=arq)

        ferr = tk.Menu(barra, tearoff=False)
        n = 1
        for grupo, classes in GRUPOS:
            sub = tk.Menu(ferr, tearoff=False)
            for cls in classes:
                atalho = f"Ctrl+{n}" if n <= 9 else None
                sub.add_radiobutton(label=cls.TITULO, variable=self.atual, value=self.caminho[cls],
                                    command=lambda c=cls: self.mostrar(c), accelerator=atalho)
                if n <= 9:
                    self.raiz.bind_all(f"<Control-Key-{n}>", lambda e, c=cls: self.mostrar(c))
                n += 1
            ferr.add_cascade(label=grupo, menu=sub)
        barra.add_cascade(label="Ferramentas", menu=ferr)

        exp = tk.Menu(barra, tearoff=False)
        exp.add_command(label="Layout e parâmetros S…", command=lambda: self.mostrar(AbaExportar),
                        accelerator="Ctrl+E")
        self.raiz.bind_all("<Control-e>", lambda e: self.mostrar(AbaExportar))
        barra.add_cascade(label="Exportar", menu=exp)

        aj = tk.Menu(barra, tearoff=False)
        aj.add_command(label="Material didático (PDF)",
                       command=lambda: self._abrir(RAIZ / "doc" / "Manual_LTA_Linhas_Periodicas.pdf"))
        aj.add_command(label="Material didático (Word)",
                       command=lambda: self._abrir(RAIZ / "doc" / "Manual_LTA_Linhas_Periodicas.docx"))
        aj.add_command(label="Manual técnico", command=lambda: self._abrir(RAIZ / "doc" / "01_manual.md"))
        aj.add_separator()
        aj.add_command(label="Atalhos de teclado", command=self._atalhos)
        aj.add_command(label="Sobre", command=self._sobre)
        barra.add_cascade(label="Ajuda", menu=aj)

    # ------------------------------------------------------------------ ações
    def mostrar(self, cls):
        """Traz a ferramenta para frente e atualiza título e menu."""
        self.ferramentas[cls].tkraise()
        self.atual.set(self.caminho[cls])
        self.titulo.set(self.caminho[cls])

    def ferramenta(self, cls):
        return self.ferramentas[cls]

    def _abrir(self, caminho):
        caminho = Path(caminho)
        if not caminho.exists():
            messagebox.showinfo("Arquivo não encontrado", f"{caminho.name} ainda não foi gerado.")
        elif hasattr(os, "startfile"):
            os.startfile(caminho)

    def _abrir_saidas(self):
        PASTA_SAIDAS.mkdir(exist_ok=True)
        self._abrir(PASTA_SAIDAS)

    def _atalhos(self):
        linhas, n = [], 1
        for grupo, classes in GRUPOS:
            for cls in classes:
                linhas.append(f"Ctrl+{n}\t{grupo} › {cls.TITULO}")
                n += 1
        linhas.append("Ctrl+E\tExportar")
        messagebox.showinfo("Atalhos de teclado", "\n".join(linhas))

    def _sobre(self):
        messagebox.showinfo("Sobre", f"LTA – Linhas de Transmissão Artificiais, versão {VERSAO}\n"
                                     "Base teórica: F. Martín, Artificial Transmission Lines for RF and "
                                     "Microwave Applications, Wiley, 2015.\n\n"
                                     "Resultados de modelos de circuito: validar por simulação EM e medição.")

    def sair(self):
        self.raiz.quit()
        self.raiz.destroy()


def main():
    raiz = tk.Tk()
    raiz.geometry("1350x860")
    app = Aplicacao(raiz)
    raiz.protocol("WM_DELETE_WINDOW", app.sair)
    raiz.mainloop()
