"""Aba Exportar: gera arquivos de layout e modelos 3D para ferramentas externas."""
import os
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, ttk

from ..exportacao import FORMATOS, exportar, simuladores
from .comum import RAIZ, AbaBase, AreaTexto, Formulario, protegido, relativo
from .previa import PreviaEstrutura


class AbaExportar(AbaBase):
    TITULO = "Exportar"
    PREVIA = PreviaEstrutura
    RESULTADOS = False

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        caixa = ttk.LabelFrame(self.painel, text="Estrutura", padding=6)
        caixa.pack(fill="x")
        self.sel = tk.StringVar()
        self.combo = ttk.Combobox(caixa, textvariable=self.sel, state="readonly")
        self.combo.pack(fill="x")
        self.combo.bind("<<ComboboxSelected>>", lambda e: self._previa())
        ttk.Label(caixa, text="Use 'Enviar p/ exportação' nas ferramentas de projeto (menu Ferramentas).",
                  foreground="#555", wraplength=300).pack(anchor="w")

        caixa = ttk.LabelFrame(self.painel, text="Formatos", padding=6)
        caixa.pack(fill="x", pady=(8, 0))
        self.marcados = {}
        for rotulo, (_, _, destino) in FORMATOS.items():
            var = tk.BooleanVar(value=True)
            ttk.Checkbutton(caixa, text=rotulo, variable=var).pack(anchor="w")
            ttk.Label(caixa, text=f"   → {destino}", foreground="#555").pack(anchor="w")
            self.marcados[rotulo] = var

        self.form = Formulario(self.painel, [
            ("fmin", "f mín", 0.5, "GHz"),
            ("fmax", "f máx", 10.0, "GHz"),
            ("fsol", "f adaptativa (HFSS)", 10.0, "GHz"),
        ], titulo="Simulação (HFSS / CST)")
        self.form.pack(fill="x", pady=(8, 0))

        caixa = ttk.LabelFrame(self.painel, text="Pasta de saída (relativa ao projeto)", padding=6)
        caixa.pack(fill="x", pady=(8, 0))
        self.pasta = tk.StringVar(value="saidas")
        ttk.Entry(caixa, textvariable=self.pasta).pack(fill="x")
        self.botoes(("Exportar layout", self.exportar), ("Escolher pasta", self._escolher),
                    ("Abrir pasta", self._abrir))

        caixa = ttk.LabelFrame(self.painel, text="Parâmetros S calculados (ADS, AWR, etc.)", padding=6)
        caixa.pack(fill="x", pady=(8, 0))
        self.sel_s = tk.StringVar()
        self.combo_s = ttk.Combobox(caixa, textvariable=self.sel_s, state="readonly")
        self.combo_s.pack(fill="x")
        self.s2p = tk.BooleanVar(value=True)
        self.csv = tk.BooleanVar(value=True)
        ttk.Checkbutton(caixa, text="Touchstone .s2p (ADS: componente S2P)", variable=self.s2p).pack(anchor="w")
        ttk.Checkbutton(caixa, text="CSV (dB e fase)", variable=self.csv).pack(anchor="w")
        ttk.Label(caixa, text="Resultados vêm do último 'Calcular' de cada aba.",
                  foreground="#555", wraplength=300).pack(anchor="w")
        self.botoes(("Exportar parâmetros S", self.exportar_s))

        self.saida = AreaTexto(self.painel, 12)
        self.saida.pack(fill="both", expand=True)
        estado.ouvir(lambda ev: ev in ("layouts", "substrato") and self._atualizar())
        estado.ouvir(lambda ev: ev == "resultados" and self._atualizar_s())
        self._atualizar()
        self._atualizar_s()

    def _caminho_pasta(self):
        p = Path(self.pasta.get())
        return p if p.is_absolute() else RAIZ / p

    def _escolher(self):
        d = filedialog.askdirectory(initialdir=self._caminho_pasta())
        if d:
            self.pasta.set(relativo(d))

    def _abrir(self):
        p = self._caminho_pasta()
        p.mkdir(parents=True, exist_ok=True)
        if hasattr(os, "startfile"):
            os.startfile(p)  # Windows
        else:
            self.estado.msg(str(p))

    def _atualizar(self):
        nomes = list(self.estado.layouts)
        self.combo.configure(values=nomes)
        if nomes and self.sel.get() not in nomes:
            self.sel.set(nomes[-1])
        self._previa()

    def _atualizar_s(self):
        nomes = list(self.estado.resultados)
        self.combo_s.configure(values=nomes)
        if nomes and self.sel_s.get() not in nomes:
            self.sel_s.set(nomes[-1])

    @protegido
    def exportar_s(self):
        nome = self.sel_s.get()
        if nome not in self.estado.resultados:
            raise ValueError("Nenhum resultado calculado ainda. Use 'Calcular' em uma aba de projeto.")
        if not (self.s2p.get() or self.csv.get()):
            raise ValueError("Selecione .s2p e/ou CSV.")
        f, s, z0 = self.estado.resultados[nome]
        pasta = self._caminho_pasta()
        pasta.mkdir(parents=True, exist_ok=True)
        arqs = []
        if self.s2p.get():
            arqs += simuladores.touchstone(pasta / f"{nome}.s2p", f, s, z0, f"Resultado: {nome}")
        if self.csv.get():
            arqs += simuladores.csv_s(pasta / f"{nome}_S.csv", f, s)
        self.saida.escrever(f"Parâmetros S (Z0 = {z0:g} Ω):\n" + "\n".join(relativo(a) for a in arqs))
        self.estado.msg(f"Parâmetros S de '{nome}' exportados.")

    def _previa(self):
        lay = self.estado.layouts.get(self.sel.get())
        if lay:
            self.previa.mostrar(lay)
            info = "\n".join(f"  {k}: {v}" for k, v in lay.info.items())
            s = lay.substrato
            self.saida.escrever(f"{lay.nome}\nSubstrato: {s.nome}, εr={s.er:g}, h={s.h_mm:g} mm\n{info}")

    @protegido
    def exportar(self):
        lay = self.estado.layouts.get(self.sel.get())
        if not lay:
            raise ValueError("Nenhuma estrutura selecionada. Gere uma nas ferramentas de projeto (menu Ferramentas).")
        formatos = [r for r, v in self.marcados.items() if v.get()]
        if not formatos:
            raise ValueError("Selecione pelo menos um formato.")
        v = self.form.get()
        if v["fmax"] <= v["fmin"]:
            raise ValueError("Faixa de frequência inválida.")
        arqs = exportar(lay, self._caminho_pasta(), formatos,
                        dict(f_min_ghz=v["fmin"], f_max_ghz=v["fmax"], f_sol_ghz=v["fsol"]))
        self.saida.escrever("Arquivos gerados:\n" + "\n".join(relativo(a) for a in arqs))
        self.estado.msg(f"{len(arqs)} arquivos exportados em {relativo(self._caminho_pasta())}")
