"""Aba Ressonadores de linha (Pozar, seções 6.2 e 6.6): Q0, RLC equivalente e acoplamento por gap."""
import tkinter as tk
from tkinter import ttk

import numpy as np

from ..core import ressonadores as rs
from ..core.descontinuidades import texto_limites
from ..geometria import circuitos
from .comum import AbaBase, AreaTexto, Formulario, desenhar_smith, protegido
from .previa import PreviaEstrutura

ACOPLAMENTO = ("Crítico (Qe = Q0)", "Capacitância dada")


class AbaRessonadores(AbaBase):
    TITULO = "Ressonadores de linha (Pozar, cap. 6)"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("tipo", "Ressonador", rs.TIPOS, ""),
            ("f0", "Ressonância f0", 5.0, "GHz"),
            ("z0", "Z0 da linha", 50.0, "Ω"),
            ("acop", "Acoplamento por gap", ACOPLAMENTO, ""),
            ("c", "Capacitância do gap (se dada)", 0.05, "pF"),
            ("rvia", "Raio da via (λ/4 em curto)", 0.2, "mm"),
            ("faixa", "Faixa do gráfico (±)", 4.0, "%"),
        ], titulo="Ressonador em microstrip")
        self.form.pack(fill="x")
        self.corrigir = tk.BooleanVar(value=True)
        ttk.Checkbutton(self.painel, text="Descontar a ponta aberta (Kirschning-Jansen)",
                        variable=self.corrigir).pack(anchor="w", pady=(6, 0))
        self.botoes(("Calcular", self.calcular), ("Enviar p/ exportação", self.publicar),
                    ("Ex. 6.2: λ/2 em Teflon", self.ex_62))
        self.saida = AreaTexto(self.painel, 22)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Q0 = β/(2α) com as perdas do substrato atual (condutor e\n"
            "dielétrico, eq. 6.35). O gap (Garg-Bahl) acopla o ressonador à\n"
            "linha de alimentação; acoplamento crítico: bc = √(π/(2Q0)) para\n"
            "o λ/2 aberto (eq. 6.82) e √(π/(4Q0)) para o λ/4 em curto.")
        self._ultimo = None
        self.observar_previa(self.form)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    def ex_62(self):
        from ..core.microstrip import Substrato
        self.estado.definir_substrato(Substrato(2.08, 1.59, 35, 0.0004, "Teflon (Pozar Ex. 6.2)"))
        self.form.set(dict(tipo=rs.TIPOS[0], f0=5, z0=50, acop=ACOPLAMENTO[0]))
        self.corrigir.set(False)
        self.calcular()

    def _projeto(self):
        v = self.form.get()
        if v["f0"] <= 0 or v["z0"] <= 0:
            raise ValueError("f0 e Z0 devem ser positivos.")
        sub = self.estado.substrato
        r = rs.projeto(v["tipo"], v["f0"] * 1e9, v["z0"], sub, self.corrigir.get())
        c = None if v["acop"] == ACOPLAMENTO[0] else v["c"] * 1e-12
        if c is not None and c <= 0:
            raise ValueError("A capacitância do gap deve ser positiva.")
        a = rs.acoplamento_gap(r, sub, c)
        lay = None
        if a.s_gap is not None:
            info = {"tipo": f"ressonador - {v['tipo']}", "f0_GHz": v["f0"], "Q0": f"{r.q0:.1f}",
                    "C_gap_pF": f"{a.c * 1e12:.4f}", "gap_um": f"{a.s_gap * 1e3:.1f}"}
            lay = circuitos.ressonador_gap("ressonador_gap", sub, r.w, r.w, r.comp, a.s_gap,
                                           curto=v["tipo"] == rs.TIPOS[1], raio_via=v["rvia"], info=info)
        return v, sub, r, a, lay

    def atualizar_previa(self):
        try:
            *_, lay = self._projeto()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        if lay is None:
            self.previa.mostrar_erro("Gap fora do alcance de um gap simples: sem layout.")
        else:
            self.previa.mostrar(lay)

    @protegido
    def calcular(self):
        v, sub, r, a, lay = self._projeto()
        if lay is not None:
            self.previa.mostrar(lay)
        f0 = r.f0
        f = np.linspace(f0 * (1 - v["faixa"] / 100), f0 * (1 + v["faixa"] / 100), 4001)
        z = rs.zin_acoplado(f, r, a)
        g = (z - r.z0) / (z + r.z0)
        zr = r.zin(f)
        fg = f / 1e9
        ax = self.grafico.eixos(1, 2)
        ax[0].plot(fg, 20 * np.log10(np.maximum(np.abs(g), 1e-9)), label="|S11| com o gap")
        ax[0].axvline(a.f_carregada / 1e9, color="C1", ls="--", lw=0.8, label="ressonância carregada")
        ax[0].axvline(f0 / 1e9, color="gray", ls=":", lw=0.8, label="f0 sem gap")
        ax[0].set_xlabel("f [GHz]")
        ax[0].set_ylabel("|S11| [dB]")
        ax[0].set_title(f"{v['tipo']}: Q0 = {r.q0:.0f}, QL = {a.ql:.0f}", fontsize=9)
        ax[0].grid(True, alpha=0.3)
        ax[0].legend(fontsize=7)
        desenhar_smith(ax[1])
        ax[1].plot(g.real, g.imag, lw=1.2, label="Zin com o gap")
        gr = (zr - r.z0) / (zr + r.z0)
        ax[1].plot(gr.real, gr.imag, lw=0.8, ls="--", label="ressonador sozinho")
        ax[1].set_title("Carta de Smith", fontsize=9)
        ax[1].legend(fontsize=7, loc="lower left")
        self.grafico.desenhar()
        rr, ll, cc_ = r.rlc()
        i = int(np.argmin(np.abs(g)))
        linhas = r.texto + ["", f"RLC paralelo equivalente: R = {rr:.1f} Ω, L = {ll * 1e9:.4f} nH, "
                                f"C = {cc_ * 1e12:.4f} pF"] + [""] + a.texto
        linhas += [f"Mínimo de |S11| na varredura: {20 * np.log10(abs(g[i])):.1f} dB em {f[i] / 1e9:.4f} GHz"]
        linhas += ["", "Limites do substrato:"] + texto_limites(sub, r.w, f[-1])
        self.saida.escrever("\n".join(linhas))

    @protegido
    def publicar(self):
        *_, lay = self._projeto()
        if lay is None:
            raise ValueError("Gap não realizável: sem layout.")
        self.previa.mostrar(lay)
        self.estado.publicar_layout(lay)
        self.estado.msg(f"Layout '{lay.nome}' enviado para Exportar (Ctrl+E).")
