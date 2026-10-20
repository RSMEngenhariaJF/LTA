"""Aba Balun de Marchand com linhas acopladas (Garg, Bahl e Bianchi, seção 4.9), com versão LTA e
versão miniaturizada com capacitores (Steer, seção 13.2)."""
import tkinter as tk
from tkinter import ttk

import numpy as np

from ..core import balun
from ..core import linhas_acopladas as la
from ..core.microstrip import sintese
from ..geometria import circuitos
from .comum import AbaBase, AreaTexto, Formulario, protegido, salvar_csv_s
from .previa import PreviaEstrutura

IMPL = ("Convencional", "LTA: modo ímpar com carga capacitiva periódica",
        "Miniaturizado: seções θ0 < 90° + capacitores ao terra (Steer 13.2)")


class AbaBalun(AbaBase):
    TITULO = "Balun de Marchand"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("zs", "Fonte desbalanceada Zs", 50.0, "Ω"),
            ("zl", "Carga balanceada ZL (diferencial)", 200.0, "Ω"),
            ("zlin", "Z0 = √(Z0e·Z0o) das linhas (0 = Zs)", 0.0, "Ω"),
            ("f0", "Frequência central f0", 2.0, "GHz"),
            ("impl", "Implementação", IMPL, ""),
            ("ncel", "LTA: capacitores por seção", 6, ""),
            ("theta0", "Miniaturizado: θ0 das seções em f0", 45.0, "°"),
            ("modelo", "Modelo das linhas acopladas", la.MODELOS, ""),
            ("rvia", "Raio das vias de terra", 0.2, "mm"),
            ("fmax", "f máx", 4.0, "GHz"),
        ], titulo="Balun de Marchand (2 seções λ/4)")
        self.form.pack(fill="x")
        self.perdas = tk.BooleanVar(value=False)
        self.dispersao = tk.BooleanVar(value=False)
        ttk.Checkbutton(self.painel, text="Perdas por modo", variable=self.perdas).pack(anchor="w", pady=(6, 0))
        ttk.Checkbutton(self.painel, text="Dispersão (Kirschning-Jansen)", variable=self.dispersao).pack(anchor="w")
        self.botoes(("Calcular", self.calcular), ("Enviar p/ exportação", self.publicar),
                    ("Exportar S (.csv)", self.exportar))
        self.saida = AreaTexto(self.painel, 22)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Portas: 1 = entrada desbalanceada (Zs); 2 e 3 = saídas balanceadas,\n"
            "cada uma com ZL/2 ao terra (carga diferencial ZL).\n"
            "O acoplamento é escolhido para casar a entrada em f0 com linhas\n"
            "TEM ideais (com Z0 = Zs, C = 1/√(1 + ZL/Zs)) e depois W e S são\n"
            "sintetizados em microstrip. Acoplamentos fortes (ZL baixo) não\n"
            "cabem em linhas lado a lado: use ZL maior, outro Z0 das linhas,\n"
            "ou linhas pela face (broadside) / interdigitadas.\n"
            "Em microstrip, εeff par ≠ ímpar desequilibra a fase; a opção LTA\n"
            "iguala as velocidades com capacitores entre as linhas.\n"
            "A opção miniaturizada encurta as seções para θ0 e casa com\n"
            "capacitores ao terra no meio, na ponta aberta e nas saídas.")
        self._ultimo = None
        self._cache = (None, None)
        self.observar_previa(self.form)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    def _projeto(self):
        v = self.form.get()
        sub = self.estado.substrato
        lta = v["impl"] == IMPL[1]
        mini = v["impl"] == IMPL[2]
        if mini:
            chave = (v["f0"], v["zs"], v["zl"], v["zlin"], v["theta0"], v["modelo"], repr(sub))
            if self._cache[0] == chave:
                p = self._cache[1]
            else:
                p = balun.projeto_miniaturizado(v["f0"] * 1e9, v["zs"], v["zl"], v["theta0"], v["zlin"] or None,
                                                sub, v["modelo"])
                self._cache = (chave, p)
        else:
            p = balun.projeto(v["f0"] * 1e9, v["zs"], v["zl"], sub, v["zlin"] or None, v["modelo"],
                              self.dispersao.get(), lta, max(1, v["ncel"]))
        sec = p.secoes[0]
        info = {"tipo": "balun de Marchand", "C_dB": f"{-20 * np.log10(p.c):.2f}", "W_S_l_mm":
                f"{sec.w:.4f} / {sec.s:.4f} / {sec.comp:.3f}"}
        if sec.n_cel:
            info["Cc_fF"] = f"{sec.cc * 1e15:.2f} x {sec.n_cel} por secao"
        w0 = float(sintese(v["zs"], sub))
        w_out = float(sintese(v["zl"] / 2, sub))
        lay = circuitos.balun_marchand("balun_marchand" + ("_lta" if lta else "_mini" if mini else ""), sub, w0,
                                       sec.w, sec.s, sec.comp, w_out, v["rvia"], info=info, chanfro=True,
                                       capacitores=p.caps if mini else None)
        return v, sub, p, lay

    def atualizar_previa(self):
        try:
            *_, lay = self._projeto()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        self.previa.mostrar(lay)

    @protegido
    def calcular(self):
        v, sub, p, lay = self._projeto()
        self.previa.mostrar(lay)
        f = np.linspace(v["fmax"] * 1e9 / 800, v["fmax"] * 1e9, 800)
        caps = getattr(p, "caps", None)
        s = balun.s_marchand(f, *p.secoes, v["zs"], v["zl"], sub, self.perdas.get(), caps=caps)
        if caps is not None:
            s_id = balun.s_marchand(f, *p.secoes_ideais, v["zs"], v["zl"], caps=p.caps_ideais)
        else:
            ideais = balun._secoes_ideais(p.c, v["zlin"] or v["zs"], v["f0"] * 1e9)
            s_id = balun.s_marchand(f, *ideais, v["zs"], v["zl"])
        self._ultimo = (f, s)
        m, mi = balun.metricas(s), balun.metricas(s_id)
        fg = f / 1e9
        ax = self.grafico.eixos(2, 2)
        ax[0].plot(fg, m[0], label="|S11|")
        ax[0].plot(fg, m[1], label="|S21|")
        ax[0].plot(fg, m[2], label="|S31|")
        ax[0].plot(fg, mi[0], ls="--", lw=0.8, label="|S11| TEM ideal")
        ax[0].set_ylim(-40, 1)
        ax[0].set_ylabel("|S| [dB]")
        ax[0].set_title("Casamento e divisão", fontsize=9)
        ax[1].plot(fg, m[3], label="microstrip")
        ax[1].plot(fg, mi[3], ls="--", label="TEM ideal")
        ax[1].set_ylim(-3, 3)
        ax[1].set_ylabel("|S21| - |S31| [dB]")
        ax[1].set_title("Desequilíbrio de amplitude", fontsize=9)
        ax[2].plot(fg, np.abs(m[4]), label="microstrip")
        ax[2].plot(fg, np.abs(mi[4]), ls="--", label="TEM ideal")
        ax[2].set_ylim(150, 190)
        ax[2].set_ylabel("|∠S21 - ∠S31| [°]")
        ax[2].set_title("Fase entre as saídas (ideal: 180°)", fontsize=9)
        ax[3].plot(fg, 10 ** (m[1] / 10) + 10 ** (m[2] / 10), label="|S21|² + |S31|²")
        ax[3].set_ylim(0, 1.05)
        ax[3].set_title("Potência nas saídas", fontsize=9)
        for a in ax:
            a.axvline(v["f0"], color="gray", ls=":", lw=0.8)
            a.set_xlabel("f [GHz]")
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        i = int(np.argmin(np.abs(f - v["f0"] * 1e9)))
        ok = (m[0] < -15) & (np.abs(m[3]) < 1) & (np.abs(np.abs(m[4]) - 180) < 10)
        a = b = i
        if ok[i]:
            while a > 0 and ok[a - 1]:
                a -= 1
            while b < len(f) - 1 and ok[b + 1]:
                b += 1
        linhas = p.texto + ["", f"Em f0: |S11| = {m[0][i]:.1f} dB, |S21| = {m[1][i]:.2f} dB, |S31| = {m[2][i]:.2f} dB, "
                                f"fase = {m[4][i]:.1f}°",
                            "Banda com |S11| < -15 dB, desequilíbrio < 1 dB e fase 180° ± 10°: "
                            + (f"{f[a] / 1e9:.3f}-{f[b] / 1e9:.3f} GHz ({100 * (f[b] - f[a]) / f[i]:.0f}%)" if ok[i]
                               else "não atende em f0"),
                            "", "Layout: vias nas pontas externas das linhas b; saídas de ZL/2 com curvas chanfradas.",
                            "Junções, vias e a separação entre as duas seções não estão no modelo."]
        if caps is not None:
            conv = balun.projeto(v["f0"] * 1e9, v["zs"], v["zl"], None, v["zlin"] or None)
            l90 = la.comprimento_quarto(p.secoes[0].modos, v["f0"] * 1e9)
            linhas += ["", f"Miniaturizado: seções de {v['theta0']:g}° ({p.secoes[0].comp:.2f} mm em vez de "
                           f"{l90:.2f} mm); acoplamento {-20 * np.log10(p.c):.2f} dB (convencional: "
                           f"{-20 * np.log10(conv.c):.2f} dB).",
                       "Capacitores SMD ao terra (via curta) ou stubs abertos curtos; a banda diminui com θ0",
                       "e o 1º espúrio sobe para perto de 2·fr = 180°/θ0·f0."]
        self.saida.escrever("\n".join(linhas))

    @protegido
    def publicar(self):
        *_, lay = self._projeto()
        self.previa.mostrar(lay)
        self.estado.publicar_layout(lay)
        self.estado.msg(f"Layout '{lay.nome}' enviado para Exportar (Ctrl+E).")

    def exportar(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            salvar_csv_s(self, "balun", *self._ultimo)
