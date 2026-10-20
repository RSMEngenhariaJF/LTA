"""Aba Amplificador: estabilidade, ganhos e casamento conjugado (Steer, seções 11.5-11.7)."""
import tkinter as tk
from tkinter import filedialog, ttk

import numpy as np

from ..core import amplificador as am
from ..core import casamento as cs
from ..exportacao.simuladores import ler_touchstone
from .comum import PASTA_SAIDAS, AbaBase, AreaTexto, Formulario, desenhar_smith, protegido

FONTES = ("Exemplo Steer: pHEMT FPD6836P70 (Fig. 11-8)", "Arquivo .s2p do transistor")
REDES = ("Rede L (elementos concentrados)", "Stub simples em derivação (aberto)", "Stub simples em derivação (curto)")


def _db(x):
    return 10 * np.log10(np.maximum(x, 1e-30))


class AbaAmplificador(AbaBase):
    TITULO = "Amplificador: estabilidade, ganho e casamento"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("fonte", "Parâmetros S", FONTES, ""),
            ("f0", "Frequência de projeto f0", 8.0, "GHz"),
            ("z0", "Z0 do sistema", 50.0, "Ω"),
            ("rede", "Redes de casamento", REDES, ""),
            ("sol", "Solução da rede (1 ou 2)", 1, ""),
        ], titulo="Transistor")
        self.form.pack(fill="x")
        self.arq = tk.StringVar(value="")
        f = ttk.Frame(self.painel)
        f.pack(fill="x", pady=1)
        ttk.Button(f, text="Arquivo .s2p", command=self._escolher).pack(side="left")
        ttk.Label(f, textvariable=self.arq, foreground="#555").pack(side="left", padx=4)
        self.botoes(("Calcular", self.calcular), ("Steer Tab. 11-7: 8 GHz", self.ex_8ghz),
                    ("Instável: 2 GHz", self.ex_2ghz))
        self.saida = AreaTexto(self.painel, 24)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Steer, cap. 11:\n"
            "• estabilidade incondicional: k > 1 e |Δ| < 1, ou μ > 1\n"
            "  (Edwards-Sinsky, um só teste);\n"
            "• círculos de estabilidade de fonte e de carga na Smith\n"
            "  (região instável é a que não contém o centro se |S11| < 1);\n"
            "• GTU,max (unilateral), GMA (k > 1), GMS = |S21/S12| (k < 1)\n"
            "  e U de Mason;\n"
            "• casamento conjugado simultâneo ΓS, ΓL → ZS, ZL e redes de\n"
            "  casamento (rede L ou stub) para 50 Ω, com o GT resultante.")

    def _escolher(self):
        c = filedialog.askopenfilename(initialdir=PASTA_SAIDAS, filetypes=[("Touchstone 2 portas", "*.s2p")])
        if c:
            self.arq.set(c)
            self.form.set(dict(fonte=FONTES[1]))

    def ex_8ghz(self):
        self.form.set(dict(fonte=FONTES[0], f0=8, z0=50, rede=REDES[0], sol=1))
        self.calcular()

    def ex_2ghz(self):
        self.form.set(dict(fonte=FONTES[0], f0=2, z0=50, rede=REDES[0], sol=1))
        self.calcular()

    def _dados(self, v):
        if v["fonte"] == FONTES[0]:
            f, s = am.s_de_tabela(am.PHEMT_STEER)
            return f, s, 50.0, "pHEMT FPD6836P70 (Steer, Fig. 11-8)"
        if not self.arq.get():
            raise ValueError("Escolha o arquivo .s2p do transistor.")
        f, s, zr = ler_touchstone(self.arq.get())
        return f, s, zr, self.arq.get()

    def _rede(self, v, z_alvo, z0, f0):
        """(função f -> ABCD da fonte z0 até o transistor, texto) que apresenta z_alvo ao transistor."""
        zc = np.conj(z_alvo)                # casa ZS* (ou ZL*) com z0; olhando de volta vê-se z_alvo
        k = max(v["sol"], 1) - 1
        if v["rede"] == REDES[0]:
            sols = cs.rede_l(zc, z0, f0)
            s = sols[min(k, len(sols) - 1)]
            return s.abcd, f"rede L ({s.topologia}): série {cs.fmt_comp(s.serie)}, derivação {cs.fmt_comp(s.derivacao)}"
        curto = v["rede"] == REDES[2]
        sols = cs.stub_simples(zc, z0, serie=False, curto=curto)
        s = sols[min(k, len(sols) - 1)]
        return (lambda f: cs.abcd_stub_simples(f, f0, z0, s, curto, False),
                f"stub {'em curto' if curto else 'aberto'}: d = {s.d:.4f}λ da carga, ℓ = {s.l:.4f}λ (Z0 = {z0:g} Ω)")

    @protegido
    def calcular(self):
        v = self.form.get()
        f, s, zr, nome = self._dados(v)
        z0 = v["z0"]
        if abs(zr - z0) > 1e-9:
            from ..core import caracterizacao as ca
            s = ca.renormalizar(s, zr, z0)
        f0 = v["f0"] * 1e9
        if not f[0] <= f0 <= f[-1]:
            raise ValueError(f"f0 fora da faixa dos dados ({f[0] / 1e9:g}-{f[-1] / 1e9:g} GHz).")
        est = am.estabilidade(s)
        gan = am.ganhos(s)
        s0 = am.interpolar(f, s, f0)[None]
        e0 = am.estabilidade(s0)
        g0 = am.ganhos(s0)
        (cs_, rs), (cl, rl) = am.circulos(s0)
        fg = f / 1e9
        ax = self.grafico.eixos(2, 2)
        ax[0].plot(fg, est.k, "o-", ms=3, label="k (Rollett)")
        ax[0].plot(fg, est.mu, "s-", ms=3, label="μ (entrada)")
        ax[0].plot(fg, est.delta, "^-", ms=3, label="|Δ|")
        ax[0].axhline(1, color="gray", ls=":", lw=0.8)
        ax[0].axvline(v["f0"], color="gray", ls=":", lw=0.8)
        ax[0].set_ylim(0, max(3, float(np.nanmin([np.nanmax(est.k), 6]))))
        ax[0].set_title("Estabilidade", fontsize=9)
        ax[1].plot(fg, _db(gan["GTU_max"]), "o-", ms=3, label="GTU,max")
        ax[1].plot(fg, _db(gan["GMA"]), "s-", ms=3, label="GMA (k > 1)")
        ax[1].plot(fg, np.where(est.k < 1, _db(gan["GMS"]), np.nan), "v-", ms=3, label="GMS (k < 1)")
        ax[1].plot(fg, _db(gan["U"]), "--", lw=0.8, label="U de Mason")
        ax[1].plot(fg, _db(gan["S21"]), lw=0.8, label="|S21|²")
        ax[1].axvline(v["f0"], color="gray", ls=":", lw=0.8)
        ax[1].set_ylabel("[dB]")
        ax[1].set_title("Ganhos", fontsize=9)
        desenhar_smith(ax[2])
        t = np.linspace(0, 2 * np.pi, 361)
        ax[2].plot((cs_ + rs * np.exp(1j * t)).real, (cs_ + rs * np.exp(1j * t)).imag, color="C0",
                   label="círculo de fonte")
        ax[2].plot((cl + rl * np.exp(1j * t)).real, (cl + rl * np.exp(1j * t)).imag, color="C3",
                   label="círculo de carga")
        ax[2].set_title(f"Círculos de estabilidade em {v['f0']:g} GHz", fontsize=9)
        linhas = [f"Transistor: {nome}", f"Z0 = {z0:g} Ω; dados de {fg[0]:g} a {fg[-1]:g} GHz ({len(f)} pontos)", ""]
        inc = est.incondicional
        faixas = [fg[i] for i in range(len(f)) if inc[i]]
        linhas.append("Incondicionalmente estável (k > 1 e |Δ| < 1) em: "
                      + (", ".join(f"{x:g}" for x in faixas) + " GHz" if faixas else "nenhuma frequência"))
        linhas += ["", f"Em f0 = {v['f0']:g} GHz (interpolado): k = {e0.k[0]:.4f}, |Δ| = {e0.delta[0]:.4f}, "
                       f"μ = {e0.mu[0]:.4f}",
                   f"  GTU,max = {_db(g0['GTU_max'][0]):.2f} dB, GMS = {_db(g0['GMS'][0]):.2f} dB, "
                   f"U = {_db(g0['U'][0]):.2f} dB",
                   f"  Círculo de fonte: centro {cs_[0].real:.4f}{cs_[0].imag:+.4f}j, raio {rs[0]:.4f}",
                   f"  Círculo de carga: centro {cl[0].real:.4f}{cl[0].imag:+.4f}j, raio {rl[0]:.4f}"]
        ax[3].axvline(v["f0"], color="gray", ls=":", lw=0.8)
        if e0.incondicional[0]:
            gs, gl = (x[0] for x in am.conjugado_simultaneo(s0))
            zs, zl = am.z_de_gama(gs, z0), am.z_de_gama(gl, z0)
            ax[2].plot(gs.real, gs.imag, "o", color="C0", label="ΓS")
            ax[2].plot(gl.real, gl.imag, "s", color="C3", label="ΓL")
            m_in, t_in = self._rede(v, zs, z0, f0)
            m_out, t_out = self._rede(v, zl, z0, f0)
            fd = np.linspace(f[0], f[-1], 600)
            sd = np.array([am.interpolar(f, s, x) for x in fd])
            ma, mb = m_in(fd), m_out(fd)
            # impedância vista pelo transistor olhando para a fonte (e para a carga) através das redes
            z_s = (ma[:, 1, 1] * z0 + ma[:, 0, 1]) / (ma[:, 1, 0] * z0 + ma[:, 0, 0])
            z_l = (mb[:, 1, 1] * z0 + mb[:, 0, 1]) / (mb[:, 1, 0] * z0 + mb[:, 0, 0])
            gsf, glf = cs.gama(z_s, z0), cs.gama(z_l, z0)
            gtf = am.gt(sd, gsf, glf)
            ax[3].plot(fd / 1e9, _db(gtf), label="GT com as redes")
            ax[3].plot(fg, _db(gan["GMA"]), "s", ms=3, label="GMA")
            ax[3].set_ylim(max(-20, float(np.nanmin(_db(gtf)))), float(np.nanmax(_db(gan["GTU_max"]))) + 3)
            gt0 = float(am.gt(s0[0], gs, gl))
            linhas += ["", "Casamento conjugado simultâneo:",
                       f"  ΓS = {abs(gs):.4f} ∠{np.degrees(np.angle(gs)):.2f}° → ZS = {zs.real:.3f} {zs.imag:+.3f}j Ω",
                       f"  ΓL = {abs(gl):.4f} ∠{np.degrees(np.angle(gl)):.2f}° → ZL = {zl.real:.3f} {zl.imag:+.3f}j Ω",
                       f"  GT = GMA = {10 * np.log10(gt0):.2f} dB",
                       f"  Entrada ({z0:g} Ω → ZS): {t_in}",
                       f"  Saída (ZL ← {z0:g} Ω): {t_out}",
                       "  (as redes são vistas do lado de 50 Ω; a ordem é fonte → rede → transistor → rede → carga)",
                       "", "Polarização (bias tees, capacitores de bloqueio) e perdas das redes não incluídas."]
        else:
            ax[3].text(0.5, 0.5, "Sem casamento conjugado:\no transistor não é\nincondicionalmente estável em f0",
                       ha="center", va="center", transform=ax[3].transAxes, fontsize=9)
            linhas += ["", "Não é incondicionalmente estável em f0: não há casamento conjugado simultâneo.",
                       f"Use o GMS = {_db(g0['GMS'][0]):.2f} dB como teto, escolha ΓS e ΓL fora das regiões instáveis",
                       "(círculos) ou estabilize com resistor série/derivação na entrada ou na saída (Steer 11.6)."]
        ax[3].set_ylabel("[dB]")
        ax[3].set_title("Ganho de transdução projetado", fontsize=9)
        for a in ax:
            a.legend(fontsize=7)
            if a is not ax[2]:
                a.set_xlabel("f [GHz]")
                a.grid(True, alpha=0.3)
        self.grafico.desenhar()
        self.saida.escrever("\n".join(linhas))
