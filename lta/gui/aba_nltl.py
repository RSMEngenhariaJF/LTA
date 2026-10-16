"""Aba Linha não linear (NLTL): sólitons, decomposição de pulsos e multiplicação de frequência (seção 5.4)."""
import numpy as np

from ..core import nltl as nl
from ..core import sintonia as st
from .comum import AbaBase, AreaTexto, Formulario, protegido

MODOS = ("Sólitons: atraso × amplitude (eqs. 5.20-5.25)", "Decomposição de um pulso (Fig. 5.42b)",
         "Multiplicador de frequência (Figs. 5.43-5.44)")
DISP = ("HBV: CM·sech(kV)", "Varactor: Cj0/(1 + V/φ)^γ")


class AbaNLTL(AbaBase):
    TITULO = "Linha não linear (NLTL)"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("modo", "Análise", MODOS, ""),
            ("L", "L por seção", 108.7, "pH"),
            ("Co", "Co por seção", 43.5, "fF"),
            ("disp", "Dispositivo", DISP, ""),
            ("cm", "CM ou Cj0", 0.12, "pF"),
            ("k", "k (HBV)", 0.5, "1/V"),
            ("phi", "φ (varactor)", 0.7, "V"),
            ("gama", "γ (varactor)", 0.5, ""),
        ], titulo="Linha (Fig. 5.39b)")
        self.form.pack(fill="x")
        self.sim = Formulario(self.painel, [
            ("t1", "Sólitons: T de", 2.6, "ps"),
            ("t2", "até", 4.1, "ps"),
            ("amp", "Amplitude da excitação", 10.0, "V"),
            ("dur", "Pulso: duração", 40.0, "ps"),
            ("f0", "Senoide: f0", 20.0, "GHz"),
            ("nmin", "Seções: de", 9, ""),
            ("nmax", "até", 45, ""),
            ("passo", "passo", 3, ""),
            ("rf", "R fonte (0 = ideal; -1 = √(L/Co))", -1.0, "Ω"),
            ("rc", "R carga (-1 = √(L/Co))", -1.0, "Ω"),
            ("harm", "Harmônico", 3, ""),
        ], titulo="Simulação")
        self.sim.pack(fill="x", pady=(8, 0))
        self.botoes(("Calcular", self.calcular), ("Ex. sólitons (Figs. 5.40-5.41)", self.ex_soliton),
                    ("Ex. pulso 10 V, 40 ps", self.ex_pulso), ("Ex. triplicador 20 GHz", self.ex_mult))
        self.saida = AreaTexto(self.painel, 18)
        self.saida.pack(fill="both", expand=True)

    def _ex(self, modo):
        self.form.set(dict(modo=modo, L=108.7, Co=43.5, disp=DISP[0], cm=0.12, k=0.5))
        self.calcular()

    def ex_soliton(self):
        self.sim.set(dict(t1=2.6, t2=4.1))
        self._ex(MODOS[0])

    def ex_pulso(self):
        self.sim.set(dict(amp=10, dur=40, nmax=95))
        self._ex(MODOS[1])

    def ex_mult(self):
        self.sim.set(dict(amp=10, f0=20, nmin=9, nmax=45, passo=3, rf=-1, rc=-1, harm=3))
        self._ex(MODOS[2])

    def _linha(self):
        v = self.form.get()
        if v["L"] <= 0 or v["Co"] <= 0 or v["cm"] <= 0:
            raise ValueError("L, Co e CM/Cj0 devem ser positivos.")
        if v["disp"] == DISP[0]:
            cd = lambda x, a=v["cm"] * 1e-12, k=v["k"]: st.hbv(x, a, k)  # noqa: E731
        else:
            cd = lambda x, a=v["cm"] * 1e-12, p=v["phi"], g=v["gama"]: a / (1 + np.abs(x) / p) ** g  # noqa: E731
        return v, nl.LinhaNL(v["L"] * 1e-12, v["Co"] * 1e-15, cd)

    @protegido
    def calcular(self):
        v, lin = self._linha()
        s = self.sim.get()
        cab = [f"L = {lin.L * 1e12:g} pH, Co = {lin.Co * 1e15:g} fF, {v['disp']}",
               f"√(L/Co) = {lin.z0_linha:.1f} Ω; pequeno sinal em 0 V: {lin.z_pequeno_sinal():.1f} Ω, "
               f"Bragg {lin.f_bragg() / 1e9:.1f} GHz", ""]
        if v["modo"] == MODOS[0]:
            self._solitons(lin, s, cab)
        elif v["modo"] == MODOS[1]:
            self._pulso(lin, s, cab)
        else:
            self._multiplicador(lin, s, cab)

    def _r(self, lin, x):
        return lin.z0_linha if x < 0 else x

    def _solitons(self, lin, s, cab):
        tmin, tmax = nl.faixa_atraso(lin)
        ts = np.linspace(max(s["t1"] * 1e-12, tmin * 1.01), min(s["t2"] * 1e-12, tmax * 0.995), 40)
        amps = np.array([nl.amplitude_soliton(lin, t) for t in ts])
        ax = self.grafico.eixos(1, 3)
        vmax = amps.max() * 1.1
        vg, fv, _ = nl.funcoes_FG(lin, vmax)
        t_ex = ts[len(ts) // 3]
        ax[0].plot(vg, fv, label="F(V) (eq. 5.21)")
        ax[0].plot(vg, t_ex**2 * vg / lin.L, "--", label=f"T²V/L, T = {t_ex * 1e12:.2f} ps")
        ax[0].set_xlabel("V [V]")
        ax[0].set_ylabel("[C]")
        ax[0].set_title("Termos da eq. 5.20 (cf. Fig. 5.40)", fontsize=9)
        ax[1].plot(ts * 1e12, amps)
        ax[1].set_xlabel("T por seção [ps]")
        ax[1].set_ylabel("Vmax [V]")
        ax[1].set_title("Amplitude × atraso (eq. 5.24)", fontsize=9)
        linhas = cab + [f"Sólitons possíveis: T de {tmin * 1e12:.3f} a {tmax * 1e12:.3f} ps "
                        "(√(L·C(∞)) a √(L·C(0)))", "", " T[ps]   Vmax[V]   largura a meia altura [ps]"]
        for t in np.linspace(ts[0], ts[-1], 4):
            sx, vx = nl.perfil_soliton(lin, t)
            ax[2].plot(sx * 1e12, vx, label=f"T = {t * 1e12:.2f} ps")
            meia = vx > vx.max() / 2
            linhas.append(f" {t * 1e12:5.2f}   {vx.max():7.3f}   {(sx[meia].max() - sx[meia].min()) * 1e12:8.2f}")
        ax[2].set_xlabel("s = t - x/v [ps]")
        ax[2].set_ylabel("V [V]")
        ax[2].set_title("Forma do sóliton (cf. Fig. 5.41)", fontsize=9)
        for a in ax:
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        linhas += ["", "Sólitons maiores são mais rápidos (T menor) e mais",
                   "estreitos. A eq. 5.24 vem da aproximação contínua",
                   "(eq. 5.17): validada contra a simulação no tempo",
                   "com erro ≤ 1,5% nesta faixa."]
        self.saida.escrever("\n".join(linhas))

    def _pulso(self, lin, s, cab):
        n = max(int(s["nmax"]), 2)
        dur, amp = s["dur"] * 1e-12, s["amp"]
        rf, rc = self._r(lin, s["rf"]), self._r(lin, s["rc"])
        t0 = 5e-12
        sobe = 2e-12

        def pulso(t):
            u = np.clip((t - t0) / sobe, 0, 1) - np.clip((t - t0 - dur) / sobe, 0, 1)
            return amp * u * (2 if rf > 0 else 1)
        nos = sorted({0, n // 4, n // 2, 3 * n // 4, n - 1})
        t_fim = t0 + dur + 1.6 * n * nl.faixa_atraso(lin)[1] + 60e-12
        t, _, reg = nl.simular(lin, [n], pulso, t_fim, 0.025e-12, rf, rc, nos=nos)
        ax = self.grafico.eixos(1, 1)
        for k, no in enumerate(nos):
            ax[0].plot(t * 1e12, reg[:, k], label=f"seção {no + 1}")
        ax[0].set_xlabel("t [ps]")
        ax[0].set_ylabel("V [V]")
        ax[0].set_title(f"Pulso de {amp:g} V e {s['dur']:g} ps se decompondo em sólitons (cf. Fig. 5.42b)", fontsize=9)
        ax[0].grid(True, alpha=0.3)
        ax[0].legend(fontsize=7)
        self.grafico.desenhar()
        v_fim = reg[:, -1]
        picos = [i for i in range(1, len(v_fim) - 1) if v_fim[i] > v_fim[i - 1] and v_fim[i] >= v_fim[i + 1]
                 and v_fim[i] > 0.15 * amp]
        linhas = cab + [f"{n} seções; fonte {rf:.1f} Ω, carga {rc:.1f} Ω", "", "Picos na última seção:"]
        for i in picos[:8]:
            linhas.append(f"  t = {t[i] * 1e12:7.1f} ps   V = {v_fim[i]:6.2f} V   "
                          f"(T previsto pela eq. 5.24: {nl.atraso_soliton(lin, v_fim[i]) * 1e12:.3f} ps)")
        self.saida.escrever("\n".join(linhas + ["", "Os maiores sólitons chegam primeiro (mais rápidos)."]))

    def _multiplicador(self, lin, s, cab):
        ns = list(range(max(1, int(s["nmin"])), int(s["nmax"]) + 1, max(1, int(s["passo"]))))
        f0, amp, h = s["f0"] * 1e9, s["amp"], max(2, int(s["harm"]))
        rf, rc = self._r(lin, s["rf"]), self._r(lin, s["rc"])
        if h * f0 > lin.f_bragg(amp):
            aviso = f"⚠ {h}·f0 acima do corte em grande sinal ({lin.f_bragg(amp) / 1e9:.0f} GHz)."
        else:
            aviso = ""

        def fonte(t):
            return -amp * np.sin(2 * np.pi * f0 * t) * (2 if rf > 0 else 1)
        t_fim = 2 * max(ns) * nl.faixa_atraso(lin)[1] + 12 / f0
        t, sai, _ = nl.simular(lin, ns, fonte, t_fim, min(0.03e-12, 1 / (f0 * 600)), rf, rc)
        p_ref = amp**2 / (2 * (rf if rf > 0 else lin.z0_linha))
        ef = nl.eficiencia_harmonica(t, sai, f0, h, rc, p_ref)
        ef2 = nl.eficiencia_harmonica(t, sai, f0, 1, rc, p_ref)
        best = int(np.argmax(ef))
        ax = self.grafico.eixos(1, 2)
        ax[0].plot(ns, 100 * ef, "o-", label=f"{h}º harmônico")
        ax[0].plot(ns, 100 * ef2, "s--", label="fundamental")
        ax[0].set_xlabel("Nº de seções")
        ax[0].set_ylabel("P na carga / P de referência [%]")
        ax[0].set_title("Eficiência de conversão (sem perdas, cf. Fig. 5.44)", fontsize=9)
        sel = t >= t[-1] - 3 / f0
        ax[1].plot((t[sel] - t[sel][0]) * 1e12, sai[best, sel], label=f"{ns[best]} seções")
        ax[1].set_xlabel("t [ps]")
        ax[1].set_ylabel("V na carga [V]")
        ax[1].set_title("Forma de onda na saída (cf. Fig. 5.43)", fontsize=9)
        for a in ax:
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        self.saida.escrever("\n".join(cab + [
            f"Excitação: {amp:g} V em {s['f0']:g} GHz; fonte {rf:.1f} Ω, carga {rc:.1f} Ω",
            f"Máximo do {h}º harmônico: {ns[best]} seções, {100 * ef[best]:.1f}%", aviso, "",
            "Referência: potência disponível da fonte (ou A²/2Z0).",
            "Livro (Fig. 5.44, CM = 0,12 pF): máximo em 27 seções;",
            "CM maior → máximo com menos seções. O valor exato",
            "depende das terminações, que o livro não informa.",
            "Sem perdas: eficiências maiores que as reais."]))
