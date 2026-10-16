"""Aba Linhas sintonizáveis: dispositivos de sintonia, SRR com varactor e CRLH sintonizável (seções 5.2-5.3)."""
import numpy as np

from ..core import abcd
from ..core import crlh
from ..core import floquet as fl
from ..core import sintonia as st
from .comum import AbaBase, AreaTexto, Formulario, protegido

MODOS = ("Dispositivos: varactor, HBV e RF-MEMS", "SRR com varactor: rejeita-faixa (5.3.1.1)",
         "CRLH sintonizável: defasador e LWA (5.3.2)")


class AbaSintonia(AbaBase):
    TITULO = "Linhas sintonizáveis"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.modo = Formulario(self.painel, [("modo", "Ferramenta", MODOS, "")])
        self.modo.pack(fill="x")
        self.disp = Formulario(self.painel, [
            ("cj0", "Varactor Cj0", 9.0, "pF"),
            ("phi", "Varactor φ", 0.7, "V"),
            ("gama", "Varactor γ (1/2 abrupta)", 0.5, ""),
            ("vmax", "Tensão máxima", 30.0, "V"),
            ("cm", "HBV CM", 0.12, "pF"),
            ("k", "HBV k", 0.5, "1/V"),
            ("g0", "MEMS altura g0", 3.0, "µm"),
            ("kmola", "MEMS constante de mola", 10.0, "N/m"),
            ("area", "MEMS área", 10000.0, "µm²"),
        ], titulo="Dispositivos (5.2)")
        self.srr = Formulario(self.painel, [
            ("L", "L da linha (seção)", 3.3, "nH"),
            ("C", "C da linha (seção)", 1.3, "pF"),
            ("M", "Indutância mútua M", 1.2, "nH"),
            ("Lr", "Lr do anel", 3.0, "nH"),
            ("CL", "CL (borda)", 1.4, "pF"),
            ("Rs", "Rs (diodo + junções)", 4.8, "Ω"),
            ("cvmin", "Cvar mín", 1.0, "pF"),
            ("cvmax", "Cvar máx", 9.0, "pF"),
            ("n", "Nº de células", 2, ""),
        ], titulo="SRR com varactor (Fig. 5.7)")
        self.cr = Formulario(self.painel, [
            ("fcl", "f_CL (x = 1)", 1.0, "GHz"),
            ("fcr", "f_CR (x = 1)", 6.0, "GHz"),
            ("zb", "ZB", 50.0, "Ω"),
            ("xmin", "Fator x mín (CL, LL)", 0.6, ""),
            ("xmax", "Fator x máx", 1.6, ""),
            ("fop", "Frequência de operação", 2.45, "GHz"),
            ("n", "Nº de células", 6, ""),
            ("q", "Q dos elementos sintonizáveis", 40.0, ""),
            ("lcel", "Célula (LWA; 0 = ignorar)", 6.0, "mm"),
        ], titulo="CRLH balanceada sintonizável (eqs. 5.9-5.12)")
        self.estados = Formulario(self.painel, [("nest", "Nº de estados de sintonia", 5, ""),
                                                ("fmax", "f máx", 6.0, "GHz")])
        self.bot = self.botoes(("Calcular", self.calcular), ("Ex. Fig. 5.7 (VLSRR)", self.ex_vlsrr),
                               ("Ex. defasador/LWA", self.ex_crlh))
        self.saida = AreaTexto(self.painel, 18)
        self.modo.ao_mudar("modo", self._organizar)
        self._organizar()

    def _organizar(self):
        for w in (self.disp, self.srr, self.cr, self.estados, self.bot, self.saida):
            w.pack_forget()
        m = self.modo.get()["modo"]
        {MODOS[0]: self.disp, MODOS[1]: self.srr, MODOS[2]: self.cr}[m].pack(fill="x", pady=(8, 0))
        if m != MODOS[0]:
            self.estados.pack(fill="x", pady=(8, 0))
        self.bot.pack(fill="x", pady=4)
        self.saida.pack(fill="both", expand=True)

    def ex_vlsrr(self):
        """Livro, Fig. 5.7: Lr = 3 nH, CL = 1,4 pF, L = 3,3 nH, C = 1,3 pF, M = 1,2 nH, Rs = 4,8 Ω, 1-9 pF."""
        self.modo.set({"modo": MODOS[1]})
        self.srr.set(dict(L=3.3, C=1.3, M=1.2, Lr=3.0, CL=1.4, Rs=4.8, cvmin=1, cvmax=9, n=2))
        self.estados.set(dict(nest=5, fmax=5))
        self._organizar()
        self.calcular()

    def ex_crlh(self):
        self.modo.set({"modo": MODOS[2]})
        self.cr.set(dict(fcl=1, fcr=6, zb=50, xmin=0.6, xmax=1.6, fop=2.45, n=6, q=40, lcel=6))
        self.estados.set(dict(nest=5, fmax=6))
        self._organizar()
        self.calcular()

    @protegido
    def calcular(self):
        m = self.modo.get()["modo"]
        if m == MODOS[0]:
            return self._dispositivos()
        e = self.estados.get()
        if e["nest"] < 2:
            raise ValueError("Use ao menos 2 estados de sintonia.")
        f = np.linspace(e["fmax"] * 1e9 / 3000, e["fmax"] * 1e9, 3000)
        if m == MODOS[1]:
            self._vlsrr(f, e["nest"])
        else:
            self._crlh(f, e["nest"])

    # ------------------------------------------------------------------ dispositivos
    def _dispositivos(self):
        v = self.disp.get()
        tens = np.linspace(0, v["vmax"], 400)
        cvar = st.varactor(tens, v["cj0"] * 1e-12, v["phi"], v["gama"])
        vh = np.linspace(-v["vmax"] / 3, v["vmax"] / 3, 400)
        ch = st.hbv(vh, v["cm"] * 1e-12, v["k"])
        g0, area = v["g0"] * 1e-6, v["area"] * 1e-12
        g = np.linspace(0.999, 0.6667, 300) * g0
        vg = st.mems_tensao(g, g0, v["kmola"], area)
        vp = st.mems_pull_in(g0, v["kmola"], area)
        ax = self.grafico.eixos(1, 3)
        ax[0].plot(tens, cvar * 1e12)
        ax[0].set_xlabel("Tensão reversa [V]")
        ax[0].set_ylabel("C [pF]")
        ax[0].set_title(f"Varactor: Cj0/(1 + V/φ)^{v['gama']:g}", fontsize=9)
        ax[1].plot(vh, ch * 1e12)
        ax[1].set_xlabel("V [V]")
        ax[1].set_title("HBV: CM·sech(kV) (simétrico)", fontsize=9)
        ax[2].plot(vg, g / g0, label="estável (eq. 5.4)")
        ax[2].axvline(vp, color="C3", ls="--", label=f"pull-in Vp = {vp:.2f} V")
        ax[2].set_xlabel("V [V]")
        ax[2].set_ylabel("g/g0")
        ax[2].set_title("RF-MEMS: altura da ponte", fontsize=9)
        ax[2].legend(fontsize=7)
        for a in ax:
            a.grid(True, alpha=0.3)
        self.grafico.desenhar()
        c_up = st.mems_capacitancia(g0, area)
        self.saida.escrever("\n".join([
            f"Varactor: {cvar[0] * 1e12:.3g} pF (0 V) a {cvar[-1] * 1e12:.3g} pF ({v['vmax']:g} V), "
            f"razão {cvar[0] / cvar[-1]:.1f}:1",
            f"HBV: máximo {ch.max() * 1e12:.3g} pF em 0 V (C(V) par: não precisa de polarização)", "",
            "RF-MEMS (placas paralelas, eqs. 5.1-5.5):",
            f"  C(g0) = {c_up * 1e15:.3g} fF; em g = 2g0/3: {c_up * 1.5 * 1e15:.3g} fF",
            f"  Vp = √(8k·g0³/(27ε0·A)) = {vp:.3f} V (instável abaixo de 2g0/3: colapso)",
            "  No estado baixo a capacitância é limitada pelo isolante",
            "  (desprezado aqui), o que dá as razões altas típicas."]))

    # ------------------------------------------------------------------ VLSRR
    def _vlsrr(self, f, nest):
        v = self.srr.get()
        cvars = np.linspace(v["cvmin"], v["cvmax"], nest) * 1e-12
        ax = self.grafico.eixos(1, 2)
        linhas = ["SRR com varactor (eqs. 5.6): Ceq = CL em série com Cvar", "",
                  " Cvar[pF]  Ceq[pF]  f0[GHz]   Q    notch[GHz]  |S21|mín[dB]  BW-10dB[MHz]"]
        f0s, qs = [], []
        for cv in cvars:
            ceq = st.vlsrr_ceq(v["CL"] * 1e-12, cv)
            ls, cs, r, w0 = st.vlsrr_tanque(v["M"] * 1e-9, v["Lr"] * 1e-9, ceq, v["Rs"])
            mm = st.vlsrr_abcd(f, v["L"] * 1e-9, v["C"] * 1e-12, v["M"] * 1e-9, v["Lr"] * 1e-9, ceq, v["Rs"])
            s = abcd.abcd_para_s(abcd.potencia(mm, v["n"]), 50)
            s21 = abcd.db(s[:, 1, 0])
            i = int(np.argmin(s21))
            rej = np.flatnonzero(s21 < -10)
            bw = (f[rej[-1]] - f[rej[0]]) / 1e6 if rej.size else 0.0
            q = w0 * r * cs
            f0s.append(w0 / 2 / np.pi)
            qs.append(q)
            ax[0].plot(f / 1e9, s21, label=f"Cvar = {cv * 1e12:.2g} pF")
            linhas.append(f" {cv * 1e12:7.2f}  {ceq * 1e12:7.3f}  {w0 / 2e9 / np.pi:7.3f}  {q:5.1f}  "
                          f"{f[i] / 1e9:9.3f}  {s21[i]:11.1f}  {bw:11.0f}")
        ax[0].set_ylim(-60, 1)
        ax[0].set_xlabel("f [GHz]")
        ax[0].set_ylabel("|S21| [dB]")
        ax[0].set_title(f"{v['n']} células com VLSRR", fontsize=9)
        ax[0].legend(fontsize=7)
        ax[1].plot(cvars * 1e12, np.array(f0s) / 1e9, "o-", label="f0")
        ax2 = ax[1].twinx()
        ax2.plot(cvars * 1e12, qs, "s--", color="C1", label="Q = ω0·Lr/Rs")
        ax2.set_ylabel("Q", color="C1")
        ax[1].set_xlabel("Cvar [pF]")
        ax[1].set_ylabel("f0 [GHz]")
        ax[1].set_title("Sintonia e fator Q do tanque", fontsize=9)
        for a in ax:
            a.grid(True, alpha=0.3)
        self.grafico.desenhar()
        linhas += ["", f"Faixa de sintonia: {min(f0s) / 1e9:.3f} a {max(f0s) / 1e9:.3f} GHz",
                   "Q = ω0·Lr/Rs cresce com f0: rejeição mais funda e", "estreita no topo da faixa (livro, p. 373).",
                   "Ls' = 4M²/Lr não depende de Cvar."]
        self.saida.escrever("\n".join(linhas))

    # ------------------------------------------------------------------ CRLH sintonizável
    def _crlh(self, f, nest):
        v = self.cr.get()
        base = crlh.balanceada_por_banda(v["fcl"] * 1e9, v["fcr"] * 1e9, v["zb"])
        xs = np.linspace(v["xmin"], v["xmax"], nest)
        fop = v["fop"] * 1e9
        q = v["q"] if v["q"] > 0 else np.inf
        ax = self.grafico.eixos(2, 2)
        fases, perdas, linhas_tab = [], [], []
        lcel = v["lcel"] * 1e-3
        for x in xs:
            c = st.crlh_sintonizada(base, x)
            m = abcd.potencia(st.crlh_abcd_q(c, f, q), v["n"])
            s = abcd.abcd_para_s(m, v["zb"])
            bl, al = c.beta_l(f)
            zb = c.zb(f)
            ax[0].plot(bl / np.pi, f / 1e9, label=f"x = {x:.2f}")
            ax[1].plot(f / 1e9, abcd.db(s[:, 1, 0]), label=f"x = {x:.2f}")
            ax[2].plot(f / 1e9, zb.real, label=f"x = {x:.2f}")
            k = np.argmin(abs(f - fop))
            fase = np.degrees(np.angle(s[k, 1, 0]))
            fases.append(fase)
            perdas.append(abcd.db(s[k, 1, 0]))
            th = ""
            if lcel > 0:
                h = fl.harmonicos(np.array([fop]), np.array([bl[k] if al[k] < 1e-6 else np.nan]), lcel, 0)[0]
                th = f"{h.theta_graus[0]:6.1f}" if h.rapido[0] else "   —  "
            linhas_tab.append(f" {x:5.2f}  {crlh.frequencia_transicao(c) / 1e9:7.3f}  {np.degrees(bl[k]):8.1f}  "
                              f"{fase:8.1f}  {abcd.db(s[k, 1, 0]):7.2f}  {th}")
        fases = np.degrees(np.unwrap(np.radians(fases)))
        ax[3].plot(xs, fases - fases[0], "o-", label="Δ∠S21 em f_op")
        ax[3].set_xlabel("Fator de sintonia x")
        ax[3].set_ylabel("Defasagem [°]")
        ax[3].set_title(f"Defasador em {v['fop']:g} GHz", fontsize=9)
        ax[0].axhline(v["fop"], color="gray", ls=":", lw=0.8)
        ax[0].set_xlim(-1.05, 1.05)
        ax[0].set_xlabel("βl/π (LH < 0)")
        ax[0].set_ylabel("f [GHz]")
        ax[0].set_title("Dispersão por estado", fontsize=9)
        ax[1].set_ylim(-30, 1)
        ax[1].set_title(f"|S21|, {v['n']} células, Q = {v['q']:g}", fontsize=9)
        ax[2].set_ylim(0, 2 * v["zb"])
        ax[2].set_title("Re ZB: constante na transição (balanceada)", fontsize=9)
        for a in ax[1:3]:
            a.set_xlabel("f [GHz]")
        for a in ax:
            a.grid(True, alpha=0.3)
            a.legend(fontsize=6)
        self.grafico.desenhar()
        fom = st.figura_merito(fases, perdas)
        self.saida.escrever("\n".join([
            f"CRLH balanceada (x = 1): LR = {base.LR * 1e9:.3g} nH, CR = {base.CR * 1e12:.3g} pF, "
            f"LL = {base.LL * 1e9:.3g} nH, CL = {base.CL * 1e12:.3g} pF",
            "Sintonia: CL e LL multiplicados por x (eqs. 5.9-5.11): ZB = √(LR/CR)",
            "permanece igual e a transição vai para f0/√x.", "",
            f" x      f0[GHz]   βl[°]    ∠S21[°]  S21[dB]  θ LWA[°] (em {v['fop']:g} GHz)"] + linhas_tab + [
            "", f"Defasagem total: {fases.max() - fases.min():.1f}°   perda máx: {-min(perdas):.2f} dB",
            f"FoM (eq. 5.12) = {fom:.1f} °/dB" + (" (Q infinito: sem perdas)" if not np.isfinite(fom) else ""),
            "", "θ LWA: feixe do modo fundamental em frequência fixa",
            "(varredura por sintonia, livro, seção 5.3.2.2)."]))
