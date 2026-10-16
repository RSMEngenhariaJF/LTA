"""Aba Linhas dual-band e multibanda: projeto por engenharia de dispersão e impedância (seção 4.2.2)."""
import numpy as np

from ..core import abcd
from ..core import multibanda as mb
from ..core.microstrip import C0, analise, sintese
from .comum import AbaBase, AreaTexto, Formulario, protegido

MODOS = ("Linha convencional + PLH (4.20-4.27)", "CRLH canônica (4.28-4.32)", "CSRR + gap (4.33-4.34)",
         "Quad-band E-CRLH (4.37)")


class AbaMultibanda(AbaBase):
    TITULO = "Linhas dual-band e multibanda"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("modo", "Implementação", MODOS, ""),
            ("f1", "f1", 0.9, "GHz"),
            ("fase1", "φ1 (LH < 0)", -90.0, "°"),
            ("z1", "Z em f1", 35.35, "Ω"),
            ("f2", "f2", 1.8, "GHz"),
            ("fase2", "φ2", 90.0, "°"),
            ("z2", "Z em f2", 35.35, "Ω"),
            ("n", "Nº de células PLH", 1, ""),
            ("fz", "fz (CSRR)", 0.5, "GHz"),
            ("flinha", "Trechos de linha: φL(f1) (CSRR)", 0.0, "°"),
            ("f3", "f3 (quad-band)", 1.57, "GHz"),
            ("f4", "f4 (quad-band)", 1.8, "GHz"),
            ("fmax", "f máx", 3.0, "GHz"),
        ], titulo="Especificação (φ = comprimento elétrico)")
        self.form.pack(fill="x")
        self.botoes(("Calcular", self.calcular), ("Ex. CSRR (p. 258)", self.ex_csrr),
                    ("Ex. CSRR + linhas (p. 259)", self.ex_csrr_linhas), ("Ex. CRLH 2,4/3,75 GHz", self.ex_crlh),
                    ("Ex. Wilkinson 1/3,1 GHz", self.ex_plh), ("Ex. quad-band GSM/GPS", self.ex_quad),
                    ("Exportar S2P", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 20)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Inversores de impedância dual-band (λ/4 em duas\n"
            "frequências): φ = -90° em f1 (banda LH) e +90° em f2\n"
            "(banda RH), com a mesma impedância.\n\n"
            "Célula T: Xbraço = Z·tan(φ/2), Xderivação = -Z/sen φ\n"
            "(inversão das eqs. 2.30 e 2.33; eqs. 4.28).")
        self._ultimo = None

    # ------------------------------------------------------------------ exemplos
    def _ex(self, **v):
        self.form.set(v)
        self.calcular()

    def ex_csrr(self):
        """Livro, p. 258: L = 12,5 nH, C = 24,9 pF, Cg = 1,25 pF, Lc = 3,38 nH, Cc = 5,10 pF."""
        self._ex(modo=MODOS[2], f1=0.9, fase1=-90, z1=35.35, f2=1.8, fase2=90, z2=35.35, fz=0.5, flinha=0, fmax=3)

    def ex_csrr_linhas(self):
        """Livro, p. 259: φL(f1) = 15° → L = 9,45 nH, C = 17,9 pF, Cg = 1,01 pF, Lc = 4,85 nH, Cc = 2,95 pF."""
        self._ex(modo=MODOS[2], f1=0.9, fase1=-90, z1=35.35, f2=1.8, fase2=90, z2=35.35, fz=0.5, flinha=15, fmax=3)

    def ex_crlh(self):
        """Livro, p. 255 (divisor com OSRR/OCSRR), antes da correção das parasitas."""
        self._ex(modo=MODOS[1], f1=2.4, fase1=-90, z1=35.35, f2=3.75, fase2=90, z2=35.35, fmax=6)

    def ex_plh(self):
        """Livro, p. 252: Wilkinson dual-band, -90° em 1 GHz e +90° em 3,1 GHz (eqs. 4.25-4.27)."""
        self._ex(modo=MODOS[0], f1=1.0, fase1=-90, z1=70.71, f2=3.1, fase2=90, z2=70.71, n=4, fmax=5)

    def ex_quad(self):
        """Livro, p. 267: GSM 0,9/1,8 GHz e GPS 1,17/1,57 GHz, ZA = 35,35 Ω."""
        self._ex(modo=MODOS[3], f1=0.9, f2=1.17, f3=1.57, f4=1.8, z1=35.35, fmax=3)

    # ------------------------------------------------------------------ projeto
    def _projeto(self, v):
        """Retorna (função ABCD(f), linhas de texto, frequências e fases-alvo, Z de referência)."""
        sub = self.estado.substrato
        f1, f2 = v["f1"] * 1e9, v["f2"] * 1e9
        p1, p2 = np.radians(v["fase1"]), np.radians(v["fase2"])
        if v["modo"] == MODOS[0]:
            if abs(v["z1"] - v["z2"]) > 1e-9:
                raise ValueError("Nesta implementação a impedância é a mesma nas duas frequências (Z = √(LL/CL)).")
            r = mb.rh_plh(f1, p1, f2, p2, v["z1"], v["n"])
            w = sintese(v["z1"], sub)
            _, ee = analise(w, sub)
            comp = r.atraso * C0 / np.sqrt(ee) * 1e3

            def m(f):
                lin = abcd.linha(v["z1"], 2 * np.pi * np.asarray(f) * r.atraso)
                plh = abcd.potencia(abcd.celula(f, "T", cs=r.CL, lp=r.LL), r.n)
                return abcd.cascata(lin, plh)
            txt = [f"Linha convencional: atraso l/vp = {r.atraso * 1e12:.2f} ps "
                   f"→ {comp:.3f} mm em {sub.nome} (W = {w:.3f} mm)",
                   f"PLH com {r.n} célula(s) T: CL = {r.CL * 1e12:.4g} pF, LL = {r.LL * 1e9:.4g} nH (eqs. 4.22-4.23)",
                   "Fórmulas supõem PLH homogênea (fase pequena por célula)."]
            return m, txt, [(f1, p1, v["z1"]), (f2, p2, v["z2"])], v["z1"]
        if v["modo"] == MODOS[1]:
            c = mb.crlh_dual(f1, p1, v["z1"], f2, p2, v["z2"])
            txt = [f"LR = {c.LR * 1e9:.4g} nH, CR = {c.CR * 1e12:.4g} pF, LL = {c.LL * 1e9:.4g} nH, "
                   f"CL = {c.CL * 1e12:.4g} pF",
                   f"fs = {c.ws / 2e9 / np.pi:.4f} GHz, fp = {c.wp / 2e9 / np.pi:.4f} GHz "
                   f"(ωs·ωp = ω1·ω2, eq. 4.30)" + (" — balanceada (eq. 4.31)" if c.balanceada else "")]
            return c.abcd, txt, [(f1, p1, v["z1"]), (f2, p2, v["z2"])], v["z1"]
        if v["modo"] == MODOS[2]:
            q1, q2 = p1, p2
            if v["flinha"]:
                q1, q2 = mb.compensar_linhas(p1, p2, f1, f2, v["flinha"])
            sols = mb.csrr_dual(f1, q1, v["z1"], f2, q2, v["z2"], v["fz"] * 1e9)
            c = sols[0]
            tl = np.radians(v["flinha"])

            def m(f):
                f = np.asarray(f, dtype=float)
                if not v["flinha"]:
                    return c.abcd(f)
                lin = abcd.linha(v["z1"], tl * f / f1)
                return abcd.cascata(lin, c.abcd(f), lin)
            txt = [f"L = {c.L * 1e9:.4g} nH, Cg = {c.Cg * 1e12:.4g} pF, C = {c.C * 1e12:.4g} pF, "
                   f"Lc = {c.Lc * 1e9:.4g} nH, Cc = {c.Cc * 1e12:.4g} pF",
                   f"fz = {c.fz / 1e9:.4f} GHz (eq. 3.86), fo(CSRR) = {c.fp / 1e9:.4f} GHz"]
            if v["flinha"]:
                txt.append(f"Com trechos de linha (Z = {v['z1']:g} Ω, φL(f1) = {v['flinha']:g}°) nas pontas: "
                           f"a célula tem φ = {np.degrees(q1):.1f}° e {np.degrees(q2):.1f}° (livro, p. 259)")
            if len(sols) > 1:
                txt.append(f"({len(sols)} soluções; mostrada a 1ª)")
            return m, txt, [(f1, p1, v["z1"]), (f2, p2, v["z2"])], v["z1"]
        freqs = np.array([v["f1"], v["f2"], v["f3"], v["f4"]]) * 1e9
        sols = mb.ecrlh_quad(freqs, v["z1"])
        e = sols[0]
        par = e.params()
        txt = ["E-CRLH (Fig. 3.25c): " + ", ".join(f"{k} = {x * (1e9 if k[0] == 'L' else 1e12):.4g} "
                                                    f"{'nH' if k[0] == 'L' else 'pF'}" for k, x in par.items())]
        if len(sols) > 1:
            txt.append(f"({len(sols)} soluções; mostrada a 1ª)")
        alvos = [(fx, (-1) ** k * np.pi / 2, v["z1"]) for k, fx in enumerate(np.sort(freqs), start=1)]
        return e.abcd, txt, alvos, v["z1"]

    @protegido
    def calcular(self):
        v = self.form.get()
        m_func, txt, alvos, z = self._projeto(v)
        f = np.linspace(v["fmax"] * 1e9 / 3000, v["fmax"] * 1e9, 3000)
        m = m_func(f)
        s = abcd.abcd_para_s(m, z)
        zb, _ = abcd.impedancia_bloch(m)
        fase = -np.angle(s[:, 1, 0])                    # comprimento elétrico = -∠S21, em (-180°, 180°]
        self._ultimo = (f, s, z, v["modo"])
        self.estado.publicar_resultado("linha_multibanda", f, s, z)
        fg = f / 1e9
        ax = self.grafico.eixos(1, 3)
        ax[0].plot(fg, np.degrees(fase), ".", ms=1.5, label="-∠S21 (portas Z)")
        ax[0].set_ylim(-185, 185)
        ax[1].plot(fg, zb.real, label="Re ZB")
        ax[1].plot(fg, zb.imag, alpha=0.6, label="Im ZB")
        ax[2].plot(fg, abcd.db(s[:, 1, 0]), label="|S21|")
        ax[2].plot(fg, abcd.db(s[:, 0, 0]), label="|S11|")
        linhas_alvo = []
        for fx, px, zx in alvos:
            k = np.argmin(abs(f - fx))
            ax[0].plot(fx / 1e9, np.degrees(px), "ko", ms=5)
            ax[1].plot(fx / 1e9, zx, "ko", ms=5)
            linhas_alvo.append(f"  {fx / 1e9:7.4f}  {np.degrees(px):7.1f}  {np.degrees(fase[k]):7.2f}"
                               f"  {zx:7.2f}  {zb[k].real:7.2f}  {abcd.db(s[k, 0, 0]):7.1f}")
        ax[0].set_ylabel("Comprimento elétrico [°]")
        ax[0].set_title("Fase (pontos: alvos)", fontsize=9)
        ax[1].set_ylim(-3 * z, 3 * z)
        ax[1].set_ylabel("ZB [Ω]")
        ax[1].set_title("Impedância de Bloch", fontsize=9)
        ax[2].set_ylim(-50, 2)
        ax[2].set_title(f"S com portas de {z:g} Ω", fontsize=9)
        for a in ax:
            a.set_xlabel("f [GHz]")
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        self.saida.escrever("\n".join([v["modo"], ""] + txt + [
            "", "Verificação nas frequências de projeto:",
            "  f[GHz]   φ alvo  φ obtido  Z alvo  Re ZB  |S11|dB"] + linhas_alvo + [
            "", "Modelo de circuito ideal. A implementação física",
            "(OSRR/OCSRR, CSRR, interdigital, SMT) acrescenta",
            "parasitas: o livro recalcula os elementos com elas",
            "e fecha o layout por extração ou ASM (Ap. G e H)."]))

    def exportar_s2p(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s, z, modo = self._ultimo
            self.salvar_s2p("linha_multibanda", f, s, z, f"Linha multibanda: {modo}")
