"""Aba CRLH / D-CRLH: síntese e análise de linhas metamateriais canônicas (seções 3.4.2-3.4.3)."""
import numpy as np

from ..core import abcd
from ..core import crlh
from ..core import floquet as fl
from ..core.microstrip import C0
from .comum import AbaBase, AreaTexto, Formulario, protegido
from .previa import PreviaCelula

MODOS = ("Balanceada: banda f_CL a f_CR", "Balanceada: fases em f1 e f2", "Elementos LR, CR, LL, CL")


class AbaCRLH(AbaBase):
    TITULO = "CRLH / D-CRLH"
    PREVIA = PreviaCelula

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("tipo", "Linha", crlh.TIPOS, ""),
            ("modo", "Especificação", MODOS, ""),
            ("zb", "ZB (transição)", 50.0, "Ω"),
            ("fcl", "f_CL (borda inferior)", 1.0, "GHz"),
            ("fcr", "f_CR (borda superior)", 6.0, "GHz"),
            ("f1", "f1", 1.0, "GHz"),
            ("fase1", "βl em f1 (LH < 0)", -45.0, "°"),
            ("f2", "f2", 2.0, "GHz"),
            ("fase2", "βl em f2", 22.5, "°"),
        ], titulo="Especificação")
        self.form.pack(fill="x")
        self.el = Formulario(self.painel, [
            ("LR", "LR (série, RH)", 3.18, "nH"),
            ("CR", "CR (derivação, RH)", 1.27, "pF"),
            ("LL", "LL (derivação, LH)", 3.32, "nH"),
            ("CL", "CL (série, LH)", 1.33, "pF"),
        ], titulo="Elementos da célula (Fig. 3.20)")
        self.el.pack(fill="x", pady=(8, 0))
        self.an = Formulario(self.painel, [
            ("n", "Nº de células", 9, ""),
            ("zp", "Z das portas", 50.0, "Ω"),
            ("lcel", "Comprimento da célula (0 = ignorar)", 6.0, "mm"),
            ("alfa", "LWA: α/k0 (de EM; 0 = ignorar)", 0.0, ""),
            ("eta", "LWA: potência a irradiar", 95.0, "%"),
            ("fmax", "f máx", 8.0, "GHz"),
            ("npts", "Pontos", 2001, ""),
        ], titulo="Análise")
        self.an.pack(fill="x", pady=(8, 0))
        self.botoes(("Calcular", self.calcular), ("Ex. balanceada 1-6 GHz", self.ex_banda),
                    ("Ex. não balanceada", self.ex_desbalanceada), ("Ex. dual-band (fases)", self.ex_fases),
                    ("Exportar S2P", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 18)
        self.saida.pack(fill="both", expand=True)
        self._ultimo = None
        self.observar_previa(self.form, self.el)
        self.atualizar_previa()

    # ------------------------------------------------------------------ exemplos
    def ex_banda(self):
        self.form.set(dict(tipo="CRLH", modo=MODOS[0], zb=50, fcl=1, fcr=6))
        self.an.set(dict(n=9, zp=50, lcel=6, fmax=8))
        self.calcular()

    def ex_desbalanceada(self):
        """ωs < ωp, como na Fig. 3.21: gap entre as bandas LH e RH."""
        self.form.set(dict(tipo="CRLH", modo=MODOS[2]))
        self.el.set(dict(LR=2.5, CR=1.0, LL=3.0, CL=2.0))
        self.an.set(dict(n=7, zp=50, lcel=6, fmax=8))
        self.calcular()

    def ex_fases(self):
        """Linha dual-band: -90° em 0,9 GHz e +90° em 1,8 GHz com 2 células (-45° e +45° por célula)."""
        self.form.set(dict(tipo="CRLH", modo=MODOS[1], zb=50, f1=0.9, fase1=-45, f2=1.8, fase2=45))
        self.an.set(dict(n=2, zp=50, lcel=0, fmax=4))
        self.calcular()

    # ------------------------------------------------------------------ projeto
    def _celula(self):
        v = self.form.get()
        if v["modo"] == MODOS[0]:
            c = crlh.balanceada_por_banda(v["fcl"] * 1e9, v["fcr"] * 1e9, v["zb"])
        elif v["modo"] == MODOS[1]:
            c = crlh.balanceada_por_fases(v["f1"] * 1e9, np.radians(v["fase1"]), v["f2"] * 1e9,
                                          np.radians(v["fase2"]), v["zb"])
        else:
            e = self.el.get()
            c = crlh.CelulaCRLH(e["LR"] * 1e-9, e["CR"] * 1e-12, e["LL"] * 1e-9, e["CL"] * 1e-12)
        c = crlh.CelulaCRLH(c.LR, c.CR, c.LL, c.CL, v["tipo"])
        if v["modo"] != MODOS[2]:
            self.el.set(dict(LR=f"{c.LR * 1e9:.4g}", CR=f"{c.CR * 1e12:.4g}", LL=f"{c.LL * 1e9:.4g}",
                             CL=f"{c.CL * 1e12:.4g}"))
        return v, c

    def atualizar_previa(self):
        try:
            v, c = self._celula()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        n = self.an.get()["n"] if self.an.get()["n"] > 0 else 1
        if c.tipo == "CRLH":
            self.previa.mostrar("T", n, ls=c.LR, cs=c.CL, lp=c.LL, cp=c.CR)
        else:
            self.previa.mostrar_erro("D-CRLH (Fig. 3.23): braços série (LR/2 ‖ 2CL),\n"
                                     "derivação LL em série com CR.")

    @protegido
    def calcular(self):
        v, c = self._celula()
        a = self.an.get()
        if a["n"] < 1 or a["fmax"] <= 0 or a["npts"] < 50:
            raise ValueError("Use ao menos 1 célula, f máx > 0 e 50 pontos.")
        self.atualizar_previa()
        f = np.linspace(a["fmax"] * 1e9 / a["npts"], a["fmax"] * 1e9, a["npts"])
        bl, al = c.beta_l(f)
        zb = c.zb(f)
        s = abcd.abcd_para_s(abcd.potencia(c.abcd(f), a["n"]), a["zp"])
        self._ultimo = (f, s, a["zp"], c, a)
        self.estado.publicar_resultado("crlh", f, s, a["zp"])
        fg = f / 1e9
        passa = al < 1e-6
        lh = passa & (bl < 0)
        rh = passa & (bl >= 0)
        lcel = a["lcel"] * 1e-3

        ax = self.grafico.eixos(2, 3)
        if lcel > 0:
            k0l = 2 * np.pi * f / C0 * lcel
            ax[0].fill_betweenx(fg, -k0l / np.pi, k0l / np.pi, color="#cfe3f5", label="Onda rápida (radia)")
        ax[0].plot(np.where(lh, bl, np.nan) / np.pi, fg, color="C3", lw=1.8, label="LH (βl < 0)")
        ax[0].plot(np.where(rh, bl, np.nan) / np.pi, fg, color="C0", lw=1.8, label="RH (βl > 0)")
        ax[0].plot(np.where(passa, np.nan, al) / np.pi, fg, color="gray", ls="--", lw=1, label="αl (bandas rejeitadas)")
        ax[0].set_xlim(-1.05, 1.05)
        ax[0].set_ylim(0, fg[-1])
        ax[0].set_xlabel("βl/π, αl/π")
        ax[0].set_ylabel("f [GHz]")
        ax[0].set_title(f"Dispersão {c.tipo} (eq. {'3.51' if c.tipo == 'CRLH' else '3.64'})", fontsize=9)

        ax[1].plot(fg, zb.real, label="Re ZB")
        ax[1].plot(fg, zb.imag, alpha=0.6, label="Im ZB")
        ax[1].set_ylim(-3 * v["zb"], 3 * v["zb"])
        ax[1].set_ylabel("ZB [Ω]")
        ax[1].set_title("Impedância de Bloch", fontsize=9)

        ax[2].plot(fg, abcd.db(s[:, 1, 0]), label="|S21|")
        ax[2].plot(fg, abcd.db(s[:, 0, 0]), label="|S11|")
        ax[2].set_ylim(-60, 2)
        ax[2].set_ylabel("|S| [dB]")
        ax[2].set_title(f"{a['n']} células, portas {a['zp']:g} Ω", fontsize=9)

        if c.tipo == "CRLH" and lcel > 0:
            eps, mu = c.meio_efetivo(f, lcel)
            ax[3].plot(fg, eps, label="εr,eff (eq. 3.61)")
            ax[3].plot(fg, mu, label="μr,eff (eq. 3.62)")
            ax[3].axhline(0, color="k", lw=0.6)
            lim = 3 * max(abs(eps[-1]), abs(mu[-1]), 1)
            ax[3].set_ylim(-lim, lim)
            ax[3].set_title("Parâmetros constitutivos efetivos", fontsize=9)
        else:
            ax[3].text(0.5, 0.5, "Informe o comprimento da célula\n(só CRLH)", ha="center", va="center",
                       transform=ax[3].transAxes, fontsize=8)
            ax[3].set_title("Parâmetros constitutivos efetivos", fontsize=9)

        fase = np.degrees(np.unwrap(np.angle(s[:, 1, 0])))
        ax[4].plot(fg, np.where(passa, -bl * a["n"], np.nan) * 180 / np.pi, label="-N·βl (Bloch)")
        ax[4].plot(fg, fase - fase[np.argmin(np.abs(bl))] , ls=":", label="∠S21 (desenrolada)")
        ax[4].set_ylabel("Fase [°]")
        ax[4].set_title("Fase de transmissão (avanço na banda LH)", fontsize=9)

        texto_leaky = []
        if lcel > 0:
            h = fl.harmonicos(f, np.where(passa, bl, np.nan), lcel, 0)[0]
            ax[5].plot(fg, h.theta_graus, color="C2")
            ax[5].axhline(90, color="gray", ls=":", lw=0.8)
            ax[5].set_ylim(0, 180)
            ax[5].set_ylabel("θ [°] (cos θ = β/k0)")
            ax[5].set_title("Feixe leaky-wave do modo fundamental (eq. 2.11)", fontsize=9)
            for fa, fb, ta, tb in fl.faixas_radiacao(f, h):
                texto_leaky.append(f"  {fa / 1e9:.3f} a {fb / 1e9:.3f} GHz: θ de {ta:.0f}° a {tb:.0f}°")
            if a["alfa"] > 0 and texto_leaky:
                texto_leaky += self._texto_lwa(c, f, h, lcel, a)
        else:
            ax[5].text(0.5, 0.5, "Informe o comprimento da célula\npara a análise leaky-wave", ha="center",
                       va="center", transform=ax[5].transAxes, fontsize=8)
        for x in ax[1:]:
            x.set_xlabel("f [GHz]")
        for x in ax:
            x.grid(True, alpha=0.3)
            if x.get_legend_handles_labels()[0]:
                x.legend(fontsize=7)
        self.grafico.desenhar()
        self.ver_resultados()
        self.saida.escrever(self._texto(c, f, bl, al, zb, s, a, texto_leaky))

    def _texto_lwa(self, c, f, h, lcel, a):
        """Projeto da antena: nº de células (eq. 4.51) e largura de feixe."""
        f0 = crlh.frequencia_transicao(c) if c.tipo == "CRLH" else f[np.nanargmin(np.abs(h.theta_graus - 90))]
        n = fl.celulas_lwa(a["alfa"], f0, lcel, a["eta"] / 100)
        comp = n * lcel
        linhas = ["", f"Antena leaky-wave (α/k0 = {a['alfa']:g}, {a['eta']:g}% irradiado):",
                  f"  N = {n} células (eq. 4.51), comprimento {comp * 1e3:.1f} mm "
                  f"= {comp * f0 / C0:.2f} λ0 em {f0 / 1e9:.3f} GHz"]
        for ang in (135, 90, 45):
            i = np.nanargmin(np.abs(np.nan_to_num(h.theta_graus, nan=1e9) - ang))
            if abs(h.theta_graus[i] - ang) < 5:
                linhas.append(f"  θ = {ang}° em {f[i] / 1e9:.3f} GHz: feixe ≈ "
                              f"{fl.largura_feixe(comp, f[i], ang):.1f}° (Jackson-Oliner)")
        linhas.append("  α/k0 depende da geometria (simular N células em EM).")
        return linhas

    def _texto(self, c, f, bl, al, zb, s, a, texto_leaky):
        g = lambda w: w / (2 * np.pi * 1e9)  # noqa: E731
        linhas = [f"{c.tipo}: LR = {c.LR * 1e9:.4g} nH, CR = {c.CR * 1e12:.4g} pF, "
                  f"LL = {c.LL * 1e9:.4g} nH, CL = {c.CL * 1e12:.4g} pF",
                  f"fR = {g(c.wR):.4f}  fL = {g(c.wL):.4f}  fs = {g(c.ws):.4f}  fp = {g(c.wp):.4f} GHz",
                  f"√(LR/CR) = {c.z_r:.2f} Ω   √(LL/CL) = {c.z_l:.2f} Ω",
                  "Balanceada (ωs = ωp)" if c.balanceada else
                  f"Não balanceada: gap de {g(c.gap[0]):.4f} a {g(c.gap[1]):.4f} GHz (eqs. 3.57-3.58)"]
        if c.tipo == "CRLH":
            b = c.bandas()
            linhas.append(f"Banda total: {b['f_cl'] / 1e9:.4f} a {b['f_cr'] / 1e9:.4f} GHz"
                          + (" (eq. 3.66)" if c.balanceada else ""))
            f0 = crlh.frequencia_transicao(c)
            if c.balanceada:
                linhas.append(f"Transição f0 = {f0 / 1e9:.4f} GHz, ZB(f0) = {c.z_r:.2f} Ω (eq. 3.63)")
        passa = al < 1e-6
        for nome, sel in (("LH", passa & (bl < 0)), ("RH", passa & (bl >= 0))):
            if sel.any():
                linhas.append(f"Banda {nome} (calculada): {f[sel][0] / 1e9:.4f} a {f[sel][-1] / 1e9:.4f} GHz")
        linhas += ["", abcd.descrever_banda(f, s)]
        if texto_leaky:
            linhas += ["", f"Radiação leaky-wave (célula de {a['lcel']:g} mm, |β| < k0):"] + texto_leaky
            linhas.append("  LH: para trás (θ > 90°); transição: broadside; RH: para a frente.")
        linhas += ["", "Modelo de circuito ideal (sem perdas, sem parasitas",
                   "de implementação). A aproximação de grande λ",
                   "(eq. 3.59) vale para |βl| pequeno por célula."]
        return "\n".join(linhas)

    def exportar_s2p(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s, z, c, a = self._ultimo
            self.salvar_s2p(f"{c.tipo.lower()}_{a['n']}cel", f, s, z,
                            f"{c.tipo}: LR={c.LR * 1e9:.4g}nH CR={c.CR * 1e12:.4g}pF LL={c.LL * 1e9:.4g}nH "
                            f"CL={c.CL * 1e12:.4g}pF N={a['n']}")
