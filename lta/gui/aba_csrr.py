"""Aba Linha CSRR + gap (microstrip): projeto do circuito e geometria inicial (seção 3.5.2.2)."""
import numpy as np

from ..core import abcd
from ..core import csrr as cs
from ..geometria import construtores
from .aba_substrato import PREDEFINIDOS
from .comum import AbaBase, AreaTexto, Formulario, protegido
from .previa import PreviaEstrutura

MODOS = ("Balanceada: f0, ZB e fz", "Elementos do modelo (Fig. 3.38b)")
SERIE = ("Gap simples", "Interdigital")


class AbaCSRR(AbaBase):
    TITULO = "Linha CSRR + gap (microstrip)"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("modo", "Especificação", MODOS, ""),
            ("f0", "f0 (transição LH/RH)", 2.0, "GHz"),
            ("zb", "ZB(f0)", 50.0, "Ω"),
            ("fz", "fz (zero de transmissão)", 1.4, "GHz"),
            ("periodo", "Período da célula", 8.0, "mm"),
            ("zlinha", "Z0 da linha", 50.0, "Ω"),
        ], titulo="Especificação elétrica")
        self.form.pack(fill="x")
        self.el = Formulario(self.painel, [
            ("L", "L (linha)", 3.48, "nH"),
            ("Cg", "Cg (série)", 1.82, "pF"),
            ("C", "C (acoplamento)", 1.45, "pF"),
            ("Lc", "Lc (CSRR)", 4.55, "nH"),
            ("Cc", "Cc (CSRR)", 1.39, "pF"),
        ], titulo="Elementos do modelo (Fig. 3.38b)")
        self.el.pack(fill="x", pady=(8, 0))
        self.geo = Formulario(self.painel, [
            ("c", "Fendas c", 0.3, "mm"),
            ("d", "Distância d", 0.3, "mm"),
            ("abertura", "Corte dos anéis", 0.3, "mm"),
            ("calib", "Calibração fo (real/estimada)", 1.0, ""),
            ("serie", "Elemento série", SERIE, ""),
            ("gap", "Gap", 0.2, "mm"),
            ("ndedos", "Nº de dedos", 6, ""),
            ("comp", "Comprimento dos dedos", 1.0, "mm"),
            ("larg", "Largura/espaço dos dedos", 0.1, "mm"),
            ("n", "Nº de células", 3, ""),
            ("lalim", "Alimentação", 5.0, "mm"),
            ("fmax", "f máx", 5.0, "GHz"),
        ], titulo="Geometria inicial e análise")
        self.geo.pack(fill="x", pady=(8, 0))
        self.botoes(("Calcular", self.calcular), ("Ex. Tabela 3.3 (Fig. 3.39)", self.ex_tabela),
                    ("Ex. balanceada 2 GHz", self.ex_balanceada), ("Enviar p/ exportação", self.publicar),
                    ("Exportar S2P", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 20)
        self.saida.pack(fill="both", expand=True)
        self._ultimo = None
        self.observar_previa(self.form, self.el, self.geo)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    # ------------------------------------------------------------------ exemplos
    def _ro3010(self):
        self.estado.definir_substrato(PREDEFINIDOS["Rogers RO3010 (livro)"])

    def ex_tabela(self):
        """Fig. 3.39a / Tabela 3.3: elementos extraídos e fendas do livro (rext real = 3,661 mm)."""
        self._ro3010()
        self.form.set(dict(modo=MODOS[1], periodo=8.0, zlinha=50))
        self.el.set(dict(L=5.78, Cg=0.31, C=6.43, Lc=1.49, Cc=3.10))
        self.geo.set(dict(c=0.365, d=0.24, abertura=0.3, calib=1.0, serie=SERIE[0], gap=0.16, n=1, fmax=4))
        self.calcular()

    def ex_balanceada(self):
        self._ro3010()
        self.form.set(dict(modo=MODOS[0], f0=2.0, zb=50, fz=1.4, periodo=8.0, zlinha=50))
        self.geo.set(dict(c=0.3, d=0.3, abertura=0.3, calib=1.0, serie=SERIE[1], ndedos=8, comp=2.0, larg=0.1,
                          n=3, fmax=5))
        self.calcular()

    # ------------------------------------------------------------------ projeto
    def _gerar(self):
        v, g = self.form.get(), self.geo.get()
        sub = self.estado.substrato
        l_linha, w = cs.indutancia_linha(v["zlinha"], sub, v["periodo"])
        if v["modo"] == MODOS[0]:
            circ = cs.projeto_balanceado(v["f0"] * 1e9, v["zb"], v["fz"] * 1e9, l_linha)
            self.el.set({k: f"{x * (1e9 if k[0] == 'L' else 1e12):.4g}" for k, x in circ.params().items()})
        else:
            e = self.el.get()
            if min(e.values()) <= 0:
                raise ValueError("Todos os elementos devem ser positivos.")
            circ = cs.CircuitoCSRR(e["L"] * 1e-9, e["Cg"] * 1e-12, e["C"] * 1e-12, e["Lc"] * 1e-9, e["Cc"] * 1e-12)
        if g["calib"] <= 0 or g["n"] < 1:
            raise ValueError("Calibração > 0 e ao menos 1 célula.")
        r = cs.raio_para_fo(circ.fp, g["c"], g["d"], sub, g["abertura"], g["calib"])
        geo = cs.GeometriaCSRR(r, g["c"], g["d"], g["abertura"])
        cs_alvo, cpar = circ.modelo_melhorado()
        info = {"tipo": "linha CSRR + elemento série", "L_nH": f"{circ.L * 1e9:.4g}", "Cg_pF": f"{circ.Cg * 1e12:.4g}",
                "C_pF": f"{circ.C * 1e12:.4g}", "Lc_nH": f"{circ.Lc * 1e9:.4g}", "Cc_pF": f"{circ.Cc * 1e12:.4g}",
                "fo_CSRR_GHz": f"{circ.fp / 1e9:.4f}", "r_ext_mm": f"{r:.4f}", "c_mm": g["c"], "d_mm": g["d"],
                "periodo_mm": v["periodo"], "W_mm": f"{w:.4f}", "serie": g["serie"], "celulas": g["n"]}
        lay = construtores.linha_csrr(g["n"], v["periodo"], w, geo, sub, g["serie"], g["gap"], g["ndedos"],
                                      g["comp"], g["larg"], g["lalim"], nome=f"csrr_{g['n']}cel", info=info)
        return v, g, circ, geo, w, (cs_alvo, cpar), lay

    def atualizar_previa(self):
        try:
            *_, lay = self._gerar()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        self.previa.mostrar(lay)

    @protegido
    def calcular(self):
        v, g, circ, geo, w, (cs_alvo, cpar), lay = self._gerar()
        self.previa.mostrar(lay)
        sub = self.estado.substrato
        f = np.linspace(g["fmax"] * 1e9 / 2000, g["fmax"] * 1e9, 2000)
        m = circ.abcd(f)
        al, bl = abcd.dispersao(m)
        wf = 2 * np.pi * f
        braco = wf * circ.L / 2 - 1 / (2 * wf * circ.Cg)
        bl = np.where(braco < 0, -bl, bl)                        # banda LH: braço série capacitivo
        zb = circ.zb(f)
        s = abcd.abcd_para_s(abcd.potencia(m, g["n"]), v["zlinha"])
        s1 = abcd.abcd_para_s(m, v["zlinha"])
        self._ultimo = (f, s, v["zlinha"], g)
        self.estado.publicar_resultado("linha_csrr", f, s, v["zlinha"])
        self.estado.publicar_resultado("celula_csrr", f, s1, v["zlinha"])
        fg = f / 1e9
        passa = al < 1e-6
        ax = self.grafico.eixos(2, 2)
        ax[0].plot(np.where(passa & (bl < 0), bl, np.nan) / np.pi, fg, color="C3", lw=1.8, label="LH")
        ax[0].plot(np.where(passa & (bl >= 0), bl, np.nan) / np.pi, fg, color="C0", lw=1.8, label="RH")
        ax[0].plot(np.where(passa, np.nan, al) / np.pi, fg, color="gray", ls="--", lw=1, label="αl")
        ax[0].set_xlim(-1.05, 1.05)
        ax[0].set_xlabel("βl/π, αl/π")
        ax[0].set_ylabel("f [GHz]")
        ax[0].set_title("Dispersão (eq. 3.89)", fontsize=9)
        ax[1].plot(fg, zb.real, label="Re ZB")
        ax[1].plot(fg, zb.imag, alpha=0.6, label="Im ZB")
        ax[1].set_ylim(-150, 150)
        ax[1].set_ylabel("ZB [Ω]")
        ax[1].set_title("Impedância de Bloch (eq. 3.90)", fontsize=9)
        ax[2].plot(fg, abcd.db(s[:, 1, 0]), label="|S21|")
        ax[2].plot(fg, abcd.db(s[:, 0, 0]), label="|S11|")
        ax[2].set_ylim(-60, 2)
        ax[2].set_title(f"{g['n']} células (modelo de circuito)", fontsize=9)
        ax[3].plot(fg, abcd.db(s1[:, 1, 0]), label="|S21| 1 célula")
        for nome, fx in (("fz", circ.fz), ("fo", circ.fp), ("fs", circ.fs)):
            if fx < f[-1]:
                ax[3].axvline(fx / 1e9, color="gray", ls=":", lw=0.8)
                ax[3].annotate(nome, (fx / 1e9, 1), fontsize=7, ha="center", color="gray")
        ax[3].set_ylim(-60, 3)
        ax[3].set_title("Célula: zero de transmissão e ressonâncias", fontsize=9)
        for a in ax[1:]:
            a.set_xlabel("f [GHz]")
        for a in ax:
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        self.ver_resultados()
        self.saida.escrever(self._texto(v, g, circ, geo, w, cs_alvo, cpar, sub, f, s))

    def _texto(self, v, g, circ, geo, w, cs_alvo, cpar, sub, f, s):
        lc_e, cc_e, fo_e = cs.estimativa_csrr(geo, sub)
        linhas = ["Modelo (Fig. 3.38b):",
                  f"  L = {circ.L * 1e9:.4g} nH  Cg = {circ.Cg * 1e12:.4g} pF  C = {circ.C * 1e12:.4g} pF",
                  f"  Lc = {circ.Lc * 1e9:.4g} nH  Cc = {circ.Cc * 1e12:.4g} pF",
                  f"  fz = {circ.fz / 1e9:.4f} GHz (eq. 3.86)  fo(CSRR) = {circ.fp / 1e9:.4f} GHz  "
                  f"fs = {circ.fs / 1e9:.4f} GHz",
                  "  Balanceada (eq. 3.91)" if circ.balanceada else
                  f"  Não balanceada: fs ≠ fo (gap entre LH e RH)",
                  "", "Modelo melhorado (Fig. 3.40, eqs. 3.87-3.88):",
                  f"  Cs (série física) = {cs_alvo * 1e12:.4g} pF   Cpar = Cf + CL = {cpar * 1e12:.4g} pF", "",
                  f"Geometria inicial ({sub.nome}), linha W = {w:.3f} mm, período {v['periodo']:g} mm:",
                  f"  CSRR: rext = {geo.r_ext:.3f} mm, c = {geo.c:g}, d = {geo.d:g} mm "
                  f"(estimativa fo = {fo_e / 1e9 * g['calib']:.4f} GHz)"]
        if g["serie"] == SERIE[1]:
            try:
                c_idc = cs.capacitancia_interdigital(g["ndedos"], g["comp"], g["larg"], sub)
                linhas.append(f"  Interdigital: {g['ndedos']} dedos de {g['comp']:g} mm → ≈ {c_idc * 1e12:.3g} pF "
                              f"(Bahl); para Cs = {cs_alvo * 1e12:.3g} pF use ≈ "
                              f"{g['comp'] * cs_alvo / c_idc:.3g} mm")
            except ValueError as e:
                linhas.append(f"  Interdigital: {e}")
        else:
            linhas.append(f"  Gap de {g['gap']:g} mm: capacitância a confirmar por EM "
                          f"(alvo Cs = {cs_alvo * 1e12:.3g} pF)")
        linhas += ["", abcd.descrever_banda(f, s), "",
                   "⚠ Geometria = estimativa quase-estática (±20% em fo).",
                   "Fluxo do livro para fechar o projeto:",
                   " 1. Enviar p/ exportação → HFSS/CST, simular 1 célula;",
                   " 2. aba Extração, modelo 'CSRR + gap', com o .s2p;",
                   " 3. calibração = fo extraído / fo estimado; recalcular.",
                   "Livro, Fig. 3.39: rext = 3,661 mm (c = 0,365; d = 0,24)."]
        return "\n".join(linhas)

    @protegido
    def publicar(self):
        *_, lay = self._gerar()
        self.previa.mostrar(lay)
        self.estado.publicar_layout(lay)
        self.estado.msg(f"Layout '{lay.nome}' enviado para Exportar (Ctrl+E).")

    def exportar_s2p(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s, z, g = self._ultimo
            self.salvar_s2p(f"csrr_{g['n']}cel", f, s, z, "Linha CSRR + elemento serie (modelo Fig. 3.38b)")
