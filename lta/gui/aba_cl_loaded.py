"""Aba CRLH microstrip (CL-loaded): geometria inicial com interdigitais e stubs com via (seção 3.5.1)."""
import numpy as np

from ..core import abcd
from ..core import cl_loaded as cl
from ..core import crlh
from ..core.csrr import capacitancia_interdigital
from ..geometria import construtores
from .aba_substrato import PREDEFINIDOS
from .comum import AbaBase, AreaTexto, Formulario, protegido
from .previa import PreviaEstrutura

MODOS = ("Balanceada: banda f_CL a f_CR", "Elementos LR, CR, LL, CL")


class AbaCLLoaded(AbaBase):
    TITULO = "CRLH microstrip (CL-loaded)"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("modo", "Especificação", MODOS, ""),
            ("zb", "ZB", 50.0, "Ω"),
            ("fcl", "f_CL", 2.0, "GHz"),
            ("fcr", "f_CR", 7.0, "GHz"),
        ], titulo="Especificação elétrica")
        self.form.pack(fill="x")
        self.el = Formulario(self.painel, [
            ("LR", "LR", 2.45, "nH"),
            ("CR", "CR", 0.5, "pF"),
            ("LL", "LL", 3.38, "nH"),
            ("CL", "CL", 0.68, "pF"),
        ], titulo="Elementos da célula (Fig. 3.20)")
        self.el.pack(fill="x", pady=(8, 0))
        self.geo = Formulario(self.painel, [
            ("ndedos", "Nº de dedos do interdigital", 14, ""),
            ("larg", "Largura/espaço dos dedos", 0.1, "mm"),
            ("zstub", "Z0 do stub", 100.0, "Ω"),
            ("rvia", "Raio da via", 0.2, "mm"),
            ("n", "Nº de células", 9, ""),
            ("zp", "Z das portas (0 = √(LR/CR))", 50.0, "Ω"),
            ("lalim", "Alimentação", 5.0, "mm"),
            ("fmax", "f máx", 10.0, "GHz"),
        ], titulo="Implementação microstrip")
        self.geo.pack(fill="x", pady=(8, 0))
        self.botoes(("Calcular", self.calcular), ("Ex. f0 = 3,9 GHz (Fig. 3.28)", self.exemplo),
                    ("Enviar p/ exportação", self.publicar), ("Exportar S2P", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 20)
        self.saida.pack(fill="both", expand=True)
        self._ultimo = None
        self.observar_previa(self.form, self.el, self.geo)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    def exemplo(self):
        """Valores de Caloz e Itoh: transição em ≈ 3,9 GHz, como a linha de 9 células da Fig. 3.28."""
        self.estado.definir_substrato(PREDEFINIDOS["Rogers RT5880 0.787 mm"])
        self.form.set(dict(modo=MODOS[1]))
        self.el.set(dict(LR=2.45, CR=0.5, LL=3.38, CL=0.68))
        self.geo.set(dict(ndedos=14, larg=0.1, zstub=100, rvia=0.2, n=9, fmax=10))
        self.calcular()

    def _gerar(self):
        v, g = self.form.get(), self.geo.get()
        if v["modo"] == MODOS[0]:
            c = crlh.balanceada_por_banda(v["fcl"] * 1e9, v["fcr"] * 1e9, v["zb"])
            self.el.set(dict(LR=f"{c.LR * 1e9:.4g}", CR=f"{c.CR * 1e12:.4g}", LL=f"{c.LL * 1e9:.4g}",
                             CL=f"{c.CL * 1e12:.4g}"))
        else:
            e = self.el.get()
            c = crlh.CelulaCRLH(e["LR"] * 1e-9, e["CR"] * 1e-12, e["LL"] * 1e-9, e["CL"] * 1e-12)
        if g["n"] < 1:
            raise ValueError("Use ao menos 1 célula.")
        sub = self.estado.substrato
        geo, avisos = cl.projeto(c, sub, g["ndedos"], g["larg"], g["zstub"], g["rvia"])
        info = {"tipo": "CRLH CL-loaded", "LR_nH": f"{c.LR * 1e9:.4g}", "CR_pF": f"{c.CR * 1e12:.4g}",
                "LL_nH": f"{c.LL * 1e9:.4g}", "CL_pF": f"{c.CL * 1e12:.4g}", "periodo_mm": f"{geo.periodo:.4f}",
                "dedos": g["ndedos"], "comp_dedo_mm": f"{geo.comp_dedo:.4f}", "stub_mm": f"{geo.comp_stub:.4f}",
                "W_stub_mm": f"{geo.w_stub:.4f}", "raio_via_mm": g["rvia"], "celulas": g["n"]}
        lay = construtores.linha_cl_loaded(geo, g["n"], sub, g["lalim"], nome=f"crlh_cl_{g['n']}cel", info=info)
        return v, g, c, geo, avisos, lay

    def atualizar_previa(self):
        try:
            *_, lay = self._gerar()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        self.previa.mostrar(lay)

    @protegido
    def calcular(self):
        v, g, c, geo, avisos, lay = self._gerar()
        self.previa.mostrar(lay)
        f = np.linspace(g["fmax"] * 1e9 / 2000, g["fmax"] * 1e9, 2000)
        bl, al = c.beta_l(f)
        zp = g["zp"] if g["zp"] > 0 else np.sqrt(c.LR / c.CR)
        s = abcd.abcd_para_s(abcd.potencia(c.abcd(f), g["n"]), zp)
        s1 = abcd.abcd_para_s(c.abcd(f), zp)
        self._ultimo = (f, s, zp, g)
        self.estado.publicar_resultado("crlh_cl_loaded", f, s, zp)
        self.estado.publicar_resultado("celula_crlh_cl", f, s1, zp)
        fg = f / 1e9
        passa = al < 1e-6
        ax = self.grafico.eixos(1, 2)
        ax[0].plot(np.where(passa & (bl < 0), bl, np.nan) / np.pi, fg, color="C3", lw=1.8, label="LH")
        ax[0].plot(np.where(passa & (bl >= 0), bl, np.nan) / np.pi, fg, color="C0", lw=1.8, label="RH")
        ax[0].set_xlim(-1.05, 1.05)
        ax[0].set_xlabel("βl/π")
        ax[0].set_ylabel("f [GHz]")
        ax[0].set_title("Dispersão do circuito-alvo (eq. 3.51)", fontsize=9)
        ax[1].plot(fg, abcd.db(s[:, 1, 0]), label="|S21|")
        ax[1].plot(fg, abcd.db(s[:, 0, 0]), label="|S11|")
        ax[1].set_ylim(-60, 2)
        ax[1].set_xlabel("f [GHz]")
        ax[1].set_title(f"{g['n']} células, portas {zp:.0f} Ω, circuito ideal (cf. Fig. 3.28)", fontsize=9)
        for a in ax:
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        self.ver_resultados()
        sub = self.estado.substrato
        b = c.bandas()
        c_idc = capacitancia_interdigital(g["ndedos"], geo.comp_dedo, g["larg"], sub)
        linhas = [f"Célula: LR = {c.LR * 1e9:.4g} nH, CR = {c.CR * 1e12:.4g} pF, LL = {c.LL * 1e9:.4g} nH, "
                  f"CL = {c.CL * 1e12:.4g} pF",
                  f"Transição f0 = {crlh.frequencia_transicao(c) / 1e9:.3f} GHz, "
                  + ("balanceada" if c.balanceada else f"gap {c.gap[0] / 2e9 / np.pi:.3f}-{c.gap[1] / 2e9 / np.pi:.3f} GHz"),
                  f"Banda total {b['f_cl'] / 1e9:.3f} a {b['f_cr'] / 1e9:.3f} GHz", "",
                  f"Geometria inicial ({sub.nome}):",
                  f"  Linha hospedeira: Z = {np.sqrt(c.LR / c.CR):.1f} Ω, W = {geo.w_linha:.3f} mm, "
                  f"período p = {geo.periodo:.3f} mm (LR)",
                  f"  Interdigital CL: {g['ndedos']} dedos de {geo.comp_dedo:.3f} mm "
                  f"(≈ {c_idc * 1e12:.3g} pF, Bahl); nas pontas 2CL: {geo.comp_dedo_ext:.3f} mm",
                  f"  Stub LL: W = {geo.w_stub:.3f} mm ({g['zstub']:g} Ω), {geo.comp_stub:.3f} mm "
                  f"({geo.l_stub * 1e9:.3g} nH) + via r = {g['rvia']:g} mm ({geo.l_via * 1e9:.3g} nH)",
                  f"  Comprimento total {lay.comprimento:.2f} mm"]
        linhas += [""] + [f"⚠ {a}" for a in avisos]
        linhas += ["", "⚠ Estimativa sem parasitas do interdigital e da junção",
                   "(o livro não dá síntese desta linha). Fluxo para fechar:",
                   " 1. Exportar → HFSS/CST e simular 1 célula;",
                   " 2. aba Extração, modelo 'CRLH (CL-loaded)';",
                   " 3. comparar com LR, CR, LL, CL e corrigir dedos/stub."]
        self.saida.escrever("\n".join(linhas))

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
            self.salvar_s2p(f"crlh_cl_{g['n']}cel", f, s, z, "CRLH CL-loaded (circuito-alvo, Fig. 3.20)")
