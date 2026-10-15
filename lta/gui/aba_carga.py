"""Aba Linha com carga capacitiva (slow-wave), seção 2.5.2."""
import numpy as np

from ..core import abcd
from ..core import carga_capacitiva as cc
from ..geometria import construtores
from .aba_substrato import PREDEFINIDOS
from .comum import AbaBase, AreaTexto, Formulario, protegido
from .previa import PreviaEstrutura

MODOS = ("fc, ZB, fs", "Fase/célula, ZB, fs")


class AbaCarga(AbaBase):
    TITULO = "Carga capacitiva"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("modo", "Especificação", MODOS, ""),
            ("zb", "ZB (Bloch)", 50.0, "Ω"),
            ("fc", "fc (corte)", 5.0, "GHz"),
            ("fase", "Fase por célula", 35.0, "°"),
            ("f0", "em f0", 1.4, "GHz"),
            ("fs", "fs (1ª espúria)", 20.0, "GHz"),
            ("n", "Nº de células", 4, ""),
        ], titulo="Especificação elétrica")
        self.form.pack(fill="x")
        self.fis = Formulario(self.painel, [
            ("zpatch", "Z do patch", 20.0, "Ω"),
            ("zporta", "Z das portas", 50.0, "Ω"),
            ("lalim", "Alimentação", 5.0, "mm"),
            ("fmax", "f máx (gráficos)", 25.0, "GHz"),
        ], titulo="Implementação microstrip")
        self.fis.pack(fill="x", pady=(8, 0))
        self.botoes(("Calcular", self.calcular), ("Exemplo Fig. 2.34", self.exemplo),
                    ("Enviar p/ exportação", self.publicar), ("Exportar S2P", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 16)
        self.saida.pack(fill="both", expand=True)
        self._ultimo = None
        self._lay = None
        self.observar_previa(self.form, self.fis)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    def exemplo(self):
        """Filtro UWB do livro: 35° a 1,4 GHz, ZB = 47,5 Ω, fs = 22 GHz, substrato Arlon."""
        self.estado.definir_substrato(PREDEFINIDOS["Arlon er=2.4 h=0.675 (livro, Fig. 2.33)"])
        self.form.set(dict(modo=MODOS[1], zb=47.5, fase=35, f0=1.4, fs=22, n=2))
        self.fis.set(dict(zpatch=10, zporta=47.5, lalim=3, fmax=25))
        self.calcular()
        self.estado.msg("Exemplo carregado (substrato alterado para Arlon εr=2,4, h=0,675 mm).")

    def _gerar(self):
        """Projeto elétrico + geometria + layout a partir dos campos."""
        v, fv = self.form.get(), self.fis.get()
        if v["n"] < 1:
            raise ValueError("Número de células deve ser ≥ 1.")
        if v["modo"] == MODOS[0]:
            p = cc.projeto_por_fc(v["fc"] * 1e9, v["zb"], v["fs"] * 1e9)
        else:
            p = cc.projeto_por_fase(v["fase"], v["f0"] * 1e9, v["zb"], v["fs"] * 1e9)
        geo = cc.geometria(p, self.estado.substrato, fv["zpatch"], fv["zporta"])
        info = {"tipo": "linha com carga capacitiva", "L_nH": f"{p.L * 1e9:.4f}",
                "C_pF": f"{p.C * 1e12:.4f}", "Cls_pF": f"{p.Cls * 1e12:.4f}", "ZB_ohm": f"{p.zb:.2f}",
                "fc_GHz": f"{p.fc / 1e9:.3f}", "fs_GHz": f"{p.fs / 1e9:.3f}", "celulas": v["n"],
                "periodo_mm": f"{geo.periodo:.4f}", "W_linha_mm": f"{geo.w_linha:.4f}",
                "W_patch_mm": f"{geo.w_patch:.4f}", "l_patch_mm": f"{geo.l_patch:.4f}"}
        lay = construtores.linha_carga_capacitiva(geo, v["n"], self.estado.substrato, fv["lalim"],
                                                  nome=f"carga_cap_{v['n']}cel", info=info)
        return v, fv, p, geo, lay

    def atualizar_previa(self):
        try:
            *_, lay = self._gerar()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        self._lay = lay
        self.previa.mostrar(lay)

    @protegido
    def calcular(self):
        v, fv, p, geo, lay = self._gerar()
        self._lay = lay
        self.previa.mostrar(lay)
        f = np.linspace(0.01, fv["fmax"], 1500) * 1e9

        md = cc.celula_distribuida(f, p)
        mc = cc.celula_concentrada(f, p)
        _, bl_d = abcd.dispersao(md)
        _, bl_c = abcd.dispersao(mc)
        zb_d, _ = abcd.impedancia_bloch(md)
        zb_c, _ = abcd.impedancia_bloch(mc)
        s = abcd.abcd_para_s(abcd.potencia(md, v["n"]), fv["zporta"])
        self._ultimo = (f, s, fv["zporta"], v, p)
        self.estado.publicar_resultado("carga_capacitiva", f, s, fv["zporta"])

        ax = self.grafico.eixos(2, 2)
        ax[0].plot(bl_d / np.pi, f / 1e9, label="Distribuído (2.116)")
        ax[0].plot(bl_c / np.pi, f / 1e9, "--", label="Concentrado (Fig. 2.27b)")
        ax[0].plot(np.clip(2 * f * np.sqrt(p.L * p.C), 0, None), f / 1e9, ":", color="gray",
                   label="Linha sem carga")
        ax[0].set_xlim(0, 1)
        ax[0].set_xlabel("βl/π")
        ax[0].set_ylabel("f [GHz]")
        ax[0].set_title("Diagrama de Brillouin", fontsize=9)
        ax[1].plot(f / 1e9, zb_d.real, label="Distribuído")
        ax[1].plot(f / 1e9, zb_c.real, "--", label="Concentrado")
        ax[1].set_ylim(0, 2 * v["zb"])
        ax[1].set_ylabel("Re ZB [Ω]")
        ax[1].set_title("Impedância de Bloch", fontsize=9)
        ax[2].plot(f / 1e9, abcd.db(s[:, 1, 0]), label="S21")
        ax[2].plot(f / 1e9, abcd.db(s[:, 0, 0]), label="S11")
        ax[2].set_ylim(-60, 2)
        ax[2].set_ylabel("|S| [dB]")
        ax[2].set_title(f"{v['n']} células (modelo distribuído)", fontsize=9)
        passa = f < p.fc
        ax[3].plot(f[passa] / 1e9, np.rad2deg(bl_d[passa]), label="Distribuído")
        ax[3].plot(f[passa] / 1e9, np.rad2deg(bl_c[passa]), "--", label="Concentrado")
        ax[3].set_ylabel("Fase por célula [°]")
        ax[3].set_title("Defasagem por célula (1ª banda)", fontsize=9)
        for a in ax:
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        for a in ax[1:]:
            a.set_xlabel("f [GHz]")
        self.grafico.desenhar()
        self.ver_resultados()

        self.saida.escrever(
            f"{abcd.descrever_banda(f, s)}\n(regra: ZB estável até ~0,55·fc = {0.55 * p.fc / 1e9:.2f} GHz)\n\n"
            "Elementos por célula:\n"
            f"  L   = {p.L * 1e9:.4f} nH\n  C   = {p.C * 1e12:.4f} pF\n  Cls = {p.Cls * 1e12:.4f} pF\n"
            f"  Zo (linha sem carga) = {p.zo:.2f} Ω\n"
            f"  ZB (baixa freq.)     = {p.zb:.2f} Ω\n"
            f"  fc = {p.fc / 1e9:.3f} GHz   fs = {p.fs / 1e9:.3f} GHz\n"
            f"  Fator slow-wave      = {p.fator_slow_wave:.2f}x\n"
            f"  Fase/célula em f0    = {np.rad2deg(p.fase_por_celula(v['f0'] * 1e9)):.2f}°\n\n"
            f"Geometria inicial ({self.estado.substrato.nome}):\n"
            f"  Período l   = {geo.periodo:.3f} mm\n"
            f"  W linha     = {geo.w_linha:.4f} mm (εeff {geo.eeff_linha:.3f})\n"
            f"  Patch W x l = {geo.w_patch:.3f} x {geo.l_patch:.3f} mm\n"
            f"  W portas    = {geo.w_porta:.4f} mm\n"
            f"  Comprimento total = {lay.comprimento:.2f} mm\n\n"
            "Geometria é estimativa quase-estática (sem efeito de\n"
            "degrau/borda): otimizar em simulador EM.")

    @protegido
    def publicar(self):
        *_, lay = self._gerar()
        self._lay = lay
        self.previa.mostrar(lay)
        self.estado.publicar_layout(lay)
        self.estado.msg(f"Layout '{lay.nome}' enviado para Exportar (Ctrl+E).")

    def exportar_s2p(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s, z, v, p = self._ultimo
            self.salvar_s2p(f"carga_cap_{v['n']}cel", f, s, z,
                            f"Linha com carga capacitiva: L={p.L * 1e9:.4f}nH C={p.C * 1e12:.4f}pF "
                            f"Cls={p.Cls * 1e12:.4f}pF N={v['n']}")
