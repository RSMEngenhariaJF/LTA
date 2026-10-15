"""Aba Aplicações EBG: linha de atraso chirpada (CDL) e ressonador de Bragg (seção 2.5.1)."""
import numpy as np

from ..core import abcd
from ..core import aplicacoes_ebg as ap
from ..core import modos_acoplados as cm
from ..geometria import construtores
from .aba_substrato import PREDEFINIDOS
from .comum import AbaBase, AreaTexto, Formulario, OpcoesModelo, protegido
from .previa import PreviaEstrutura

MODOS = ("Linha de atraso chirpada (CDL)", "Ressonador de Bragg")


class AbaAplicacoesEBG(AbaBase):
    TITULO = "Aplicações EBG (CDL, ressonador)"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.geral = Formulario(self.painel, [
            ("modo", "Aplicação", MODOS, ""),
            ("zref", "Z referência", 50.0, "Ω"),
            ("lalim", "Alimentação", 5.0, "mm"),
            ("fmin", "f mín", 1.0, "GHz"),
            ("fmax", "f máx", 17.0, "GHz"),
            ("npts", "Pontos", 1601, ""),
        ], titulo="Geral")
        self.geral.pack(fill="x")
        self.cdl = Formulario(self.painel, [
            ("f0", "f0 (central)", 9.0, "GHz"),
            ("df", "Δf (banda)", 12.0, "GHz"),
            ("psi", "ψ (inclinação de τg)", -0.5, "ns/GHz"),
            ("amp", "Amplitude A", 0.4, ""),
            ("nper", "Nº de períodos (0 = auto)", 0, ""),
            ("janela", "Janela W(z)", ap.JANELAS_CDL, ""),
        ], titulo="CDL (eqs. 2.109-2.115)")
        self.res = Formulario(self.painel, [
            ("f0", "f0 (ressonância)", 4.2, "GHz"),
            ("k", "|K| dos refletores", 0.64, "1/cm"),
            ("nper", "Períodos por refletor", 3, ""),
        ], titulo="Ressonador de Bragg (eq. 2.106)")
        self.modelo_fis = OpcoesModelo(self.painel)
        self.botoes_frame = self.botoes(("Calcular", self.calcular), ("Ex. CDL Fig. 2.26", self.ex_cdl),
                                        ("Ex. ressonador Fig. 2.21", self.ex_res),
                                        ("Enviar p/ exportação", self.publicar), ("Exportar S2P", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 16)
        self._organizar()
        self.geral.ao_mudar("modo", self._trocar_modo)
        self._ultimo = None
        self.observar_previa(self.geral, self.cdl, self.res)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    # ------------------------------------------------------------------ layout do painel
    def _organizar(self):
        for w in (self.cdl, self.res, self.modelo_fis, self.botoes_frame, self.saida):
            w.pack_forget()
        especifico = self.cdl if self._modo_cdl() else self.res
        especifico.pack(fill="x", pady=(8, 0))
        self.modelo_fis.pack(fill="x", pady=(8, 0))
        self.botoes_frame.pack(fill="x", pady=4)
        self.saida.pack(fill="both", expand=True)

    def _modo_cdl(self):
        return self.geral.get()["modo"] == MODOS[0]

    def _trocar_modo(self):
        if self._modo_cdl():
            self.geral.set(dict(fmin=1, fmax=17, npts=1601))
        else:
            self.geral.set(dict(fmin=3.5, fmax=5, npts=3001))
        self._organizar()
        self.atualizar_previa()

    def _ro3010(self):
        self.estado.definir_substrato(PREDEFINIDOS["Rogers RO3010 (livro)"])

    def ex_cdl(self):
        """Livro, p. 101: f0 = 9 GHz, Δf = 12 GHz, ψ = -0,5 ns/GHz, L = 50·l0, A = 0,4."""
        self._ro3010()
        self.geral.set(dict(modo=MODOS[0], zref=50, fmin=1, fmax=17, npts=1601))
        self.cdl.set(dict(f0=9, df=12, psi=-0.5, amp=0.4, nper=50, janela=ap.JANELAS_CDL[0]))
        self.modelo_fis.set(True, True)
        self._organizar()
        self.calcular()

    def ex_res(self):
        """Livro, p. 93: f0 = 4,2 GHz em RO3010 (lá com furos no terra; aqui largura modulada)."""
        self._ro3010()
        self.geral.set(dict(modo=MODOS[1], zref=50, fmin=3.5, fmax=5, npts=3001))
        self.res.set(dict(f0=4.2, k=0.64, nper=3))
        self.modelo_fis.set(True, True)
        self._organizar()
        self.calcular()

    # ------------------------------------------------------------------ projeto
    def _gerar(self):
        g = self.geral.get()
        sub = self.estado.substrato
        if self._modo_cdl():
            v = self.cdl.get()
            proj = ap.projeto_cdl(v["f0"] * 1e9, v["df"] * 1e9, v["psi"] * 1e-18, sub, g["zref"], v["amp"],
                                  v["nper"] or None)
            if proj.n_periodos > 2000:
                raise ValueError(f"{proj.n_periodos} períodos: estrutura longa demais. Reduza |ψ| ou Δf.")
            p = ap.perfil_cdl(proj, g["zref"], v["janela"])
            info = {"tipo": "CDL (linha de atraso chirpada)", "f0_GHz": v["f0"], "df_GHz": v["df"],
                    "psi_ns_GHz": f"{proj.psi * 1e18:.4f}", "l0_mm": f"{proj.l0 * 1e3:.4f}",
                    "C_1_m2": f"{proj.chirp:.1f}", "L_mm": f"{proj.comprimento * 1e3:.2f}",
                    "periodos": proj.n_periodos, "A": v["amp"], "janela": v["janela"]}
            nome = "cdl"
        else:
            v = self.res.get()
            if v["nper"] < 1:
                raise ValueError("Use ao menos 1 período por refletor.")
            perdas, dispersivo = self.modelo_fis.get()
            proj = ap.projeto_ressonador(v["f0"] * 1e9, v["k"] * 100, v["nper"], sub, g["zref"], perdas, dispersivo)
            p = proj.perfil
            info = {"tipo": "Ressonador de Bragg", "f0_GHz": v["f0"], "K_1_cm": v["k"],
                    "periodos_refletor": v["nper"], "periodo_mm": f"{p.periodo * 1e3:.4f}",
                    "Lm_mm": f"{proj.lm * 1e3:.4f}", "n_eq_2_106": proj.n}
            nome = "ressonador_bragg"
        w, ee = cm.larguras(p, sub)
        lay = construtores.linha_largura_modulada(p.z * 1e3, w, sub, g["lalim"], nome=nome, info=info)
        return g, v, proj, p, w, ee, lay

    def atualizar_previa(self):
        try:
            *_, lay = self._gerar()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        self.previa.mostrar(lay)

    @protegido
    def calcular(self):
        g, v, proj, p, w, ee, lay = self._gerar()
        self.previa.mostrar(lay)
        if g["fmin"] <= 0 or g["fmax"] <= g["fmin"] or g["npts"] < 50:
            raise ValueError("Faixa de frequência inválida (use ao menos 50 pontos).")
        f = np.linspace(g["fmin"], g["fmax"], g["npts"]) * 1e9
        sub = self.estado.substrato
        perdas, dispersivo = self.modelo_fis.get()
        z0, eef, al = cm.parametros_linha(f, p, sub, perdas, dispersivo)
        s11, s21 = cm.resolver_cmt(f, p.z, z0, eef, al)
        s22, _ = cm.resolver_cmt(f, p.z, z0[:, ::-1], eef[:, ::-1], al[:, ::-1])
        s = np.empty((f.size, 2, 2), dtype=complex)
        s[:, 0, 0], s[:, 1, 0], s[:, 0, 1], s[:, 1, 1] = s11, s21, s21, s22
        self._ultimo = (f, s, g["zref"], lay.nome)
        self.estado.publicar_resultado(lay.nome, f, s, g["zref"])
        if self._modo_cdl():
            texto = self._graficos_cdl(f, s11, s21, proj, p, w, v)
        else:
            texto = self._graficos_res(f, s11, s21, proj, p, w)
        self.grafico.desenhar()
        self.ver_resultados()
        self.saida.escrever(texto + "\n\n" + self.modelo_fis.descrever(sub)
                            + "\nCMT: despreza descontinuidades e radiação.")

    def _perfil(self, ax, p, w):
        zmm = p.z * 1e3
        ax.plot(zmm, p.z0)
        ax.set_xlabel("z [mm]")
        ax.set_ylabel("Z0 [Ω]")
        axw = ax.twinx()
        axw.plot(zmm, w, color="C1", lw=0.8)
        axw.set_ylabel("W [mm]", color="C1")
        ax.set_title("Perfil de impedância e largura", fontsize=9)

    def _graficos_cdl(self, f, s11, s21, proj, p, w, v):
        ax = self.grafico.eixos(2, 2)
        self._perfil(ax[0], p, w)
        zc = p.z - proj.comprimento / 2
        ax[1].plot(p.z * 1e3, proj.f_local(zc) / 1e9)
        ax[1].set_xlabel("z [mm]")
        ax[1].set_ylabel("f refletida [GHz]")
        ax[1].set_title("Frequência refletida em cada posição (eq. 2.111)", fontsize=9)
        fg = f / 1e9
        ax[2].plot(fg, abcd.db(s11), label="|S11|")
        ax[2].plot(fg, abcd.db(s21), label="|S21|")
        ax[2].set_ylim(-40, 2)
        ax[2].set_ylabel("|S| [dB]")
        ax[2].set_title("Resposta (CMT)", fontsize=9)
        tau = abcd.atraso_grupo(f, s11)
        f_ini, f_fim = proj.f0 - 0.4 * proj.df, proj.f0 + 0.4 * proj.df
        psi, tc = ap.inclinacao_atraso(f, tau, f_ini, f_fim)
        banda = (f >= proj.f0 - proj.df / 2) & (f <= proj.f0 + proj.df / 2)
        ax[3].plot(fg, tau * 1e9, lw=0.8, label="τg de S11")
        ax[3].plot(fg[banda], (tc + proj.psi * (f[banda] - proj.f0)) * 1e9, "k--", label="ψ de projeto")
        ax[3].set_ylabel("τg [ns]")
        ax[3].set_ylim(min(0, np.nanmin(tau[banda]) * 1e9 - 0.5), np.nanmax(tau[banda]) * 1e9 + 1)
        ax[3].set_title("Atraso de grupo em reflexão", fontsize=9)
        for a in ax[2:]:
            a.set_xlabel("f [GHz]")
            a.legend(fontsize=7)
        for a in ax:
            a.grid(True, alpha=0.3)
        s11b = abcd.db(s11)[banda]
        return "\n".join([
            "CDL (eqs. 2.109-2.115):",
            f"  l0 = {proj.l0 * 1e3:.3f} mm   C = {proj.chirp:.1f} 1/m²",
            f"  L = {proj.n_periodos}·l0 = {proj.comprimento * 1e3:.2f} mm   εeff = {proj.eeff:.3f}",
            f"  ψ de projeto = {proj.psi * 1e18:.4f} ns/GHz"
            + (" (L arredondado para N·l0)" if abs(proj.psi * 1e18 - v["psi"]) > 1e-6 else ""),
            f"  Faixa refletida: {proj.f_local(proj.comprimento / 2) / 1e9:.2f} a "
            f"{proj.f_local(-proj.comprimento / 2) / 1e9:.2f} GHz", "",
            f"Simulado (ajuste em {f_ini / 1e9:.2f}-{f_fim / 1e9:.2f} GHz):",
            f"  ψ = {psi * 1e18:.4f} ns/GHz ({100 * (psi / proj.psi - 1):+.1f}%)",
            f"  τg em f0 = {tc * 1e9:.3f} ns",
            f"  |S11| na banda: média {s11b.mean():.2f} dB, de {s11b.min():.2f} a {s11b.max():.2f} dB",
            "", "ψ < 0: altas frequências refletem no início da linha,",
            "o que equaliza as perdas (livro, p. 101)."]
            + (["", "Dispersão ligada: as eqs. 2.111-2.115 usam εeff fixo;",
                "como εeff cresce com f, a banda real se comprime e ψ",
                "muda. Desligue a dispersão para conferir o projeto."]
               if self.modelo_fis.get()[1] else []))

    def _graficos_res(self, f, s11, s21, proj, p, w):
        ax = self.grafico.eixos(2, 2)
        self._perfil(ax[0], p, w)
        fg = f / 1e9
        ax[1].plot(fg, abcd.db(s21), label="|S21|")
        ax[1].plot(fg, abcd.db(s11), label="|S11|")
        ax[1].axvline(proj.f0 / 1e9, color="gray", ls=":", lw=0.8)
        ax[1].set_ylabel("|S| [dB]")
        ax[1].set_title("Resposta do ressonador (CMT)", fontsize=9)
        try:
            fp, il, ql, q0 = ap.fator_q(f, s21)
            meia = max(3 * fp / ql, 5 * (f[1] - f[0]))
            perto = np.abs(f - fp) < meia
            ax[2].plot(fg[perto], abcd.db(s21[perto]))
            ax[2].axhline(-il - 3, color="gray", ls=":", lw=0.8)
            ax[2].set_title(f"Ressonância: QL = {ql:.1f}", fontsize=9)
            qtxt = [f"  Pico em {fp / 1e9:.4f} GHz ({100 * (fp / proj.f0 - 1):+.2f}% de f0)",
                    f"  IL = {il:.2f} dB   QL = {ql:.1f} (eq. 1.75)   "
                    f"Q0 = {'∞ (sem perdas)' if np.isinf(q0) else f'{q0:.1f}'} (eq. 1.74)"]
        except ValueError as e:
            ax[2].set_title("Ressonância não resolvida", fontsize=9)
            qtxt = [f"  {e}"]
        ax[2].set_ylabel("|S21| [dB]")
        # fases dos refletores (vistas da cavidade)
        ref = proj.refletor
        sub = self.estado.substrato
        perdas, dispersivo = self.modelo_fis.get()
        z0, eef, al = cm.parametros_linha(f, ref, sub, perdas, dispersivo)
        r11, _ = cm.resolver_cmt(f, ref.z, z0, eef, al)
        r22, _ = cm.resolver_cmt(f, ref.z, z0[:, ::-1], eef[:, ::-1], al[:, ::-1])
        ax[3].plot(fg, np.degrees(np.angle(r22)), label="φref1 (S22 refletor 1)")
        ax[3].plot(fg, np.degrees(np.angle(r11)), label="φref2 (S11 refletor 2)")
        ax[3].axvline(proj.f0 / 1e9, color="gray", ls=":", lw=0.8)
        ax[3].set_ylabel("Fase [°]")
        ax[3].set_title("Fase da reflexão dos refletores", fontsize=9)
        for a in ax[1:]:
            a.set_xlabel("f [GHz]")
        ax[1].legend(fontsize=7)
        ax[3].legend(fontsize=7)
        for a in ax:
            a.grid(True, alpha=0.3)
        return "\n".join([
            "Ressonador de Bragg (eq. 2.106):",
            f"  Refletores: senoidais, l = {ref.periodo * 1e3:.3f} mm (Bragg em f0), {ref.n_periodos} períodos",
            f"  φref1 = {np.degrees(proj.phi1):.2f}°   φref2 = {np.degrees(proj.phi2):.2f}°",
            f"  β na cavidade = {proj.beta:.2f} rad/m",
            f"  n = {proj.n}  →  βLm = {np.degrees(proj.beta * proj.lm):.2f}°,  Lm = {proj.lm * 1e3:.3f} mm",
            "", "Simulado:"] + qtxt + [
            "", "Livro (furos no terra): φref = -95,84°, βLm = 84,15°,",
            "Lm = 6,25 mm, Q medido = 129,2."])

    # ------------------------------------------------------------------ exportação
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
            f, s, z, nome = self._ultimo
            self.salvar_s2p(nome, f, s, z, f"{nome} (CMT)")
