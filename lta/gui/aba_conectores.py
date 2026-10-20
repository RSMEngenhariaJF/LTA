"""Aba Conectores SMA e transições coaxial-microstrip (Garg, Bahl e Bianchi, seção 1.4.1)."""
from pathlib import Path
from tkinter import filedialog

import numpy as np

from ..core import abcd
from ..core import conectores as co
from ..core.microstrip import sintese
from ..exportacao.simuladores import ler_touchstone
from ..geometria import circuitos
from .comum import PASTA_SAIDAS, AbaBase, AreaTexto, Formulario, protegido, relativo
from .previa import PreviaEstrutura

FONTES = ("Predefinido", "Cs e Ls informados", "Cs pelo VSWR medido")
COMPENSAR = ("Sim", "Não")


class AbaConectores(AbaBase):
    TITULO = "Conectores SMA e transições"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("fonte", "Modelo da transição", FONTES, ""),
            ("predef", "Predefinido", tuple(co.PREDEFINIDOS), ""),
            ("cs", "Cs (derivação)", 30.0, "fF"),
            ("ls", "Ls (série)", 0.0, "pH"),
            ("vswr", "VSWR medido de uma transição", 1.10, ""),
            ("fv", "… na frequência", 10.0, "GHz"),
            ("z0", "Z0", 50.0, "Ω"),
            ("comp", "Linha de teste entre conectores", 30.0, "mm"),
            ("comp_on", "Trecho de compensação", COMPENSAR, ""),
            ("zc", "Z do trecho de compensação (0 = auto)", 0.0, "Ω"),
            ("alim", "De-embedding: alimentação até o DUT", 5.0, "mm"),
            ("fmax", "f máx", 12.0, "GHz"),
        ], titulo="Transição coaxial-microstrip (Fig. 1.20)")
        self.form.pack(fill="x")
        self.botoes(("Calcular", self.calcular), ("Enviar p/ exportação", self.publicar),
                    ("De-embedding de .s2p…", self.de_embedding))
        self.saida = AreaTexto(self.painel, 22)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Transição modelada como Ls em série e Cs em derivação (Fig. 1.20).\n"
            "Boas práticas (seção 1.4.1):\n"
            "• pino central soldado sobre a trilha;\n"
            "• sem folga entre o terra da placa e a base do fixture, nem entre o\n"
            "  flange do conector e a parede (a corrente de terra faria um laço\n"
            "  longo e somaria indutância);\n"
            "• minimizar Ls e Cs (pino deslocado/coaxial afilado) ou compensar\n"
            "  com √(Ls/Cs) = Z0; aqui, um trecho curto de alta ou baixa\n"
            "  impedância junto ao pad cancela a reatância líquida.\n"
            "Os valores de Ls e Cs dependem do conector e do pad: obtenha-os em\n"
            "EM ou medindo (aba Caracterização) antes de confiar no resultado.")
        self._ultimo = None
        self.observar_previa(self.form)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    def _transicao(self, v):
        if v["fonte"] == FONTES[0]:
            return co.PREDEFINIDOS[v["predef"]]
        if v["fonte"] == FONTES[1]:
            return co.Transicao("informada", v["ls"] * 1e-12, v["cs"] * 1e-15)
        return co.Transicao("vswr", 0.0, co.cs_de_vswr(v["vswr"], v["fv"] * 1e9, v["z0"]))

    def _projeto(self):
        v = self.form.get()
        sub = self.estado.substrato
        t = self._transicao(v)
        trecho = None
        if v["comp_on"] == COMPENSAR[0]:
            trecho = co.compensacao(t, sub, v["z0"], v["zc"] or None)
        w0 = float(sintese(v["z0"], sub))
        info = {"tipo": "fixture com conectores", "Ls_pH": f"{t.ls * 1e12:.2f}", "Cs_fF": f"{t.cs * 1e15:.2f}",
                "Z0_ohm": v["z0"]}
        if trecho:
            info["compensacao_Z_W_l"] = f"{trecho[0]:.1f} ohm / {trecho[1]:.4f} mm / {trecho[2]:.4f} mm"
        lay = circuitos.fixture_conector("fixture_conectores", sub, w0, v["comp"],
                                         (trecho[1], trecho[2]) if trecho else None, info)
        return v, sub, t, trecho, lay

    def atualizar_previa(self):
        try:
            *_, lay = self._projeto()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        self.previa.mostrar(lay)

    @protegido
    def calcular(self):
        v, sub, t, trecho, lay = self._projeto()
        self.previa.mostrar(lay)
        f = np.linspace(v["fmax"] * 1e9 / 600, v["fmax"] * 1e9, 600)
        s_t = abcd.abcd_para_s(co.abcd_transicao(f, t), v["z0"])
        s_sem = abcd.abcd_para_s(co.dispositivo_de_teste(f, t, sub, v["comp"], v["z0"]), v["z0"])
        s_com = abcd.abcd_para_s(co.dispositivo_de_teste(f, t, sub, v["comp"], v["z0"], comp_trecho=trecho),
                                 v["z0"]) if trecho else None
        self._ultimo = (f, s_com if s_com is not None else s_sem, v["z0"])
        self.estado.publicar_resultado("fixture_conectores", f, self._ultimo[1], v["z0"])
        fg = f / 1e9
        ax = self.grafico.eixos(1, 2)
        ax[0].plot(fg, abcd.db(s_t[:, 0, 0]), ls="--", label="1 transição")
        ax[0].plot(fg, abcd.db(s_sem[:, 0, 0]), label="fixture sem compensação")
        if s_com is not None:
            ax[0].plot(fg, abcd.db(s_com[:, 0, 0]), label="fixture compensado")
        ax[0].set_ylim(-60, 0)
        ax[0].set_ylabel("|S11| [dB]")
        ax[0].set_title("Perda de retorno", fontsize=9)
        ax[1].plot(fg, co.vswr(s_t[:, 0, 0]), ls="--", label="1 transição")
        ax[1].plot(fg, co.vswr(s_sem[:, 0, 0]), label="fixture sem compensação")
        if s_com is not None:
            ax[1].plot(fg, co.vswr(s_com[:, 0, 0]), label="fixture compensado")
        ax[1].axhline(1.065, color="gray", ls=":", lw=0.8)
        ax[1].axhline(1.10, color="gray", ls=":", lw=0.8)
        ax[1].set_ylim(1, 1.6)
        ax[1].set_ylabel("VSWR")
        ax[1].set_title("VSWR (linhas: 1,065 e 1,10 da Tabela 2.5)", fontsize=9)
        for a in ax:
            a.set_xlabel("f [GHz]")
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        i = -1
        linhas = [f"Transição: Ls = {t.ls * 1e12:.2f} pH, Cs = {t.cs * 1e15:.2f} fF "
                  f"(excesso capacitivo líquido {(t.cs - t.ls / v['z0'] ** 2) * 1e15:+.2f} fF)",
                  f"√(Ls/Cs) = {np.sqrt(t.ls / t.cs) if t.cs > 0 and t.ls > 0 else 0:.1f} Ω (compensada se = Z0)", "",
                  f"Em {fg[i]:g} GHz: VSWR de uma transição = {co.vswr(s_t[i, 0, 0]):.3f}; fixture: "
                  f"{co.vswr(s_sem[i, 0, 0]):.3f}" + (f" → {co.vswr(s_com[i, 0, 0]):.3f} compensado" if s_com is not None
                                                      else "")]
        if trecho:
            linhas.append(f"Compensação: trecho de {trecho[0]:.1f} Ω (W = {trecho[1]:.4f} mm) com {trecho[2]:.4f} mm "
                          "junto a cada pad")
        elif v["comp_on"] == COMPENSAR[0]:
            linhas.append("Transição já compensada (√(Ls/Cs) = Z0): nenhum trecho necessário.")
        linhas += ["", "Perda de inserção do fixture em f máx: "
                       f"{-abcd.db(self._ultimo[1][i, 1, 0]):.3f} dB (sem perdas na linha)",
                   "Use 'De-embedding' para remover transições e alimentações de uma medida .s2p."]
        self.saida.escrever("\n".join(linhas))

    @protegido
    def publicar(self):
        *_, lay = self._projeto()
        self.previa.mostrar(lay)
        self.estado.publicar_layout(lay)
        self.estado.msg(f"Layout '{lay.nome}' enviado para Exportar (Ctrl+E).")

    @protegido
    def de_embedding(self):
        caminho = filedialog.askopenfilename(initialdir=PASTA_SAIDAS, filetypes=[("Touchstone 2 portas", "*.s2p")])
        if not caminho:
            return
        v, sub, t, trecho, _ = self._projeto()
        f, s, z_ref = ler_touchstone(caminho)
        s_dut = co.de_embutir(f, s, t, sub, v["alim"], z_ref, comp_trecho=(None, trecho[1], trecho[2])
                              if trecho else None)
        fg = f / 1e9
        ax = self.grafico.eixos(1, 2)
        for j, nome in ((0, "S11"), (1, "S21")):
            ax[j].plot(fg, abcd.db(s[:, j, 0]), ls="--", label=f"{nome} medido")
            ax[j].plot(fg, abcd.db(s_dut[:, j, 0]), label=f"{nome} do DUT (de-embedded)")
            ax[j].set_xlabel("f [GHz]")
            ax[j].set_ylabel("[dB]")
            ax[j].grid(True, alpha=0.3)
            ax[j].legend(fontsize=7)
        self.grafico.desenhar()
        PASTA_SAIDAS.mkdir(exist_ok=True)
        destino = Path(caminho).with_name(Path(caminho).stem + "_deembedded.s2p")
        from ..exportacao.simuladores import touchstone
        touchstone(destino, f, s_dut, z_ref, "de-embedding das transicoes coaxial-microstrip (LTA)")
        self.saida.escrever(f"De-embedding de {relativo(caminho)}:\n  removidas 2 transições (Ls = {t.ls * 1e12:.2f} pH, "
                            f"Cs = {t.cs * 1e15:.2f} fF) e {v['alim']:g} mm de alimentação por lado.\n"
                            f"  Resultado salvo em {relativo(destino)}.")
