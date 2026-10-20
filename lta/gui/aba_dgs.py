"""Aba DGS: célula em haltere, linha periódica e passa-baixas (Garg, Bahl e Bianchi, cap. 6)."""
from tkinter import filedialog

import numpy as np

from ..core import abcd
from ..core import dgs
from ..core.microstrip import C0, analise, sintese
from ..exportacao.simuladores import ler_touchstone
from ..geometria import circuitos
from .comum import PASTA_SAIDAS, AbaBase, AreaTexto, Formulario, protegido, relativo
from .previa import PreviaEstrutura

MODOS = ("Célula: LC a partir de f0 e fc (eq. 6.8)", "Linha periódica com DGS (rejeição e onda lenta)",
         "Passa-baixas com DGS e stubs (seção 6.4.1)")


class AbaDGS(AbaBase):
    TITULO = "DGS: terra com defeitos (Garg-Bahl, cap. 6)"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("modo", "Projeto", MODOS, ""),
            ("z0", "Z0 da linha", 50.0, "Ω"),
            ("f0", "Rejeição f0 da célula (EM)", 8.0, "GHz"),
            ("fc", "Corte de 3 dB fc da célula (EM)", 3.5, "GHz"),
            ("n", "Nº de células", 3, ""),
            ("p", "Período", 10.0, "mm"),
            ("fcf", "Passa-baixas: corte", 1.3, "GHz"),
            ("ripple", "Passa-baixas: ripple", 0.01, "dB"),
            ("a", "Haltere: cabeça a (ao longo da linha)", 5.0, "mm"),
            ("b", "Haltere: cabeça b (transversal)", 5.0, "mm"),
            ("g", "Haltere: largura do pescoço g", 0.5, "mm"),
            ("lp", "Haltere: comprimento do pescoço", 2.4, "mm"),
            ("fmax", "f máx", 12.0, "GHz"),
        ], titulo="DGS em haltere")
        self.form.pack(fill="x")
        self.botoes(("Calcular", self.calcular), ("Enviar p/ exportação", self.publicar),
                    ("Ex. Fig. 6.1 (8 / 3,5 GHz)", self.ex_celula), ("Ex. passa-baixas 1,3 GHz", self.ex_lpf),
                    ("Extrair modelo π de .s2p…", self.extrair))
        self.saida = AreaTexto(self.painel, 22)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "A DGS em haltere é um LC paralelo em série com a linha: C e L saem\n"
            "da rejeição f0 e do corte fc obtidos em EM (eq. 6.8). As dimensões\n"
            "a, b, g e pescoço só desenham o layout: o efeito delas (L cresce\n"
            "com as cabeças, C com o pescoço estreito) deve ser ajustado em EM\n"
            "e confirmado pela extração do modelo π (eqs. 6.10-6.13).\n"
            "A linha periódica é uma linha artificial: banda de rejeição de\n"
            "Bragg/ressonância e fator de onda lenta abaixo dela.")
        self._ultimo = None
        self.observar_previa(self.form)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    def ex_celula(self):
        from ..core.microstrip import Substrato
        self.estado.definir_substrato(Substrato(2.2, 0.787, 35, 0.0009, "Fig. 6.1: er 2,2, h 31 mil"))
        self.form.set(dict(modo=MODOS[0], z0=50, f0=8, fc=3.5, n=1, a=5, b=5, g=0.5, lp=2.4, fmax=12))
        self.calcular()

    def ex_lpf(self):
        self.form.set(dict(modo=MODOS[2], z0=50, fcf=1.3, ripple=0.01, f0=4.0, fmax=6))
        self.calcular()

    def _projeto(self):
        v = self.form.get()
        sub = self.estado.substrato
        w = float(sintese(v["z0"], sub))
        info = {"tipo": f"DGS - {v['modo']}", "Z0_ohm": v["z0"]}
        L, C = dgs.lc_de_resposta(v["f0"] * 1e9, v["fc"] * 1e9, v["z0"])
        p = dict(v=v, sub=sub, w=w, L=L, C=C, texto=[])
        if v["modo"] == MODOS[2]:
            el = dgs.passa_baixas(v["fcf"] * 1e9, v["z0"], 3, v["ripple"])
            p["el"] = el
            stubs, desenho = {}, []
            xs = [v["p"] * (k + 0.5) for k in range(3)]
            for k, (tipo, val) in enumerate(el):
                if tipo == "C":
                    comp, wst = dgs.stub_para_capacitor(val, v["fcf"] * 1e9, sub, v["z0"], True)
                    stubs[k] = (comp, wst, True)
                    desenho.append((xs[k], comp, wst))
            p["stubs"] = stubs
            p["texto"] += [f"Protótipo Chebyshev N = 3, {v['ripple']:g} dB, fc = {v['fcf']:g} GHz:"] + \
                [f"  {t}{k + 1} = " + (f"{x * 1e9:.4f} nH (DGS)" if t == "L" else f"{x * 1e12:.4f} pF "
                                       f"(par de stubs de {stubs[k][0]:.3f} mm, W = {stubs[k][1]:.3f} mm)")
                 for k, (t, x) in enumerate(el)]
            cd = 1 / ((2 * np.pi * v["f0"] * 1e9) ** 2 * el[0][1])
            p["texto"].append(f"  Rejeição das DGS em {v['f0']:g} GHz: Cd = {cd * 1e12:.4f} pF")
            lay = circuitos.linha_dgs("passa_baixas_dgs", sub, w, 3, v["p"], v["a"], v["b"], v["g"], v["lp"],
                                      stubs=desenho, info=info)
            lay.fendas_terra = lay.fendas_terra[:3] + lay.fendas_terra[6:]   # DGS só nos indutores (1 e 3)
            p["lay"] = lay
        else:
            n = 1 if v["modo"] == MODOS[0] else max(1, v["n"])
            info.update(L_nH=f"{L * 1e9:.4f}", C_pF=f"{C * 1e12:.4f}")
            p["lay"] = circuitos.linha_dgs("linha_dgs" if n > 1 else "celula_dgs", sub, w, n, v["p"], v["a"], v["b"],
                                           v["g"], v["lp"], info=info)
            p["n"] = n
        return p

    def atualizar_previa(self):
        try:
            p = self._projeto()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        self.previa.mostrar(p["lay"])

    @protegido
    def calcular(self):
        p = self._projeto()
        v, sub = p["v"], p["sub"]
        self.previa.mostrar(p["lay"])
        f = np.linspace(v["fmax"] * 1e9 / 800, v["fmax"] * 1e9, 800)
        fg = f / 1e9
        cel = dgs.CelulaDGS(p["L"], p["C"])
        ax = self.grafico.eixos(1, 2)
        linhas = list(p["texto"])
        if v["modo"] == MODOS[2]:
            el = p["el"]
            s_id = abcd.abcd_para_s(dgs.abcd_passa_baixas(f, el, v["z0"]), v["z0"])
            s_dgs = abcd.abcd_para_s(dgs.abcd_passa_baixas(f, el, v["z0"], sub, p["stubs"], v["fcf"] * 1e9,
                                                           v["f0"] * 1e9), v["z0"])
            ax[0].plot(fg, abcd.db(s_id[:, 1, 0]), ls="--", label="|S21| concentrado")
            ax[0].plot(fg, abcd.db(s_dgs[:, 1, 0]), label="|S21| DGS + stubs")
            ax[0].plot(fg, abcd.db(s_dgs[:, 0, 0]), lw=0.8, label="|S11| DGS + stubs")
            ax[1].plot(fg, np.degrees(np.angle(s_dgs[:, 1, 0])))
            ax[1].set_title("Fase de S21", fontsize=9)
            self._ultimo = (f, s_dgs)
            k = int(np.argmin(np.abs(f - v["fcf"] * 1e9)))
            linhas += ["", f"|S21| em fc: {abcd.db(s_dgs[k, 1, 0]):.2f} dB; a DGS também dá um zero em "
                           f"{v['f0']:g} GHz, que melhora a rejeição (seção 6.4.1)."]
        else:
            n = p["n"]
            m = dgs.celula_periodica(f, cel, v["p"], p["w"], sub)
            s = abcd.abcd_para_s(abcd.potencia(m, n), v["z0"])
            ax[0].plot(fg, abcd.db(s[:, 1, 0]), label=f"|S21|, {n} célula(s)")
            ax[0].plot(fg, abcd.db(s[:, 0, 0]), lw=0.8, label="|S11|")
            al, bl = abcd.dispersao(m)
            k0l = 2 * np.pi * f / C0 * v["p"] * 1e-3
            ax[1].plot(fg, np.degrees(bl), label="βl da célula (Bloch)")
            ax[1].plot(fg, np.degrees(k0l * np.sqrt(sub.er)), ls=":", lw=0.8, label="k0·l·√εr")
            ax[1].set_ylabel("[°]")
            ax[1].set_title("Dispersão da linha periódica", fontsize=9)
            self._ultimo = (f, s)
            i = len(f) // 20
            swf = bl[i] / (k0l[i]) if k0l[i] > 0 else np.nan
            linhas += [f"Célula: L = {p['L'] * 1e9:.4f} nH, C = {p['C'] * 1e12:.4f} pF (eq. 6.8), ressonância "
                       f"{cel.f0 / 1e9:.3f} GHz",
                       f"Linha periódica: fator de onda lenta βl/(k0·l) ≈ {swf:.3f} em {fg[i]:.2f} GHz "
                       f"(√εeff da linha sem DGS ≈ {np.sqrt(float(analise(p['w'], sub)[1])):.3f})"]
            rej = abcd.db(s[:, 1, 0]) < -20
            if rej.any():
                linhas.append(f"|S21| < -20 dB de {fg[rej].min():.3f} a {fg[rej].max():.3f} GHz")
        ax[0].set_ylim(-60, 2)
        ax[0].set_ylabel("|S| [dB]")
        ax[0].set_title(v["modo"], fontsize=9)
        for a in ax:
            a.set_xlabel("f [GHz]")
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        linhas += ["", "Geometria: estimativa a ajustar em EM; extraia L e C do .s2p simulado com o botão"
                       " 'Extrair modelo π'."]
        self.saida.escrever("\n".join(linhas))

    @protegido
    def publicar(self):
        p = self._projeto()
        self.previa.mostrar(p["lay"])
        self.estado.publicar_layout(p["lay"])
        self.estado.msg(f"Layout '{p['lay'].nome}' enviado para Exportar (Ctrl+E).")

    @protegido
    def extrair(self):
        caminho = filedialog.askopenfilename(initialdir=PASTA_SAIDAS, filetypes=[("Touchstone 2 portas", "*.s2p")])
        if not caminho:
            return
        f, s, z0 = ler_touchstone(caminho)
        m = dgs.extrair_pi(f, s, z0)
        s_m = abcd.abcd_para_s(dgs.abcd_pi(f, m), z0)
        ax = self.grafico.eixos(1, 1)
        ax[0].plot(f / 1e9, abcd.db(s[:, 1, 0]), label="|S21| medido/EM")
        ax[0].plot(f / 1e9, abcd.db(s_m[:, 1, 0]), ls="--", label="|S21| modelo π")
        ax[0].plot(f / 1e9, abcd.db(s[:, 0, 0]), label="|S11| medido/EM")
        ax[0].plot(f / 1e9, abcd.db(s_m[:, 0, 0]), ls="--", label="|S11| modelo π")
        ax[0].set_xlabel("f [GHz]")
        ax[0].grid(True, alpha=0.3)
        ax[0].legend(fontsize=7)
        self.grafico.desenhar()
        self.saida.escrever(f"Modelo π extraído de {relativo(caminho)} (eqs. 6.10-6.13):\n"
                            f"  f0 = {m.f0 / 1e9:.4f} GHz\n  R = {m.R:.1f} Ω, L = {m.L * 1e9:.4f} nH, "
                            f"C = {m.C * 1e12:.4f} pF\n  Rp = {m.Rp:.1f} Ω, Cp = {m.Cp * 1e12:.4f} pF")
