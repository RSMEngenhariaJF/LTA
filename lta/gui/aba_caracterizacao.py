"""Aba Caracterização por medidas: εr, εeff, α, Z0 e RLGC (Garg, Bahl e Bianchi, seções 1.5 e 2.4.7),
calibração TRL e renormalização de S (Steer, Microwave and RF Design, seções 6.10-6.11)."""
import tkinter as tk
from tkinter import filedialog, ttk

import numpy as np

from ..core import caracterizacao as ca
from ..core import conectores as co
from ..core.microstrip import C0, NP_DB, analise, parametros_linha, sintese
from ..exportacao.simuladores import ler_touchstone, touchstone
from .comum import PASTA_SAIDAS, AbaBase, AreaTexto, Formulario, protegido, relativo

METODOS = ("Duas linhas (.s2p): εeff, α e εr sem os conectores", "Uma linha (.s2p): Zc, γ, RLGC e εr",
           "Dois ressonadores λ/2: εeff e ℓe", "Q0 de um ressonador → α", "VSWR máx./mín. → Z0",
           "Calibração TRL (.s2p): Thru, Line, Reflect → DUT corrigido", "Renormalizar .s2p para outra Zref")
REFLECT = ("Aberto (Γ ≈ +1)", "Curto (Γ ≈ -1)")
IMPEDANCIA = ("Linha de alta impedância (Z0 > Zref)", "Linha de baixa impedância (Z0 < Zref)")


class AbaCaracterizacao(AbaBase):
    TITULO = "Caracterização por medidas"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("metodo", "Método", METODOS, ""),
            ("l1", "Comprimento da linha A / ressonador 1", 30.0, "mm"),
            ("l2", "Comprimento da linha B / ressonador 2", 60.0, "mm"),
            ("w", "Largura W das linhas (para εr)", 1.8, "mm"),
            ("f1", "Ressonância 1", 5.0, "GHz"),
            ("f2", "Ressonância 2", 4.8, "GHz"),
            ("n", "Modo n (nº de meias ondas)", 1, ""),
            ("q0", "Q0 medido", 250.0, ""),
            ("lg", "λg na ressonância", 30.0, "mm"),
            ("vmax", "VSWR máximo", 1.5, ""),
            ("vmin", "VSWR mínimo", 1.0, ""),
            ("zref", "Zref do sistema", 50.0, "Ω"),
            ("imp", "Linha medida", IMPEDANCIA, ""),
            ("dl", "TRL: comprimento extra da Line", 10.0, "mm"),
            ("eest", "TRL: εeff estimado da Line", 2.7, ""),
            ("refl", "TRL: padrão Reflect", REFLECT, ""),
            ("znova", "Renormalizar: nova Zref", 75.0, "Ω"),
        ], titulo="Medida")
        self.form.pack(fill="x")
        self.arq = [tk.StringVar(value="") for _ in range(4)]
        for i, rot in enumerate(("Arquivo A / Thru (.s2p)", "Arquivo B / Line (.s2p)", "TRL: Reflect (S11, S22)",
                                 "TRL: DUT medido")):
            f = ttk.Frame(self.painel)
            f.pack(fill="x", pady=1)
            ttk.Button(f, text=rot, command=lambda k=i: self._escolher(k)).pack(side="left")
            ttk.Label(f, textvariable=self.arq[i], foreground="#555").pack(side="left", padx=4)
        self.botoes(("Calcular", self.calcular), ("Gerar exemplo sintético", self.exemplo))
        self.saida = AreaTexto(self.painel, 20)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Métodos da seção 1.5 (e 2.4.7):\n"
            "• duas linhas iguais de comprimentos diferentes: M1⁻¹·M2 tem\n"
            "  autovalores e^(±γΔL), sem influência de conectores e transições;\n"
            "• uma linha: Zc e γ pelas eqs. 2.147-2.148 (inclui conectores);\n"
            "• dois ressonadores λ/2 (eqs. 1.94-1.95): cancela ponta e gap;\n"
            "• Q0 → α (eq. 1.96); VSWR → Z0 (eq. 1.91).\n"
            "εr sai da inversão de Hammerstad-Jensen com o W e o h do\n"
            "substrato atual. 'Gerar exemplo' cria duas medidas sintéticas\n"
            "(com transições SMA) no substrato atual para testar o fluxo.\n\n"
            "TRL (Steer 6.11, Engen-Hoer): Thru de comprimento zero, Line\n"
            "com ΔL (fase entre 20° e 160° na banda) e Reflect igual nas duas\n"
            "portas (aberto ou curto, só o sinal precisa ser conhecido);\n"
            "remove as caixas de erro e grava o DUT corrigido.\n"
            "Renormalizar: muda a impedância de referência do arquivo A.")

    def _escolher(self, k):
        c = filedialog.askopenfilename(initialdir=PASTA_SAIDAS, filetypes=[("Touchstone 2 portas", "*.s2p")])
        if c:
            self.arq[k].set(c)

    @protegido
    def exemplo(self):
        """Duas linhas de 50 Ω (30 e 60 mm) no substrato atual, com perdas e transições típicas."""
        v = self.form.get()
        if v["metodo"] in METODOS[5:]:
            return self._exemplo_trl(v)
        sub = self.estado.substrato
        w = float(sintese(50, sub))
        f = np.linspace(0.2e9, 10e9, 300)
        z, ee, a = (x[:, 0] for x in parametros_linha(w, sub, f, True, True))
        t = co.PREDEFINIDOS["Típica, VSWR 1,10 em 10 GHz (Tabela 2.5)"]
        PASTA_SAIDAS.mkdir(exist_ok=True)
        nomes = []
        for comp in (v["l1"], v["l2"]):
            from ..core import abcd
            th = 2 * np.pi * f * np.sqrt(ee) / C0 * comp * 1e-3 - 1j * a * comp * 1e-3
            m = abcd.cascata(co.abcd_transicao(f, t, True), abcd.linha(z, th), co.abcd_transicao(f, t, False))
            caminho = PASTA_SAIDAS / f"linha_teste_{comp:g}mm.s2p"
            touchstone(caminho, f, abcd.abcd_para_s(m, 50.0), 50.0,
                       f"linha de teste sintetica {comp:g} mm em {sub.nome} (LTA)")
            nomes.append(caminho)
        self.arq[0].set(str(nomes[0]))
        self.arq[1].set(str(nomes[1]))
        self.form.set(dict(metodo=METODOS[0], w=f"{w:.4f}"))
        self.estado.msg(f"Exemplos salvos em {relativo(nomes[0])} e {relativo(nomes[1])}")
        self.calcular()

    @protegido
    def calcular(self):
        v = self.form.get()
        sub = self.estado.substrato
        met = v["metodo"]
        linhas = [met, ""]
        if met in METODOS[:2]:
            if not self.arq[0].get():
                raise ValueError("Escolha o arquivo A (ou gere o exemplo).")
            f, s1, zr = ler_touchstone(self.arq[0].get())
            if met == METODOS[0]:
                if not self.arq[1].get():
                    raise ValueError("Escolha o arquivo B.")
                f2, s2, _ = ler_touchstone(self.arq[1].get())
                if len(f2) != len(f) or not np.allclose(f2, f):
                    raise ValueError("As duas medidas precisam ter as mesmas frequências.")
                g = ca.duas_linhas(f, s1, s2, v["l1"], v["l2"], zr)
                zc = None
            else:
                ln = ca.linha_de_s(f, s1, v["l1"], zr)
                g, zc = ln.gama, ln.zc
            eeff = (g.imag * C0 / (2 * np.pi * f)) ** 2
            alfa = g.real * NP_DB / 100
            er = np.array([_er_seguro(e, v["w"], sub.h_mm) for e in eeff])
            fg = f / 1e9
            ax = self.grafico.eixos(2, 2)
            ax[0].plot(fg, eeff)
            ax[0].set_ylabel("εeff")
            ax[1].plot(fg, alfa)
            ax[1].set_ylabel("α [dB/cm]")
            ax[2].plot(fg, er)
            ax[2].axhline(sub.er, color="gray", ls=":", lw=0.8, label=f"εr do substrato ({sub.er:g})")
            ax[2].set_ylabel("εr (inversão de Hammerstad-Jensen)")
            ax[2].legend(fontsize=7)
            if zc is not None:
                r, l_, gg, c = ca.Linha(f, zc, g).rlgc()
                ax[3].plot(fg, zc.real, label="Re Zc [Ω]")
                ax[3].plot(fg, zc.imag, label="Im Zc [Ω]")
                ax[3].set_ylabel("Zc [Ω]")
                ax[3].legend(fontsize=7)
                k = len(f) // 2
                linhas += [f"Em {fg[k]:.3f} GHz: Zc = {zc[k].real:.2f} {zc[k].imag:+.2f}j Ω",
                           f"  R = {r[k]:.3f} Ω/m, L = {l_[k] * 1e9:.3f} nH/m, G = {gg[k] * 1e3:.4f} mS/m, "
                           f"C = {c[k] * 1e12:.3f} pF/m",
                           "  (inclui transições e conectores; para excluí-los use o método de duas linhas)"]
            else:
                ax[3].plot(fg, np.full_like(fg, float(analise(v["w"], sub)[1])), ls="--",
                           label="εeff do modelo (estático)")
                ax[3].plot(fg, eeff, label="εeff medido")
                ax[3].set_ylabel("εeff")
                ax[3].legend(fontsize=7)
            for a in ax:
                a.set_xlabel("f [GHz]")
                a.grid(True, alpha=0.3)
            self.grafico.desenhar()
            sel = slice(len(f) // 10, None)
            linhas += [f"εeff médio = {np.mean(eeff[sel]):.4f} (de {eeff[sel].min():.4f} a {eeff[sel].max():.4f})",
                       f"εr médio = {np.nanmean(er[sel]):.4f} (substrato atual: {sub.er:g}, h = {sub.h_mm:g} mm, "
                       f"W = {v['w']:g} mm)",
                       f"α em {fg[-1]:g} GHz = {alfa[-1]:.4f} dB/cm",
                       "A inversão para εr é quase-estática: com dispersão, o εr aparente sobe com f; use as",
                       "frequências baixas (fh pequeno) para estimar o εr do substrato."]
        elif met == METODOS[5]:
            linhas += self._trl(v)
        elif met == METODOS[6]:
            linhas += self._renormalizar(v)
        elif met == METODOS[2]:
            ee, le = ca.dois_ressonadores(v["f1"] * 1e9, v["f2"] * 1e9, v["l1"], v["l2"], v["n"])
            linhas += [f"εeff = {ee:.4f} (eq. 1.94); ℓe = Δℓg + Δℓoc = {le:.4f} mm (eq. 1.95)",
                       f"εr = {_er_seguro(ee, v['w'], sub.h_mm):.4f} para W = {v['w']:g} mm, h = {sub.h_mm:g} mm"]
        elif met == METODOS[3]:
            a = ca.alfa_de_q0(v["q0"], v["lg"])
            linhas += [f"α = 27,3/(Q0·λg) = {a:.4f} dB/cm (eq. 1.96)",
                       "Inclui todas as perdas do ressonador (condutor, dielétrico e radiação)."]
        else:
            z = ca.z0_de_vswr(v["zref"], v["vmax"], v["vmin"], v["imp"] == IMPEDANCIA[0])
            linhas += [f"Z0 = {z:.2f} Ω (eq. 1.91; o sinal do expoente vem de saber se a linha é de alta ou baixa "
                       "impedância, por exemplo por TDR)"]
        self.saida.escrever("\n".join(linhas))


    def _exemplo_trl(self, v):
        """Medidas sintéticas de TRL: caixas de erro = transição SMA + 10 mm de linha de 50 Ω; DUT =
        trecho de 8 mm de 30 Ω; Line com ΔL = λ/4 no centro da banda; Reflect aberto ideal."""
        from ..core import abcd
        sub = self.estado.substrato
        w = float(sintese(50, sub))
        f = np.linspace(1e9, 8e9, 281)
        z, ee, a = (x[:, 0] for x in parametros_linha(w, sub, f, True, True))
        t = co.PREDEFINIDOS["Típica, VSWR 1,10 em 10 GHz (Tabela 2.5)"]

        def lin(zz, eeff, alfa, comp):
            return abcd.linha(zz, 2 * np.pi * f * np.sqrt(eeff) / C0 * comp * 1e-3 - 1j * alfa * comp * 1e-3)
        caixa_a = abcd.cascata(co.abcd_transicao(f, t, True), lin(z, ee, a, 10.0))
        caixa_b = abcd.cascata(lin(z, ee, a, 10.0), co.abcd_transicao(f, t, False))
        dl = C0 / (4 * 4.5e9 * np.sqrt(np.mean(ee))) * 1e3
        w30 = float(sintese(30, sub))
        z3, e3, a3 = (x[:, 0] for x in parametros_linha(w30, sub, f, True, True))
        dut = lin(z3, e3, a3, 8.0)
        sa, sb = abcd.abcd_para_s(caixa_a, 50.0), abcd.abcd_para_s(caixa_b, 50.0)
        gm1 = sa[:, 0, 0] + sa[:, 0, 1] * sa[:, 1, 0] / (1 - sa[:, 1, 1])
        gm2 = sb[:, 1, 1] + sb[:, 1, 0] * sb[:, 0, 1] / (1 - sb[:, 0, 0])
        refl = np.zeros((len(f), 2, 2), dtype=complex)
        refl[:, 0, 0], refl[:, 1, 1] = gm1, gm2
        PASTA_SAIDAS.mkdir(exist_ok=True)
        medidas = {"thru": abcd.abcd_para_s(abcd.cascata(caixa_a, caixa_b), 50.0),
                   "line": abcd.abcd_para_s(abcd.cascata(caixa_a, lin(z, ee, a, dl), caixa_b), 50.0),
                   "reflect": refl, "dut": abcd.abcd_para_s(abcd.cascata(caixa_a, dut, caixa_b), 50.0)}
        nomes = []
        for k, sm in medidas.items():
            caminho = PASTA_SAIDAS / f"trl_{k}.s2p"
            touchstone(caminho, f, sm, 50.0, f"TRL sintetico ({k}) em {sub.nome}; DUT: 8 mm de 30 ohm")
            nomes.append(caminho)
        # o TRL refere o DUT à Z0 da Line (aqui Z(f) da microstrip dispersiva, ~50-53 Ω), não a 50 Ω
        touchstone(PASTA_SAIDAS / "trl_dut_verdadeiro.s2p", f, abcd.abcd_para_s(dut, z), 50.0,
                   "DUT sem as caixas de erro, referido a Z0(f) da Line (referencia do exemplo TRL)")
        for var, n in zip(self.arq, nomes):
            var.set(str(n))
        self.form.set(dict(metodo=METODOS[5], dl=f"{dl:.4f}", eest=f"{np.mean(ee):.4f}", refl=REFLECT[0]))
        self.estado.msg(f"Medidas TRL sintéticas salvas em {relativo(PASTA_SAIDAS)} (trl_*.s2p)")
        self.calcular()

    def _ler(self, k, nome):
        if not self.arq[k].get():
            raise ValueError(f"Escolha o arquivo {nome} (ou gere o exemplo).")
        return ler_touchstone(self.arq[k].get())

    def _trl(self, v):
        f, st, zr = self._ler(0, "Thru")
        dados = [self._ler(k, n) for k, n in ((1, "Line"), (2, "Reflect"), (3, "DUT"))]
        for fk, _, _ in dados:
            if len(fk) != len(f) or not np.allclose(fk, f):
                raise ValueError("Thru, Line, Reflect e DUT precisam ter as mesmas frequências.")
        sl, sr, sd = (d[1] for d in dados)
        r = ca.trl(f, st, sl, sr[:, 0, 0], sr[:, 1, 1], sd, v["dl"], v["eest"], v["refl"] == REFLECT[0])
        fase = np.degrees(r.gama.imag * v["dl"] * 1e-3)
        eeff = (r.gama.imag * C0 / (2 * np.pi * f)) ** 2
        fg = f / 1e9
        ax = self.grafico.eixos(2, 2)
        for (i, j), rot in (((0, 0), "S11"), ((1, 0), "S21")):
            ax[0].plot(fg, 20 * np.log10(np.abs(r.s_dut[:, i, j]) + 1e-15), label=f"|{rot}| corrigido")
            ax[0].plot(fg, 20 * np.log10(np.abs(sd[:, i, j]) + 1e-15), ls="--", lw=0.8, label=f"|{rot}| medido")
        ax[0].set_ylabel("[dB]")
        ax[0].set_title("DUT: medido × corrigido (TRL)", fontsize=9)
        ax[1].plot(fg, np.degrees(np.angle(r.s_dut[:, 1, 0])), label="∠S21 corrigido")
        ax[1].plot(fg, np.degrees(np.angle(sd[:, 1, 0])), ls="--", lw=0.8, label="∠S21 medido")
        ax[1].set_ylabel("[°]")
        ax[2].plot(fg, fase, label="fase da Line (βΔL)")
        ax[2].axhspan(20, 160, color="green", alpha=0.08, label="faixa útil 20°-160°")
        ax[2].set_ylabel("[°]")
        ax[3].plot(fg, eeff, label="εeff da Line")
        ax[3].plot(fg, r.gama.real * NP_DB / 100, label="α [dB/cm]")
        for a in ax:
            a.set_xlabel("f [GHz]")
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        PASTA_SAIDAS.mkdir(exist_ok=True)
        saida = PASTA_SAIDAS / "dut_corrigido_trl.s2p"
        touchstone(saida, f, r.s_dut, zr, "DUT corrigido por TRL (LTA)")
        ruim = (fase < 20) | (fase > 160)
        k = len(f) // 2
        linhas = [f"DUT corrigido gravado em {relativo(saida)} (referido à Z0 da Line; gravado como {zr:g} Ω nominal)",
                  f"Em {fg[k]:.3f} GHz: |S11| = {20 * np.log10(abs(r.s_dut[k, 0, 0]) + 1e-15):.2f} dB, "
                  f"|S21| = {20 * np.log10(abs(r.s_dut[k, 1, 0]) + 1e-15):.2f} dB ∠{np.degrees(np.angle(r.s_dut[k, 1, 0])):.2f}°",
                  f"Reflect estimado: Γ = {abs(r.gama_r[k]):.4f} ∠{np.degrees(np.angle(r.gama_r[k])):.2f}° "
                  f"(sinal escolhido: {v['refl']})",
                  f"εeff da Line: {np.mean(eeff):.4f} (média); α em {fg[-1]:g} GHz = "
                  f"{r.gama.real[-1] * NP_DB / 100:.4f} dB/cm"]
        if ruim.any():
            linhas.append(f"Atenção: fase da Line fora de 20°-160° em {100 * ruim.mean():.0f}% dos pontos "
                          "(TRL mal condicionado; use outra Line ou divida a banda).")
        verdadeiro = PASTA_SAIDAS / "trl_dut_verdadeiro.s2p"
        if verdadeiro.exists() and self.arq[3].get().endswith("trl_dut.s2p"):
            fv, sv, _ = ler_touchstone(verdadeiro)
            if len(fv) == len(f):
                linhas.append(f"Exemplo: erro máximo |S corrigido - S verdadeiro| = {np.max(np.abs(r.s_dut - sv)):.2e}")
        linhas += ["", "O plano de referência fica no meio do Thru; a Zref é a impedância característica da Line",
                   "(renormalize para 50 Ω se a Line não for de 50 Ω)."]
        return linhas

    def _renormalizar(self, v):
        f, s, zr = self._ler(0, "A")
        z1 = v["znova"]
        if z1 <= 0:
            raise ValueError("A nova Zref deve ser positiva.")
        s1 = ca.renormalizar(s, zr, z1)
        fg = f / 1e9
        ax = self.grafico.eixos(1, 2)
        for (i, j), rot in (((0, 0), "S11"), ((1, 0), "S21"), ((1, 1), "S22")):
            l, = ax[0].plot(fg, 20 * np.log10(np.abs(s1[:, i, j]) + 1e-15), label=f"|{rot}| em {z1:g} Ω")
            ax[0].plot(fg, 20 * np.log10(np.abs(s[:, i, j]) + 1e-15), ls="--", lw=0.8, color=l.get_color(),
                       label=f"|{rot}| em {zr:g} Ω")
        ax[0].set_ylabel("[dB]")
        ax[1].plot(fg, np.degrees(np.angle(s1[:, 1, 0])), label=f"∠S21 em {z1:g} Ω")
        ax[1].plot(fg, np.degrees(np.angle(s[:, 1, 0])), ls="--", lw=0.8, label=f"∠S21 em {zr:g} Ω")
        ax[1].set_ylabel("[°]")
        for a in ax:
            a.set_xlabel("f [GHz]")
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        from pathlib import Path
        orig = Path(self.arq[0].get())
        saida = PASTA_SAIDAS / f"{orig.stem}_z{z1:g}.s2p"
        PASTA_SAIDAS.mkdir(exist_ok=True)
        touchstone(saida, f, s1, z1, f"{orig.name} renormalizado de {zr:g} para {z1:g} ohm (LTA)")
        return [f"{orig.name}: Zref {zr:g} Ω → {z1:g} Ω (Steer 6.10: via matriz Z)",
                f"Gravado em {relativo(saida)}",
                f"Em {fg[len(f) // 2]:.3f} GHz: |S11| {20 * np.log10(abs(s[len(f) // 2, 0, 0]) + 1e-15):.2f} → "
                f"{20 * np.log10(abs(s1[len(f) // 2, 0, 0]) + 1e-15):.2f} dB"]


def _er_seguro(eeff, w, h):
    try:
        return ca.er_de_eeff(float(eeff), w, h)
    except ValueError:
        return np.nan
