"""Aba Linhas acopladas e acopladores direcionais (Pozar, seções 7.6-7.7), com versão LTA."""
import tkinter as tk
from tkinter import ttk

import numpy as np

from ..core import abcd
from ..core import linhas_acopladas as la
from ..core.descontinuidades import texto_limites
from ..core.microstrip import sintese
from ..geometria import circuitos
from .comum import AbaBase, AreaTexto, Formulario, protegido, salvar_csv_s
from .previa import PreviaEstrutura

MODOS = ("Acoplador λ/4 (1 seção)", "Acoplador multisseção binomial", "Lange (4 dedos)",
         "Análise: W, S → Z0e, Z0o")
IMPL = ("Convencional", "LTA: modo ímpar com carga capacitiva periódica",
        "Capacitores entre as linhas nas pontas (Steer 9.10)")


class AbaAcopladas(AbaBase):
    TITULO = "Linhas acopladas e acopladores (Pozar, 7.6-7.7)"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("modo", "Projeto", MODOS, ""),
            ("tipo", "Linha", la.TIPOS, ""),
            ("impl", "Implementação", IMPL, ""),
            ("modelo", "Modelo da microstrip acoplada", la.MODELOS, ""),
            ("z0", "Z0 do sistema", 50.0, "Ω"),
            ("c_db", "Acoplamento", 20.0, "dB"),
            ("n", "Seções (ímpar)", 3, ""),
            ("f0", "Frequência central f0", 3.0, "GHz"),
            ("ncel", "LTA: capacitores por seção", 8, ""),
            ("w", "Análise: W", 1.0, "mm"),
            ("s", "Análise: S", 0.3, "mm"),
            ("b", "Stripline: distância b", 3.2, "mm"),
            ("fmax", "f máx", 6.0, "GHz"),
        ], titulo="Linhas acopladas")
        self.form.pack(fill="x")
        self.perdas = tk.BooleanVar(value=False)
        self.dispersao = tk.BooleanVar(value=False)
        self.chanfro = tk.BooleanVar(value=True)
        ttk.Checkbutton(self.painel, text="Perdas por modo (dielétrico; condutor por Wheeler)",
                        variable=self.perdas).pack(anchor="w", pady=(6, 0))
        ttk.Checkbutton(self.painel, text="Dispersão dos modos (Kirschning-Jansen)",
                        variable=self.dispersao).pack(anchor="w")
        ttk.Checkbutton(self.painel, text="Chanfrar as curvas das alimentações (eq. 3.115)",
                        variable=self.chanfro, command=self.atualizar_previa).pack(anchor="w")
        self.botoes(("Calcular", self.calcular), ("Enviar p/ exportação", self.publicar),
                    ("Ex. 7.7: 20 dB stripline", self.ex_77), ("Ex. 7.8: 3 seções", self.ex_78),
                    ("LTA: 10 dB microstrip", self.ex_lta), ("Caps. nas pontas: 10 dB", self.ex_pontas),
                    ("Exportar S (.csv)", self.exportar))
        self.saida = AreaTexto(self.painel, 22)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Portas (Pozar, Fig. 7.31): 1 entrada, 2 direta, 3 acoplada,\n4 isolada.\n\n"
            "Microstrip acoplada: Kirschning-Jansen (padrão; com espessura,\ndispersão e perdas de Wheeler) "
            "ou Garg-Bahl; stripline:\nCohn (exata, t = 0).\n\n"
            "Na microstrip, εeff do modo par > ímpar: a diretividade cai.\n"
            "A opção LTA carrega só o modo ímpar com capacitores Cc entre\n"
            "as linhas (célula slow-wave do Martín, seção 2.5.2) e iguala\n"
            "as velocidades de fase, recuperando a diretividade.\n"
            "A opção de capacitores nas pontas (Steer 9.10) faz o mesmo\n"
            "com dois capacitores, mas só em torno de f0.")
        self._ultimo = None
        self.observar_previa(self.form)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    def ex_77(self):
        from ..core.microstrip import Substrato
        self.estado.definir_substrato(Substrato(2.2, 1.6, 50.8, 0.0, "Pozar Ex. 7.7 (er = 2,2)"))
        self.form.set(dict(modo=MODOS[0], tipo=la.TIPOS[1], impl=IMPL[0], z0=50, c_db=20, f0=3, b=3.2, fmax=6))
        self.calcular()

    def ex_78(self):
        self.form.set(dict(modo=MODOS[1], impl=IMPL[0], z0=50, c_db=20, n=3, f0=3, fmax=6))
        self.calcular()

    def ex_pontas(self):
        self.form.set(dict(modo=MODOS[0], tipo=la.TIPOS[0], impl=IMPL[2], z0=50, c_db=10, f0=2, fmax=6))
        self.calcular()

    def ex_lta(self):
        self.form.set(dict(modo=MODOS[0], tipo=la.TIPOS[0], impl=IMPL[1], z0=50, c_db=10, f0=2, ncel=8, fmax=6))
        self.calcular()

    # ------------------------------------------------------------------ projeto
    def _projeto(self):
        v = self.form.get()
        sub = self.estado.substrato
        if v["f0"] <= 0 or v["z0"] <= 0:
            raise ValueError("f0 e Z0 devem ser positivos.")
        disp = self.dispersao.get()
        if v["modo"] == MODOS[3]:
            md = la.analisar(v["tipo"], v["w"], v["s"], sub, v["b"], v["modelo"], disp)
            sec = la.Secao(md, la.comprimento_quarto(la._em_f0(md, v["f0"] * 1e9), v["f0"] * 1e9), v["w"], v["s"],
                           tipo=v["tipo"])
            p = la.ProjetoAcoplador("análise", md.z0, v["f0"] * 1e9, [md.c], [sec])
        else:
            lta, pontas = v["impl"] == IMPL[1], v["impl"] == IMPL[2]
            if (lta or pontas) and v["tipo"] != la.TIPOS[0]:
                raise ValueError("A compensação das velocidades é para microstrip (na stripline os modos já "
                                 "têm a mesma velocidade).")
            n = v["n"] if v["modo"] == MODOS[1] else 1
            p = la.projeto_acoplador(v["tipo"], v["z0"], v["c_db"], v["f0"] * 1e9, sub, n, v["b"], lta,
                                     max(1, v["ncel"]), lange=v["modo"] == MODOS[2], modelo=v["modelo"],
                                     dispersivo=disp, pontas=pontas)
        lay = None
        if v["tipo"] == la.TIPOS[0] and v["modo"] != MODOS[2]:
            w0 = float(sintese(v["z0"], sub))
            info = {"tipo": f"acoplador - {p.tipo}", "C_dB": v["c_db"], "f0_GHz": v["f0"],
                    "secoes": len(p.secoes)}
            for i, s in enumerate(p.secoes, 1):
                info[f"secao{i}_W_S_l_mm"] = f"{s.w:.4f} / {s.s:.4f} / {s.comp:.3f}"
                if s.n_cel:
                    info[f"secao{i}_Cc_fF"] = f"{s.cc * 1e15:.2f} x {s.n_cel}"
                if s.c_pontas:
                    info[f"secao{i}_C_pontas_fF"] = f"{s.c_pontas * 1e15:.2f}"
            lay = circuitos.acoplador_linhas("acoplador_" + ("lta" if p.lta else "linhas"), sub,
                                             [(s.comp, s.w, s.s) for s in p.secoes], w0, 5.0, info,
                                             capacitores=[s.n_cel for s in p.secoes] if p.lta else None,
                                             chanfro=self.chanfro.get(),
                                             pontas=any(s.c_pontas for s in p.secoes))
        return v, sub, p, lay

    def atualizar_previa(self):
        try:
            *_, lay = self._projeto()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        if lay is None:
            self.previa.mostrar_erro("Layout só para microstrip (stripline: 3 camadas; Lange: fios de solda).")
        else:
            self.previa.mostrar(lay)

    @protegido
    def calcular(self):
        v, sub, p, lay = self._projeto()
        if lay is not None:
            self.previa.mostrar(lay)
        z0 = v["z0"]
        f = np.linspace(v["fmax"] * 1e9 / 1500, v["fmax"] * 1e9, 1500)
        perdas = self.perdas.get()
        s = la.acoplador(f, p.secoes, z0, sub, perdas)
        # referência: mesmas impedâncias com modos de igual velocidade (TEM ideal)
        ideal = [la.Secao(la.Modos(x.modos.ze, x.modos.zo, 1.0, 1.0),
                          la.comprimento_quarto(la.Modos(1, 1, 1.0, 1.0), p.f0)) for x in p.secoes]
        s_id = la.acoplador(f, ideal, z0)
        self._ultimo = (f, s, p.tipo)
        c, d, t, rl = la.metricas(s)
        c_i, d_i, _, _ = la.metricas(s_id)
        fg = f / 1e9
        ax = self.grafico.eixos(2, 2)
        ax[0].plot(fg, -c, label="|S31| acoplada")
        ax[0].plot(fg, -t, label="|S21| direta")
        ax[0].plot(fg, -c_i, ls="--", lw=0.8, label="|S31| TEM ideal")
        ax[0].set_ylim(-45, 1)
        ax[0].set_ylabel("|S| [dB]")
        ax[0].set_title(f"Acoplador {p.tipo}", fontsize=9)
        ax[1].plot(fg, d, label="diretividade")
        ax[1].plot(fg, np.minimum(d_i, 80), ls="--", lw=0.8, label="TEM ideal")
        ax[1].set_ylim(0, 80)
        ax[1].set_ylabel("D = |S31|/|S41| [dB]")
        ax[1].set_title("Diretividade", fontsize=9)
        ax[2].plot(fg, -rl, label="|S11|")
        ax[2].plot(fg, abcd.db(s[:, 3, 0]), label="|S41| isolada")
        ax[2].set_ylim(-80, 0)
        ax[2].set_ylabel("|S| [dB]")
        ax[2].set_title("Casamento e isolamento", fontsize=9)
        ax[3].plot(fg, np.degrees(np.angle(s[:, 2, 0] * np.conj(s[:, 1, 0]))), label="∠S31 − ∠S21")
        ax[3].set_ylim(-180, 180)
        ax[3].set_ylabel("[°]")
        ax[3].set_title("Quadratura entre as saídas", fontsize=9)
        for a in ax:
            a.axvline(v["f0"], color="gray", ls=":", lw=0.8)
            a.set_xlabel("f [GHz]")
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        i0 = int(np.argmin(np.abs(f - p.f0)))
        linhas = [f"{v['modo']} – {v['tipo']} – {p.tipo}", f"Z0 = {p.z0:.2f} Ω, f0 = {p.f0 / 1e9:g} GHz", ""]
        linhas += p.texto
        if v["modo"] == MODOS[2]:
            md = p.secoes[0].modos
            linhas += [f"  Lange: Ze4 = {md.ze:.2f} Ω, Zo4 = {md.zo:.2f} Ω (eq. 7.97), √(Ze4·Zo4) = "
                       f"{np.sqrt(md.ze * md.zo):.2f} Ω", "  Layout não gerado: dedos interdigitados com fios de "
                       "solda (bond wires); usar as dimensões W e S do par."]
        if v["modo"] == MODOS[3]:
            md = p.secoes[0].modos
            linhas += [f"  W = {v['w']:g} mm, S = {v['s']:g} mm → Z0e = {md.ze:.2f} Ω, Z0o = {md.zo:.2f} Ω",
                       f"  εe,par = {md.eee:.4f}, εe,ímpar = {md.eeo:.4f}, C = {md.c:.4f} "
                       f"({-20 * np.log10(md.c):.2f} dB), Z0 = √(Z0e·Z0o) = {md.z0:.2f} Ω",
                       f"  Par diferencial (Steer 5.13): Z0d = 2·Z0o = {2 * md.zo:.2f} Ω, modo comum "
                       f"Z0c = Z0e/2 = {md.ze / 2:.2f} Ω"]
            if v["tipo"] == la.TIPOS[0]:
                f0 = v["f0"] * 1e9
                if md.disp is not None:
                    e_f = md.disp(f0)
                    linhas.append(f"  Em {v['f0']:g} GHz (dispersão): εe,par = {float(e_f[0]):.4f}, "
                                  f"εe,ímpar = {float(e_f[1]):.4f}")
                ace, aco = la.alfa_condutor_modos(v["w"], v["s"], sub, np.array([f0]), md)
                dle, dlo = la.ponta_aberta_acoplada(v["w"], v["s"], sub)
                linhas += [f"  Perda no condutor em f0 (Wheeler): par {ace[0] * 8.686 / 100:.4f}, ímpar "
                           f"{aco[0] * 8.686 / 100:.4f} dB/cm",
                           f"  Ponta aberta acoplada (eq. 4.155): Δl par = {dle:.4f} mm, ímpar = {dlo:.4f} mm",
                           f"  Modelo: {v['modelo']}" + (f"; espessura t = {sub.t_um:g} µm" if sub.incluir_t else "")]
        linhas += ["", f"Em f0: acoplamento {c[i0]:.2f} dB, direta {t[i0]:.2f} dB, diretividade {d[i0]:.1f} dB, "
                       f"|S11| = {-rl[i0]:.1f} dB"]
        ok = np.abs(c - c[i0]) <= 1.0
        a = b = i0
        while a > 0 and ok[a - 1]:
            a -= 1
        while b < len(f) - 1 and ok[b + 1]:
            b += 1
        linhas.append(f"Acoplamento dentro de ±1 dB: {f[a] / 1e9:.3f}-{f[b] / 1e9:.3f} GHz "
                      f"({100 * (f[b] - f[a]) / p.f0:.0f}%)")
        if any(x.c_pontas for x in p.secoes):
            linhas += ["", "Capacitores nas pontas: um capacitor entre as linhas em cada extremidade de cada seção",
                       "(entre seções eles se somam). Igualam a fase dos modos só em f0: compare a banda de",
                       "diretividade com a opção LTA, que distribui a carga ao longo da seção."]
        if p.lta:
            linhas += ["", "LTA: capacitores Cc entre as linhas, uniformemente espaçados (centro de cada célula);",
                       "realizar com SMD 0201 ou dedos interdigitais e ajustar em EM. O modo par não os vê."]
        if v["tipo"] == la.TIPOS[0]:
            linhas += ["", "Limites do substrato:"] + texto_limites(sub, p.secoes[0].w, v["fmax"] * 1e9)
        linhas += ["", "Junções e dobras das alimentações não modeladas (reduzem a diretividade)."]
        self.saida.escrever("\n".join(linhas))

    @protegido
    def publicar(self):
        *_, lay = self._projeto()
        if lay is None:
            raise ValueError("Layout disponível só para acopladores em microstrip (exceto Lange).")
        self.previa.mostrar(lay)
        self.estado.publicar_layout(lay)
        self.estado.msg(f"Layout '{lay.nome}' enviado para Exportar (Ctrl+E).")

    def exportar(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s, nome = self._ultimo
            salvar_csv_s(self, "acoplador", f, s)
