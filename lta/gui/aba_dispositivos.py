"""Aba Híbridas e divisores: versões convencionais (Pozar, cap. 7) e com linhas artificiais
(Martín, seções 4.2.1-4.2.2 e 2.5.2)."""
import tkinter as tk
from tkinter import ttk

import numpy as np

from ..core import abcd
from ..core import dispositivos as dv
from ..core.descontinuidades import texto_limites
from ..geometria import circuitos
from .comum import AbaBase, AreaTexto, Formulario, OpcoesModelo, protegido, salvar_csv_s
from .previa import PreviaEstrutura

DISPOSITIVOS = ("Branch-line (híbrida 90°)", "Rat-race (híbrida 180°)", "Wilkinson", "Divisor com inversor (Y)",
                "Wilkinson desigual", "Junção T sem perdas")
IMPL = ("Dual-band CRLH (-90° em f1, +90° em f2)", "Banda larga: rat-race com -90° LH no lugar de 270°",
        "Convencional em microstrip (Pozar)", "Slow-wave: linha com carga capacitiva (miniaturizada)",
        "Concentrada: π LC no lugar das linhas (Steer 7.13)")
SIM_NAO = ("Sim", "Não")


class AbaDispositivos(AbaBase):
    TITULO = "Híbridas e divisores"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("disp", "Dispositivo", DISPOSITIVOS, ""),
            ("impl", "Implementação", IMPL, ""),
            ("z0", "Z0 das portas", 50.0, "Ω"),
            ("f0", "f0 (convencional e slow-wave)", 1.6, "GHz"),
            ("f1", "f1 (banda LH)", 0.9, "GHz"),
            ("f2", "f2 (banda RH)", 1.8, "GHz"),
            ("nlh", "Células da linha LH (rat-race)", 1, ""),
            ("k2", "Divisão P3/P2 (desigual, T)", 2.0, ""),
            ("trafo", "Transformadores λ/4 de saída", SIM_NAO, ""),
            ("ncel", "Slow-wave: células por λ/4", 3, ""),
            ("zo", "Slow-wave: Z da linha não carregada", 100.0, "Ω"),
            ("carga", "Slow-wave: carga", dv.CARGAS_SW, ""),
            ("zc", "Slow-wave: Z do stub/patch", 60.0, "Ω"),
            ("gr", "Gap do resistor (Wilkinson)", 1.0, "mm"),
            ("fmax", "f máx", 3.0, "GHz"),
        ], titulo="Dispositivo")
        self.form.pack(fill="x")
        self.modelo = OpcoesModelo(self.painel)
        self.modelo.pack(fill="x", pady=(6, 0))
        self.chanfro = tk.BooleanVar(value=True)
        ttk.Checkbutton(self.painel, text="Chanfrar curvas no layout (Garg-Bahl eqs. 3.115-3.116)",
                        variable=self.chanfro, command=self.atualizar_previa).pack(anchor="w")
        self.botoes(("Calcular", self.calcular), ("Enviar p/ exportação", self.publicar),
                    ("Ex. branch-line GSM (Fig. 4.16)", self.ex_branch),
                    ("Ex. divisor 2,4/3,75 GHz (Fig. 4.11)", self.ex_divisor),
                    ("Ex. rat-race banda larga (Fig. 4.5)", self.ex_ratrace),
                    ("Ex. branch-line slow-wave", self.ex_slow), ("Exportar S (.csv)", self.exportar))
        self.saida = AreaTexto(self.painel, 22)
        self.saida.pack(fill="both", expand=True)
        self._ultimo = None
        self.observar_previa(self.form)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    def ex_branch(self):
        """Híbrida dual-band em 0,9 e 1,8 GHz (livro, Fig. 4.16, aqui com CRLH canônica ideal)."""
        self.form.set(dict(disp=DISPOSITIVOS[0], impl=IMPL[0], z0=50, f0=0.9, f1=0.9, f2=1.8, fmax=3))
        self.calcular()

    def ex_divisor(self):
        """Divisor dual-band com inversor de 35,35 Ω em 2,4 e 3,75 GHz (livro, Fig. 4.11)."""
        self.form.set(dict(disp=DISPOSITIVOS[3], impl=IMPL[0], z0=50, f0=2.4, f1=2.4, f2=3.75, fmax=6))
        self.calcular()

    def ex_ratrace(self):
        """Rat-race em 1,6 GHz: 270° convencional contra -90° LH (livro, Figs. 4.3 e 4.5)."""
        self.form.set(dict(disp=DISPOSITIVOS[1], impl=IMPL[1], z0=50, f0=1.6, nlh=3, fmax=3.2))
        self.calcular()

    def ex_slow(self):
        """Branch-line miniaturizada com stubs (fator slow-wave ≈ Zo/Z)."""
        self.form.set(dict(disp=DISPOSITIVOS[0], impl=IMPL[3], z0=50, f0=2.0, ncel=3, zo=100,
                           carga=dv.CARGAS_SW[0], zc=60, fmax=4.0))
        self.calcular()

    # ------------------------------------------------------------------ montagem
    def _impedancias(self, v):
        """Impedâncias dos ramos e transformadores de saída (None se não houver)."""
        z0, disp = v["z0"], v["disp"]
        trafo = v["trafo"] == SIM_NAO[0]
        if disp == DISPOSITIVOS[0]:
            return dict(ramos=(z0 / np.sqrt(2), z0))
        if disp == DISPOSITIVOS[1]:
            return dict(ramos=(z0 * np.sqrt(2),))
        if disp == DISPOSITIVOS[2]:
            return dict(ramos=(z0 * np.sqrt(2),))
        if disp == DISPOSITIVOS[3]:
            return dict(ramos=(z0 / np.sqrt(2),))
        if v["k2"] <= 0:
            raise ValueError("A razão P3/P2 deve ser positiva.")
        if disp == DISPOSITIVOS[4]:
            z2, z3, r, r2, r3 = dv.impedancias_wilkinson(z0, v["k2"])
            return dict(ramos=(z2, z3), r=r, saidas=(r2, r3),
                        trafos=(np.sqrt(z0 * r2), np.sqrt(z0 * r3)) if trafo else None)
        z2, z3 = dv.impedancias_juncao_t(z0, v["k2"])
        return dict(ramos=(), saidas=(z2, z3), trafos=(np.sqrt(z0 * z2), np.sqrt(z0 * z3)) if trafo else None)

    def _linha(self, v, z, fase, convencional_ideal=False):
        """(função f -> ABCD, descrição, dados físicos) da linha de impedância z e fase em f0."""
        f0 = v["f0"] * 1e9
        sub = self.estado.substrato
        perdas, disp_ = self.modelo.get()
        impl = v["impl"]
        if convencional_ideal:
            return dv.linha_convencional(z, fase, f0), None, None
        if impl == IMPL[2]:
            m, (w, comp, ee) = dv.linha_microstrip(z, fase, f0, sub, perdas, disp_)
            return m, f"Z = {z:.2f} Ω: W = {w:.4f} mm, ℓ = {comp:.3f} mm ({np.degrees(fase):.0f}°)", (w, comp)
        if impl == IMPL[3]:
            n = max(1, int(round(v["ncel"] * abs(fase) / (np.pi / 2))))
            p = dv.projeto_slow_wave(z, fase, f0, n, v["zo"], sub, v["zc"], v["carga"])
            g = p.geo
            txt = (f"Z = {z:.2f} Ω ({np.degrees(fase):.0f}°): {n} células de {g.periodo:.3f} mm, linha Zo = "
                   f"{v['zo']:g} Ω (W = {g.w_linha:.3f}), {g.tipo.lower()[:-1]} de {g.z_carga:g} Ω "
                   f"{g.w_carga:.3f} × {g.l_carga:.3f} mm, Cls = {p.carga.Cls * 1e12:.4f} pF; ℓ = {p.comp:.2f} mm "
                   f"({100 * p.reducao:.0f}% de {p.comp_convencional:.2f} mm)")
            return dv.linha_slow_wave(p, sub, perdas, disp_), txt, p
        if impl == IMPL[4]:
            m, (se, de) = dv.linha_concentrada(z, fase, f0)
            fmt = lambda c: f"L = {c[1] * 1e9:.4g} nH" if c[0] == "L" else f"C = {c[1] * 1e12:.4g} pF"  # noqa: E731
            return m, (f"Z = {z:.2f} Ω ({np.degrees(fase):.0f}°): π com série {fmt(se)} e derivações "
                       f"{fmt(de)} cada"), None
        raise ValueError("Implementação sem linha de fase única.")

    def _montar(self, v, f, artificial):
        """S do dispositivo; artificial=False dá a versão TEM ideal em f0 (referência)."""
        z0, f0 = v["z0"], v["f0"] * 1e9
        disp, impl = v["disp"], v["impl"]
        zz = self._impedancias(v)
        ramos = zz["ramos"]
        info, fis = [], {}
        ideal = not artificial
        l270 = None
        if ideal or impl in IMPL[2:]:
            mk = lambda z, fase: self._linha(v, z, fase, ideal)      # noqa: E731
            ls = []
            for i, z in enumerate(ramos):
                m, txt, d = mk(z, np.pi / 2)
                ls.append(m)
                fis[f"ramo{i}"] = d
                if txt and (i == 0 or z != ramos[0]):
                    info.append("  " + txt)
            if disp == DISPOSITIVOS[1]:
                l270, txt, d = mk(ramos[0], 3 * np.pi / 2)
                fis["ramo270"] = d
                if txt:
                    info.append("  " + txt)
            tr = None
            if zz.get("trafos"):
                tr = []
                for i, z in enumerate(zz["trafos"]):
                    m, txt, d = mk(z, np.pi / 2)
                    tr.append(m)
                    fis[f"trafo{i}"] = d
                    if txt:
                        info.append("  trafo: " + txt)
        elif impl == IMPL[0]:
            if disp == DISPOSITIVOS[1]:
                raise ValueError("Rat-race dual-band exige ramo de 270° nas duas bandas; use a opção "
                                 "de banda larga para o rat-race.")
            f1, f2 = v["f1"] * 1e9, v["f2"] * 1e9
            ls = []
            todos = list(ramos) + list(zz.get("trafos") or [])
            for i, z in enumerate(todos):
                func, c = dv.linha_crlh_dual(z, f1, f2)
                ls.append(func)
                info.append(f"  Z = {z:.2f} Ω: LR = {c.LR * 1e9:.4g} nH, CR = {c.CR * 1e12:.4g} pF, "
                            f"LL = {c.LL * 1e9:.4g} nH, CL = {c.CL * 1e12:.4g} pF")
            tr = ls[len(ramos):] or None
            ls = ls[:len(ramos)]
        else:
            if disp != DISPOSITIVOS[1]:
                raise ValueError("A opção de banda larga (linha LH de -90°) se aplica ao rat-race.")
            ls = [dv.linha_convencional(z, np.pi / 2, f0) for z in ramos]
            l270, (cl, ll) = dv.linha_plh(ramos[0], -np.pi / 2, f0, max(1, v["nlh"]))
            tr = None
            info.append(f"  Linha LH de -90° em f0: {max(1, v['nlh'])} célula(s) T, CL = {cl * 1e12:.4g} pF, "
                        f"LL = {ll * 1e9:.4g} nH (Z = {ramos[0]:.2f} Ω)")
        if disp == DISPOSITIVOS[0]:
            s = dv.branch_line(f, ls[0], ls[1], z0)
        elif disp == DISPOSITIVOS[1]:
            s = dv.rat_race(f, ls[0], l270, z0)
        elif disp == DISPOSITIVOS[2]:
            s = dv.wilkinson(f, ls[0], z0)
        elif disp == DISPOSITIVOS[3]:
            s = dv.divisor_inversor(f, ls[0], z0)
        elif disp == DISPOSITIVOS[4]:
            s = dv.wilkinson_desigual(f, ls[0], ls[1], v["k2"], z0, *(tr or (None, None)))
        else:
            s = dv.juncao_t(f, z0, *zz["saidas"], *(tr or (None, None)))
        return s, info, fis, zz

    def _layout(self, v):
        """Layout microstrip (convencional ou slow-wave) ou None."""
        if v["impl"] not in IMPL[2:4]:
            return None
        sub = self.estado.substrato
        f = np.array([v["f0"] * 1e9])
        _, _, fis, zz = self._montar(v, f, True)
        w0 = float(dv.sintese(v["z0"], sub))
        disp = v["disp"]
        info = {"tipo": f"{disp} - {v['impl']}", "f0_GHz": v["f0"], "Z0_ohm": v["z0"]}
        nome = {DISPOSITIVOS[0]: "branch_line", DISPOSITIVOS[1]: "rat_race", DISPOSITIVOS[2]: "wilkinson",
                DISPOSITIVOS[3]: "divisor_y", DISPOSITIVOS[4]: "wilkinson_desigual",
                DISPOSITIVOS[5]: "juncao_t"}[disp] + ("_slow_wave" if v["impl"] == IMPL[3] else "")
        slow = v["impl"] == IMPL[3]

        def wl(chave):
            d = fis[chave]
            return (d.geo.w_linha, d.comp) if slow else d

        if disp == DISPOSITIVOS[0]:
            if slow:
                return circuitos.branch_line(nome, sub, w0, 0, 0, 0, 0, info=info,
                                             sw_s=(fis["ramo0"].geo, fis["ramo0"].n),
                                             sw_p=(fis["ramo1"].geo, fis["ramo1"].n))
            (ws, ls_), (wp, lp_) = fis["ramo0"], fis["ramo1"]
            return circuitos.branch_line(nome, sub, w0, ws, ls_, wp, lp_, info=info)
        if slow:
            return None                        # outros dispositivos slow-wave: só circuito
        if disp == DISPOSITIVOS[1]:
            w_r, l90 = fis["ramo0"]
            return circuitos.rat_race(nome, sub, w0, w_r, 6 * l90 / (2 * np.pi), info=info,
                                      chanfro=self.chanfro.get())
        if disp in (DISPOSITIVOS[2], DISPOSITIVOS[4]):
            w2, l2 = wl("ramo0")
            w3, l3 = wl("ramo1") if disp == DISPOSITIVOS[4] else (w2, l2)
            wo2 = wo3 = w0
            t2 = t3 = None
            if disp == DISPOSITIVOS[4]:
                r2, r3 = zz["saidas"]
                wo2, wo3 = float(dv.sintese(r2, sub)), float(dv.sintese(r3, sub))
                if zz.get("trafos"):
                    t2, t3 = fis["trafo0"], fis["trafo1"]
                info["R_ohm"] = f"{zz['r']:.2f}"
            else:
                info["R_ohm"] = f"{2 * v['z0']:.2f}"
            return circuitos.wilkinson(nome, sub, w0, w2, l2, w3, l3, wo2, wo3, v["gr"], t2, t3, info=info,
                                       chanfro_f0=v["f0"] * 1e9 if self.chanfro.get() else None)
        if disp == DISPOSITIVOS[5]:
            z2, z3 = zz["saidas"]
            t2 = t3 = None
            if zz.get("trafos"):
                t2, t3 = fis["trafo0"], fis["trafo1"]
            return circuitos.juncao_t(nome, sub, w0, float(dv.sintese(z2, sub)), float(dv.sintese(z3, sub)),
                                      t2, t3, info=info, chanfro=self.chanfro.get())
        return None

    def atualizar_previa(self):
        try:
            v = self.form.get()
            lay = self._layout(v)
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        if lay is None:
            self.previa.mostrar_erro("Layout gerado para as versões em microstrip (convencional; slow-wave só "
                                     "na branch-line). CRLH e LH: elementos concentrados.")
        else:
            self.previa.mostrar(lay)

    def _metricas(self, v, f, s, nome, f_proj):
        """Faixas em torno de cada frequência de projeto (critérios da Tabela 4.1)."""
        disp = v["disp"]
        s11 = abcd.db(s[:, 0, 0])
        a, b = abcd.db(s[:, 1, 0]), abcd.db(s[:, 2, 0])
        linhas = [f"[{nome}]"]
        if disp in DISPOSITIVOS[4:]:
            k2 = v["k2"]
            a0, b0 = 10 * np.log10(1 / (1 + k2)), 10 * np.log10(k2 / (1 + k2))
        else:
            a0 = b0 = -3.01
        crit = {"|S11| < -10 dB": s11 < -10,
                f"saídas {a0:.2f}/{b0:.2f} ± 1 dB": (np.abs(a - a0) < 1) & (np.abs(b - b0) < 1)}
        if disp == DISPOSITIVOS[0]:
            crit["fase S21-S31 = ±90° ± 5°"] = np.abs(np.abs(dv.diferenca_fase(s[:, 1, 0], s[:, 2, 0])) - 90) < 5
            crit["isolamento |S41| < -15 dB"] = abcd.db(s[:, 3, 0]) < -15
        elif disp == DISPOSITIVOS[1]:
            crit["Δ: fase S21-S31 = 180° ± 5°"] = np.abs(np.abs(dv.diferenca_fase(s[:, 1, 0], s[:, 2, 0])) - 180) < 5
            crit["Σ: fase S24-S34 = 0° ± 5°"] = np.abs(dv.diferenca_fase(s[:, 1, 3], s[:, 2, 3])) < 5
            crit["isolamento |S41| < -15 dB"] = abcd.db(s[:, 3, 0]) < -15
        elif disp in (DISPOSITIVOS[2], DISPOSITIVOS[4]):
            crit["isolamento |S32| < -15 dB"] = abcd.db(s[:, 2, 1]) < -15
        todos = np.all(list(crit.values()), axis=0)
        for fp in f_proj:
            linhas.append(f"  em {fp / 1e9:.3f} GHz:")
            for nome_c, ok in list(crit.items()) + [("TODOS os critérios", todos)]:
                fx = dv.faixa_em_torno(f, ok, fp)
                linhas.append(f"    {nome_c:32s} " + (f"{fx[0] / 1e9:.3f}-{fx[1] / 1e9:.3f} GHz "
                                                       f"({100 * (fx[1] - fx[0]) / fp:.1f}%)" if fx else "não atende"))
        return linhas

    @protegido
    def calcular(self):
        v = self.form.get()
        f = np.linspace(v["fmax"] * 1e9 / 2000, v["fmax"] * 1e9, 2000)
        s_art, info, _, zz = self._montar(v, f, True)
        s_conv, _, _, _ = self._montar(v, f, False)
        self._ultimo = (f, s_art, v)
        try:
            lay = self._layout(v)
        except ValueError:
            lay = None
        if lay is not None:
            self.previa.mostrar(lay)
        disp = v["disp"]
        fg = f / 1e9
        ax = self.grafico.eixos(2, 2)
        n_portas = s_art.shape[1]
        for j in range(n_portas):
            ax[0].plot(fg, abcd.db(s_art[:, j, 0]), label=f"|S{j + 1}1|")
            ax[1].plot(fg, abcd.db(s_conv[:, j, 0]), ls="--", label=f"|S{j + 1}1|")
        ax[0].set_title(v["impl"].split(":")[0], fontsize=9)
        ax[1].set_title(f"TEM ideal em f0 = {v['f0']:g} GHz", fontsize=9)
        for a in ax[:2]:
            a.set_ylim(-40, 1)
            a.set_ylabel("|S| [dB]")
        if disp == DISPOSITIVOS[1]:
            pares = (("Δ: ∠S21 - ∠S31", (1, 0), (2, 0)), ("Σ: ∠S24 - ∠S34", (1, 3), (2, 3)))
        else:
            pares = (("∠S21 - ∠S31", (1, 0), (2, 0)),)
        for nome, (i1, j1), (i2, j2) in pares:
            ax[2].plot(fg, dv.diferenca_fase(s_art[:, i1, j1], s_art[:, i2, j2]), label=f"{nome}")
            ax[2].plot(fg, dv.diferenca_fase(s_conv[:, i1, j1], s_conv[:, i2, j2]), ls="--",
                       label=f"{nome} TEM ideal")
        ax[2].set_ylim(-190, 190)
        ax[2].set_ylabel("Balanço de fase [°]")
        ax[2].set_title("Balanço de fase entre as saídas", fontsize=9)
        ax[3].plot(fg, np.abs(s_art[:, 1, 0]) ** 2 + np.abs(s_art[:, 2, 0]) ** 2, label="|S21|²+|S31|²")
        ax[3].plot(fg, np.abs(s_conv[:, 1, 0]) ** 2 + np.abs(s_conv[:, 2, 0]) ** 2, ls="--",
                   label="|S21|²+|S31|² TEM ideal")
        ax[3].set_ylim(0, 1.05)
        ax[3].set_title("Potência entregue às saídas", fontsize=9)
        for a in ax:
            a.set_xlabel("f [GHz]")
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        f_proj = [v["f1"] * 1e9, v["f2"] * 1e9] if v["impl"] == IMPL[0] else [v["f0"] * 1e9]
        linhas = [disp, f"Implementação: {v['impl']}"]
        if disp == DISPOSITIVOS[4]:
            linhas.append(f"  Pozar eq. 7.37: Z02 = {zz['ramos'][0]:.2f} Ω, Z03 = {zz['ramos'][1]:.2f} Ω, "
                          f"R = {zz['r']:.2f} Ω, saídas R2 = {zz['saidas'][0]:.2f} e R3 = {zz['saidas'][1]:.2f} Ω")
        if disp == DISPOSITIVOS[5]:
            linhas.append(f"  Saídas Z2 = {zz['saidas'][0]:.2f} Ω, Z3 = {zz['saidas'][1]:.2f} Ω (Ex. 7.1); "
                          + ("transformadores λ/4 para Z0" if zz.get("trafos") else
                             f"portas 2 e 3 com referência Z2 e Z3: |S22| = {abs(s_art[0, 1, 1]):.3f}, "
                             f"|S33| = {abs(s_art[0, 2, 2]):.3f}, sem isolamento"))
        linhas += info + [""]
        linhas += self._metricas(v, f, s_art, v["impl"].split(":")[0], f_proj) + [""]
        linhas += self._metricas(v, f, s_conv, f"TEM ideal, f0 = {v['f0']:g} GHz", [v["f0"] * 1e9])
        if lay is not None:
            linhas += ["", f"Layout: {lay.comprimento:.2f} × {lay.largura:.2f} mm (placa)."]
            if disp == DISPOSITIVOS[1]:
                linhas.append("Rat-race: alimentações de comprimentos diferentes; compare fases após de-embedding.")
            if disp in (DISPOSITIVOS[2], DISPOSITIVOS[4]):
                linhas.append(f"Resistor SMD entre as pontas dos braços (gap de {v['gr']:g} mm).")
        if v["impl"] in IMPL[2:4]:
            sub = self.estado.substrato
            linhas += ["", self.modelo.descrever(sub), "Limites do substrato:"] + texto_limites(sub, None, f[-1])
        linhas += ["", "Junções T, curvas e cantos não compensados (Pozar: alongar os ramos",
                   "10-20° na branch-line). Na banda LH (f1) as fases são invertidas."]
        self.saida.escrever("\n".join(linhas))

    @protegido
    def publicar(self):
        v = self.form.get()
        lay = self._layout(v)
        if lay is None:
            raise ValueError("Sem layout para esta combinação (use a implementação em microstrip).")
        self.previa.mostrar(lay)
        self.estado.publicar_layout(lay)
        self.estado.msg(f"Layout '{lay.nome}' enviado para Exportar (Ctrl+E).")

    def exportar(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s, v = self._ultimo
            salvar_csv_s(self, "dispositivo", f, s)
