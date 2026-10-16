"""Aba Híbridas e divisores: dispositivos convencionais e com linhas artificiais (seções 4.2.1-4.2.2)."""
from pathlib import Path
from tkinter import filedialog

import numpy as np

from ..core import abcd
from ..core import dispositivos as dv
from .comum import PASTA_SAIDAS, AbaBase, AreaTexto, Formulario, protegido, relativo

DISPOSITIVOS = ("Branch-line (híbrida 90°)", "Rat-race (híbrida 180°)", "Wilkinson", "Divisor com inversor (Y)")
IMPL = ("Dual-band CRLH (-90° em f1, +90° em f2)", "Banda larga: rat-race com -90° LH no lugar de 270°")


class AbaDispositivos(AbaBase):
    TITULO = "Híbridas e divisores"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("disp", "Dispositivo", DISPOSITIVOS, ""),
            ("impl", "Linhas artificiais", IMPL, ""),
            ("z0", "Z0 das portas", 50.0, "Ω"),
            ("f0", "f0 (versão convencional)", 1.6, "GHz"),
            ("f1", "f1 (banda LH)", 0.9, "GHz"),
            ("f2", "f2 (banda RH)", 1.8, "GHz"),
            ("nlh", "Células da linha LH (rat-race)", 1, ""),
            ("fmax", "f máx", 3.0, "GHz"),
        ], titulo="Dispositivo")
        self.form.pack(fill="x")
        self.botoes(("Calcular", self.calcular), ("Ex. branch-line GSM (Fig. 4.16)", self.ex_branch),
                    ("Ex. divisor 2,4/3,75 GHz (Fig. 4.11)", self.ex_divisor),
                    ("Ex. rat-race banda larga (Fig. 4.5)", self.ex_ratrace), ("Exportar S (.csv)", self.exportar))
        self.saida = AreaTexto(self.painel, 22)
        self.saida.pack(fill="both", expand=True)
        self._ultimo = None

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

    # ------------------------------------------------------------------ montagem
    def _montar(self, v, f, artificial):
        """S do dispositivo; artificial=False dá a versão convencional em f0."""
        z0, f0 = v["z0"], v["f0"] * 1e9
        disp = v["disp"]
        z_ramo = {DISPOSITIVOS[0]: (z0 / np.sqrt(2), z0), DISPOSITIVOS[1]: (z0 * np.sqrt(2),) * 2,
                  DISPOSITIVOS[2]: (z0 * np.sqrt(2),) * 2, DISPOSITIVOS[3]: (z0 / np.sqrt(2),) * 2}[disp]
        info = []
        if not artificial:
            l = [dv.linha_convencional(z, np.pi / 2, f0) for z in z_ramo]
            l270 = dv.linha_convencional(z_ramo[0], 3 * np.pi / 2, f0)
        elif v["impl"] == IMPL[0]:
            if disp == DISPOSITIVOS[1]:
                raise ValueError("Rat-race dual-band exige ramo de 270° nas duas bandas; use a opção "
                                 "de banda larga para o rat-race.")
            f1, f2 = v["f1"] * 1e9, v["f2"] * 1e9
            l, l270 = [], None
            for z in z_ramo:
                func, c = dv.linha_crlh_dual(z, f1, f2)
                l.append(func)
                if len(l) == 1 or z != z_ramo[0]:
                    info.append(f"  Z = {z:.2f} Ω: LR = {c.LR * 1e9:.4g} nH, CR = {c.CR * 1e12:.4g} pF, "
                                f"LL = {c.LL * 1e9:.4g} nH, CL = {c.CL * 1e12:.4g} pF")
        else:
            if disp != DISPOSITIVOS[1]:
                raise ValueError("A opção de banda larga (linha LH de -90°) se aplica ao rat-race.")
            l = [dv.linha_convencional(z, np.pi / 2, f0) for z in z_ramo]
            l270, (cl, ll) = dv.linha_plh(z_ramo[0], -np.pi / 2, f0, max(1, v["nlh"]))
            info.append(f"  Linha LH de -90° em f0: {max(1, v['nlh'])} célula(s) T, CL = {cl * 1e12:.4g} pF, "
                        f"LL = {ll * 1e9:.4g} nH (Z = {z_ramo[0]:.2f} Ω)")
        if disp == DISPOSITIVOS[0]:
            s = dv.branch_line(f, l[0], l[1], z0)
        elif disp == DISPOSITIVOS[1]:
            s = dv.rat_race(f, l[0], l270, z0)
        elif disp == DISPOSITIVOS[2]:
            s = dv.wilkinson(f, l[0], z0)
        else:
            s = dv.divisor_inversor(f, l[0], z0)
        return s, info

    def _metricas(self, v, f, s, nome, f_proj):
        """Faixas em torno de cada frequência de projeto (critérios da Tabela 4.1)."""
        disp = v["disp"]
        s11 = abcd.db(s[:, 0, 0])
        a, b = abcd.db(s[:, 1, 0]), abcd.db(s[:, 2, 0])
        linhas = [f"[{nome}]"]
        crit = {"|S11| < -10 dB": s11 < -10,
                "saídas -3 ± 1 dB": (np.abs(a + 3.01) < 1) & (np.abs(b + 3.01) < 1)}
        if disp == DISPOSITIVOS[0]:
            crit["fase S21-S31 = ±90° ± 5°"] = np.abs(np.abs(dv.diferenca_fase(s[:, 1, 0], s[:, 2, 0])) - 90) < 5
            crit["isolamento |S41| < -15 dB"] = abcd.db(s[:, 3, 0]) < -15
        elif disp == DISPOSITIVOS[1]:
            crit["Δ: fase S21-S31 = 180° ± 5°"] = np.abs(np.abs(dv.diferenca_fase(s[:, 1, 0], s[:, 2, 0])) - 180) < 5
            crit["Σ: fase S24-S34 = 0° ± 5°"] = np.abs(dv.diferenca_fase(s[:, 1, 3], s[:, 2, 3])) < 5
            crit["isolamento |S41| < -15 dB"] = abcd.db(s[:, 3, 0]) < -15
        elif disp == DISPOSITIVOS[2]:
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
        s_art, info = self._montar(v, f, True)
        s_conv, _ = self._montar(v, f, False)
        self._ultimo = (f, s_art, v)
        disp = v["disp"]
        fg = f / 1e9
        ax = self.grafico.eixos(2, 2)
        n_portas = s_art.shape[1]
        for j in range(n_portas):
            ax[0].plot(fg, abcd.db(s_art[:, j, 0]), label=f"|S{j + 1}1|")
            ax[1].plot(fg, abcd.db(s_conv[:, j, 0]), ls="--", label=f"|S{j + 1}1|")
        ax[0].set_title("Com linhas artificiais", fontsize=9)
        ax[1].set_title(f"Convencional em f0 = {v['f0']:g} GHz", fontsize=9)
        for a in ax[:2]:
            a.set_ylim(-40, 1)
            a.set_ylabel("|S| [dB]")
        if disp == DISPOSITIVOS[1]:
            pares = (("Δ: ∠S21 - ∠S31", (1, 0), (2, 0), 180), ("Σ: ∠S24 - ∠S34", (1, 3), (2, 3), 0))
        else:
            pares = (("∠S21 - ∠S31", (1, 0), (2, 0), 90 if disp == DISPOSITIVOS[0] else 0),)
        for nome, (i1, j1), (i2, j2), _ in pares:
            ax[2].plot(fg, dv.diferenca_fase(s_art[:, i1, j1], s_art[:, i2, j2]), label=f"{nome} artificial")
            ax[2].plot(fg, dv.diferenca_fase(s_conv[:, i1, j1], s_conv[:, i2, j2]), ls="--",
                       label=f"{nome} convencional")
        ax[2].set_ylim(-190, 190)
        ax[2].set_ylabel("Balanço de fase [°]")
        ax[2].set_title("Balanço de fase entre as saídas", fontsize=9)
        ax[3].plot(fg, np.abs(s_art[:, 1, 0]) ** 2 + np.abs(s_art[:, 2, 0]) ** 2, label="|S21|²+|S31|² artificial")
        ax[3].plot(fg, np.abs(s_conv[:, 1, 0]) ** 2 + np.abs(s_conv[:, 2, 0]) ** 2, ls="--",
                   label="|S21|²+|S31|² convencional")
        ax[3].set_ylim(0, 1.05)
        ax[3].set_title("Potência entregue às saídas", fontsize=9)
        for a in ax:
            a.set_xlabel("f [GHz]")
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        if v["impl"] == IMPL[0]:
            f_proj = [v["f1"] * 1e9, v["f2"] * 1e9]
        else:
            f_proj = [v["f0"] * 1e9]
        linhas = [disp, f"Linhas artificiais: {v['impl']}"] + info + [""]
        linhas += self._metricas(v, f, s_art, "artificial", f_proj) + [""]
        linhas += self._metricas(v, f, s_conv, f"convencional, f0 = {v['f0']:g} GHz", [v["f0"] * 1e9])
        linhas += ["", "Circuito ideal sem perdas. Na banda LH (f1) as fases",
                   "são invertidas: na branch-line, a saída acoplada passa",
                   "a adiantar 90° (livro, Fig. 4.16: -93° e +92°)."]
        self.saida.escrever("\n".join(linhas))

    def exportar(self):
        if not self._ultimo:
            self.calcular()
        if not self._ultimo:
            return
        f, s, v = self._ultimo
        PASTA_SAIDAS.mkdir(exist_ok=True)
        caminho = filedialog.asksaveasfilename(initialdir=PASTA_SAIDAS, initialfile="dispositivo_S.csv",
                                               defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if caminho:
            n = s.shape[1]
            cab = ["f_GHz"] + [f"S{i + 1}{j + 1}_{q}" for i in range(n) for j in range(n) for q in ("dB", "graus")]
            linhas = [";".join(cab)]
            for fi, si in zip(f, s):
                vals = [f"{fi / 1e9:.9g}"]
                for i in range(n):
                    for j in range(n):
                        vals += [f"{abcd.db(si[i, j]):.6f}", f"{np.degrees(np.angle(si[i, j])):.4f}"]
                linhas.append(";".join(vals))
            Path(caminho).write_text("\n".join(linhas) + "\n", encoding="ascii")
            self.estado.msg(f"S do dispositivo salvo em {relativo(caminho)}")
