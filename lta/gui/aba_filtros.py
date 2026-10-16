"""Aba Filtros: passa-faixa, dual-band e células híbridas com CSRR (seção 4.2.3)."""
import numpy as np

from ..core import abcd
from ..core import filtros as ft
from .comum import AbaBase, AreaTexto, Formulario, protegido

MODOS = ("Passa-faixa: inversores + ressonadores (Fig. 4.28)", "Dual-band (eqs. 4.49-4.50)",
         "Híbrido com CSRR, periódico (eqs. 4.38-4.46)")
APROX = ("Chebyshev", "Butterworth")


class AbaFiltros(AbaBase):
    TITULO = "Filtros passa-faixa e dual-band"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("modo", "Filtro", MODOS, ""),
            ("aprox", "Aproximação", APROX, ""),
            ("n", "Ordem", 3, ""),
            ("ripple", "Ripple (Chebyshev)", 0.5, "dB"),
            ("f0", "f0", 1.0, "GHz"),
            ("fbw", "FBW", 10.0, "%"),
            ("f1", "f1 (dual-band)", 0.9, "GHz"),
            ("f2", "f2 (dual-band)", 1.8, "GHz"),
            ("fz", "fz/f0 (híbrido)", 2.0, ""),
            ("z0", "Z0", 50.0, "Ω"),
            ("fmax", "f máx", 3.5, "GHz"),
        ], titulo="Especificação")
        self.form.pack(fill="x")
        self.botoes(("Calcular", self.calcular), ("Ex. híbrido 1 GHz (Tab. 4.2)", self.ex_hibrido),
                    ("Ex. Chebyshev 3ª ordem", self.ex_bp), ("Ex. dual-band GSM", self.ex_dual),
                    ("Exportar S2P", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 22)
        self.saida.pack(fill="both", expand=True)
        self._ultimo = None

    def ex_hibrido(self):
        """Livro, p. 283: ordem 3, f0 = 1 GHz, FBW = 10%, g = 1,521, fz = 2f0 (Tabela 4.2)."""
        self.form.set(dict(modo=MODOS[2], n=3, f0=1.0, fbw=10, fz=2.0, z0=50, fmax=3.5))
        self.calcular()

    def ex_bp(self):
        self.form.set(dict(modo=MODOS[0], aprox="Chebyshev", n=3, ripple=0.5, f0=2.0, fbw=20, fmax=4))
        self.calcular()

    def ex_dual(self):
        self.form.set(dict(modo=MODOS[1], aprox="Chebyshev", n=3, ripple=0.1, f1=0.9, f2=1.8, fbw=10, fmax=2.6))
        self.calcular()

    @protegido
    def calcular(self):
        v = self.form.get()
        if v["n"] < 1 or v["fbw"] <= 0 or v["fmax"] <= 0:
            raise ValueError("Ordem ≥ 1, FBW > 0 e f máx > 0.")
        fbw, z0 = v["fbw"] / 100, v["z0"]
        f = np.linspace(v["fmax"] * 1e9 / 4000, v["fmax"] * 1e9, 4000)
        linhas = [v["modo"], ""]
        if v["modo"] == MODOS[2]:
            g = ft.g_periodico(v["n"])
            p = ft.celula_hibrida(v["f0"] * 1e9, fbw, g, v["fz"] * v["f0"] * 1e9, z0)
            m = abcd.potencia(ft.abcd_hibrida(f, p), v["n"])
            s = abcd.abcd_para_s(m, z0)
            s1 = abcd.abcd_para_s(ft.abcd_hibrida(f, p), z0)
            linhas += [f"g de {v['n']} elementos iguais com 3 dB em Ω = 1: g = {g:.4f} (livro: 1,521 para n = 3)",
                       f"Δ = 2·FBW·ω0/g → ressonadores com 3 dB de {2 * fbw / g * v['f0'] * 1e3:.1f} MHz",
                       "Célula (Fig. 3.60b, L desprezada):",
                       f"  Cg = {p['Cg'] * 1e12:.4g} pF (eq. 4.46)   Ld = {p['Ld'] * 1e9:.4g} nH",
                       f"  C = {p['C'] * 1e12:.4g} pF   Lc = {p['Lc'] * 1e9:.4g} nH   Cc = {p['Cc'] * 1e12:.4g} pF",
                       "Livro (Tabela 4.2): Cg 1,59 / Ld 1,33 / C 12,33 / Cc 21,7 pF / Lc 0,23 nH.",
                       "Os valores da tabela dão fz = 1,80 GHz e Zp(f0) = j33,7 Ω, e não as",
                       "condições declaradas; devem incluir o ajuste de layout."]
            extra = [("1 célula", s1)]
        else:
            g = ft.prototipo(v["aprox"], v["n"], v["ripple"])
            linhas.append("g = " + ", ".join(f"{x:.4f}" for x in g))
            if v["modo"] == MODOS[0]:
                f0 = v["f0"] * 1e9
                ceq, leq, delta = ft.ressonadores(g, f0, fbw, z0)
                ys = [ft.tanque(c, l) for c, l in zip(ceq, leq)]
                linhas.append(f"Inversores K = {z0:g} Ω (ideais) e tanques em derivação (eqs. 4.40-4.42):")
                for i, (c, l, d) in enumerate(zip(ceq, leq, delta), 1):
                    linhas.append(f"  {i}: Ceq = {c * 1e12:.4g} pF, Leq = {l * 1e9:.4g} nH, Δ = {d / 2e6 / np.pi:.1f} MHz")
            else:
                f1, f2 = v["f1"] * 1e9, v["f2"] * 1e9
                if not 0 < f1 < f2:
                    raise ValueError("É preciso 0 < f1 < f2.")
                f0 = np.sqrt(f1 * f2)
                ceq, leq, _ = ft.ressonadores(g, f0, fbw, z0)
                ys = []
                linhas += [f"Passa-faixa intermediário em f0 = √(f1·f2) = {f0 / 1e9:.4f} GHz, FBW = {v['fbw']:g}%",
                           "Ressonadores dual-band (L'‖C')‖(L''+C''), eqs. 4.50:"]
                for i, (c, l) in enumerate(zip(ceq, leq), 1):
                    cp, lp, ls, cs = ft.ressonador_dual(c, l, f1, f2)
                    ys.append(ft.admitancia_dual(cp, lp, ls, cs))
                    linhas.append(f"  {i}: C' = {cp * 1e12:.4g} pF, L' = {lp * 1e9:.4g} nH, "
                                  f"L'' = {ls * 1e9:.4g} nH, C'' = {cs * 1e12:.4g} pF")
            s = ft.resposta_inversores(f, g, ys, z0)
            extra = []
        self._ultimo = (f, s, z0)
        self.estado.publicar_resultado("filtro", f, s, z0)
        fg = f / 1e9
        ax = self.grafico.eixos(1, 2)
        ax[0].plot(fg, abcd.db(s[:, 1, 0]), label="|S21|")
        ax[0].plot(fg, abcd.db(s[:, 0, 0]), label="|S11|")
        for nome, sx in extra:
            ax[0].plot(fg, abcd.db(sx[:, 1, 0]), ":", label=f"|S21| {nome}")
        ax[0].set_ylim(-70, 2)
        ax[0].set_title("Resposta (circuito ideal)", fontsize=9)
        ax[1].plot(fg, abcd.db(s[:, 1, 0]), label="|S21|")
        ax[1].plot(fg, abcd.db(s[:, 0, 0]), label="|S11|")
        ax[1].set_ylim(-30, 1)
        ax[1].set_title("Detalhe da(s) banda(s)", fontsize=9)
        for a in ax:
            a.set_xlabel("f [GHz]")
            a.set_ylabel("|S| [dB]")
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        s21 = np.nan_to_num(abcd.db(s[:, 1, 0]), nan=-np.inf)
        ok = s21 > s21.max() - 3.0103
        idx = np.flatnonzero(ok)
        linhas += ["", "Bandas de 3 dB obtidas:"]
        for grupo in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
            if grupo.size:
                a, b = f[grupo[0]], f[grupo[-1]]
                linhas.append(f"  {a / 1e9:.4f} a {b / 1e9:.4f} GHz (centro {np.sqrt(a * b) / 1e9:.4f} GHz, "
                              f"{100 * (b - a) / np.sqrt(a * b):.1f}%)")
        linhas += ["", "Inversores ideais. Na prática são células CRLH de",
                   "-90° (livro, p. 281) ou linhas de λ/4, que limitam a",
                   "banda e criam espúrios: simule a implementação."]
        self.saida.escrever("\n".join(linhas))

    def exportar_s2p(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s, z0 = self._ultimo
            self.salvar_s2p("filtro", f, s, z0, "Filtro (circuito ideal, cap. 4)")
