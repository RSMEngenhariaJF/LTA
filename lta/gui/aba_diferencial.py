"""Aba Linha diferencial com supressão do modo comum (seção 6.3)."""
import numpy as np

from ..core import abcd
from ..core import diferencial as df
from .comum import AbaBase, AreaTexto, Formulario, protegido

PRESETS = tuple(df.TABELA_6_1) + ("Personalizado",)


class AbaDiferencial(AbaBase):
    TITULO = "Linha diferencial (filtro de modo comum)"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("preset", "Elementos", PRESETS, ""),
            ("Le", "Le (modo par)", 6.3, "nH"),
            ("Ce", "Ce (linha-CSRR)", 1.1, "pF"),
            ("Lc", "Lc (CSRR)", 2.1, "nH"),
            ("Cc", "Cc (CSRR)", 3.2, "pF"),
            ("CR", "CR (entre CSRRs)", 0.05, "pF"),
            ("n", "Nº de células", 6, ""),
            ("lcel", "Período da célula", 7.5, "mm"),
            ("zodd", "Z do modo ímpar", 50.0, "Ω"),
            ("eeff", "εeff do modo ímpar", 6.8, ""),
            ("nmax", "Curva FBW até n =", 10, ""),
        ], titulo="Par diferencial com CSRR/DS-CSRR (Fig. 6.12)")
        self.form.pack(fill="x")
        self.form.ao_mudar("preset", self._carregar)
        self.botoes(("Calcular", self.calcular), ("Exportar S (modo comum)", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 20)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Modo comum (par): o CSRR, com o plano de simetria\n"
            "alinhado ao do par, é excitado e cria um zero de\n"
            "transmissão (eq. 6.11). Modo diferencial (ímpar):\n"
            "o CSRR não é excitado e a linha é transparente.\n\n"
            "Elementos da Tabela 6.1: escolha em 'Elementos'.")
        self._ultimo = None

    def _carregar(self):
        p = self.form.get()["preset"]
        if p in df.TABELA_6_1:
            f0, le, ce, lc, cc, cr, _ = df.TABELA_6_1[p]
            self.form.set(dict(Le=le, Ce=ce, Lc=lc, Cc=cc, CR=cr, lcel=df.PERIODO_FIG_6_13[p]))

    @protegido
    def calcular(self):
        v = self.form.get()
        n_, p = 1e-9, 1e-12
        el = (v["Le"] * n_, v["Ce"] * p, v["Lc"] * n_, v["Cc"] * p, v["CR"] * p)
        if min(el) <= 0 or v["n"] < 1 or v["lcel"] <= 0 or v["zodd"] <= 0 or v["eeff"] < 1:
            raise ValueError("Elementos, nº de células, período e Z ímpar devem ser positivos; εeff ≥ 1.")
        fz = df.zero_comum(*el[:4])
        fa, fb, fbw = df.banda_rejeicao(*el)
        f = np.linspace(0.2 * fz, 2.5 * fz, 3000)
        scc = df.s_comum(f, v["n"], *el)
        lo, co = df.linha_impar(v["zodd"], v["eeff"], v["lcel"])
        sdd = df.s_diferencial(f, v["n"], lo, co)
        self._ultimo = (f, scc)
        self.estado.publicar_resultado("modo_comum", f, scc, 50.0)
        fg = f / 1e9
        ax = self.grafico.eixos(1, 3)
        ax[0].plot(fg, abcd.db(scc[:, 1, 0]), label="|Scc21|")
        ax[0].plot(fg, abcd.db(scc[:, 0, 0]), label="|Scc11|")
        ax[0].plot(fg, abcd.db(sdd[:, 1, 0]), label="|Sdd21|")
        ax[0].axvspan(fa / 1e9, fb / 1e9, color="#f4d6d6", alpha=0.5, label="banda máx. (∞ células)")
        ax[0].set_ylim(-60, 2)
        ax[0].set_ylabel("|S| [dB]")
        ax[0].set_title(f"{v['n']} células (circuito, cf. Figs. 6.11 e 6.14)", fontsize=9)
        c1, c2 = df.dispersao_comum(f, *el)
        for c, nome in ((c1, "modo 1"), (c2, "modo 2")):
            g = np.arccosh(c.astype(complex))
            ax[1].plot(fg, np.abs(g.real), label=f"αl {nome}")
        ax[1].set_ylabel("αl [Np]")
        ax[1].set_ylim(0, 3)
        ax[1].set_title("Atenuação de Bloch do modo comum (eq. 3.96)", fontsize=9)
        ns = np.arange(1, max(2, v["nmax"]) + 1)
        fbw_n = []
        for k in ns:
            s21 = abcd.db(df.s_comum(f, int(k), *el)[:, 1, 0])
            rej = s21 < -20
            i = int(np.argmin(np.abs(f - fz)))
            if rej[i]:
                a = b = i
                while a > 0 and rej[a - 1]:
                    a -= 1
                while b < len(f) - 1 and rej[b + 1]:
                    b += 1
                fbw_n.append(100 * (f[b] - f[a]) / np.sqrt(f[a] * f[b]))
            else:
                fbw_n.append(0.0)
        ax[2].plot(ns, fbw_n, "o-", label="FBW a -20 dB")
        ax[2].axhline(100 * fbw, color="C3", ls="--", label="∞ células")
        ax[2].set_xlabel("Nº de células")
        ax[2].set_ylabel("FBW [%]")
        ax[2].set_title("Rejeição × ordem (cf. Fig. 6.13)", fontsize=9)
        for a in ax[:2]:
            a.set_xlabel("f [GHz]")
        for a in ax:
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        banda = (f >= fa) & (f <= fb)
        linhas = [f"fz (modo comum) = {fz / 1e9:.4f} GHz (eq. 6.11: 1/(2π√(2Lc(Ce + Cc/2))))",
                  f"Banda máxima (estrutura infinita): {fa / 1e9:.3f} a {fb / 1e9:.3f} GHz, FBW = {100 * fbw:.1f}%"]
        if v["preset"] in df.TABELA_6_1:
            linhas.append(f"  Livro (Tabela 6.1): {df.TABELA_6_1[v['preset']][6]:.1f}%")
        linhas += ["", " n   FBW a -20 dB"] + [f" {k:2d}  {x:6.1f}%" for k, x in zip(ns, fbw_n)]
        linhas += ["", f"Modo diferencial: |Sdd21| na banda ≥ {abcd.db(sdd[banda, 1, 0]).min():.2f} dB "
                   f"(Lo = {lo * 1e9:.3g} nH, Co = {co * 1e12:.3g} pF por célula de {v['lcel']:g} mm, "
                   f"Z = {v['zodd']:g} Ω)", "",
                   "Ressonadores com fendas estreitas (c = d = 0,2 mm)",
                   "mantêm o modelo válido; o livro extrai Le, Ce, Lc e",
                   "Cc de 1 célula e ajusta CR com 2 células (Ap. G)."]
        self.saida.escrever("\n".join(linhas))

    def exportar_s2p(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s = self._ultimo
            self.salvar_s2p("modo_comum", f, s, 50.0, "Modo comum (meio circuito) de linha diferencial com CSRR")
