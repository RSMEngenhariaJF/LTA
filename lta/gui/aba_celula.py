"""Aba Célula ABCD: dispersão e impedância de Bloch de células concentradas (seção 2.3)."""
import numpy as np

from ..core import abcd
from .comum import AbaBase, AreaTexto, Formulario, protegido
from .previa import PreviaCelula


class AbaCelula(AbaBase):
    TITULO = "Célula ABCD"
    PREVIA = PreviaCelula

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("topologia", "Topologia", abcd.TOPOLOGIAS, ""),
            ("comparar", "Comparar com", ("Nenhuma",) + abcd.TOPOLOGIAS, ""),
            ("ls", "Ls (série)", 2.5, "nH"),
            ("cs", "Cs (série)", 0.0, "pF"),
            ("lp", "Lp (paralelo)", 0.0, "nH"),
            ("cp", "Cp (paralelo)", 1.0, "pF"),
            ("n", "Nº de células", 20, ""),
            ("zp", "Z das portas", 50.0, "Ω"),
            ("fmin", "f mín", 0.01, "GHz"),
            ("fmax", "f máx", 8.0, "GHz"),
            ("npts", "Pontos", 1001, ""),
            ("s11max", "Banda útil: |S11| ≤", -15.0, "dB"),
            ("s21min", "Banda útil: |S21| ≥", -1.0, "dB"),
        ], titulo="Célula unitária (0 = elemento ausente)")
        self.form.pack(fill="x")
        self.botoes(("Calcular", self.calcular), ("Exemplo Fig. 2.8", self.exemplo),
                    ("Ex. passa-tudo (treliça)", self.exemplo_trelica), ("Ex. T em ponte", self.exemplo_ponte),
                    ("Exportar S2P", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 12)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Ramo série: Ls em série com Cs (Cs = 0 → curto).\n"
            "Ramo paralelo: Cp em paralelo com Lp (Lp = 0 → aberto).\n\n"
            "Exemplos:\n• Linha convencional: Ls, Cp\n• Linha LH: Cs, Lp\n• CRLH: Ls, Cs, Lp, Cp\n\n"
            "Banda larga (S21 e S11 planos): topologias passa-tudo\n"
            "• Treliça (X): braços = ramo série, diagonais =\n  ramo paralelo; ZB = √(Zs·Zp) constante.\n"
            "• T em ponte: Z1 = ramo série, Z2 = ramo paralelo,\n  ponte Z3 = 4·Z2; ZB = √(2·Z1·Z2).\n"
            "(Martín 2015, seção 6.4, pp. 429-439)\n\n"
            "Dispersão: cosh(γl) = (A+D)/2 (eq. 2.21)\nZB: eq. 2.26")
        self._ultimo = None
        self.observar_previa(self.form)
        self.atualizar_previa()

    def _elementos(self, v):
        el = dict(ls=v["ls"] * 1e-9, cs=v["cs"] * 1e-12, lp=v["lp"] * 1e-9, cp=v["cp"] * 1e-12)
        if any(x < 0 for x in el.values()):
            raise ValueError("Valores de elementos não podem ser negativos.")
        if not any(el.values()):
            raise ValueError("Defina pelo menos um elemento.")
        if {v["topologia"], v.get("comparar")} & {"Treliça (X)", "T em ponte"}:
            if not (el["ls"] or el["cs"]) or not (el["lp"] or el["cp"]):
                raise ValueError("Treliça e T em ponte precisam dos dois ramos (série e paralelo).")
        return el

    def atualizar_previa(self):
        try:
            v = self.form.get()
            self.previa.mostrar(v["topologia"], v["n"], **self._elementos(v))
        except ValueError as e:
            self.previa.mostrar_erro(str(e))

    def exemplo(self):
        self.form.set(dict(topologia="T", comparar="L", ls=2.5, cs=0, lp=0, cp=1.0, n=20, zp=50,
                           fmin=0.01, fmax=8, npts=1001))
        self.calcular()

    def exemplo_trelica(self):
        """Mesmos L e C da Fig. 2.8, em treliça: passa-tudo (Martín, Fig. 6.25a)."""
        self.form.set(dict(topologia="Treliça (X)", comparar="T", ls=2.5, cs=0, lp=0, cp=1.0, n=20, zp=50,
                           fmin=0.01, fmax=20, npts=2001))
        self.calcular()

    def exemplo_ponte(self):
        """T em ponte com Z1 = 1/(jwC), Z2 = jwL, Z3 = 4Z2 (Martín, Fig. 6.31): ZB = sqrt(2L/C) = 50 ohm."""
        self.form.set(dict(topologia="T em ponte", comparar="T", ls=0, cs=1.0, lp=1.25, cp=0, n=10, zp=50,
                           fmin=0.01, fmax=20, npts=2001))
        self.calcular()

    @protegido
    def calcular(self):
        v = self.form.get()
        if v["fmin"] <= 0 or v["fmax"] <= v["fmin"]:
            raise ValueError("Faixa de frequência inválida.")
        f = np.linspace(v["fmin"], v["fmax"], v["npts"]) * 1e9
        el = self._elementos(v)
        self.previa.mostrar(v["topologia"], v["n"], **el)
        tops = [v["topologia"]] + ([v["comparar"]] if v["comparar"] not in ("Nenhuma", v["topologia"]) else [])

        ax = self.grafico.eixos(2, 2)
        texto = []
        for i, top in enumerate(tops):
            m = abcd.celula(f, top, **el)
            al, bl = abcd.dispersao(m)
            zb, _ = abcd.impedancia_bloch(m)
            s = abcd.abcd_para_s(abcd.potencia(m, v["n"]), v["zp"])
            ls = "-" if i == 0 else "--"
            ax[0].plot(bl / np.pi, f / 1e9, ls, color=f"C{i}", label=f"{top}: βl")
            ax[0].plot(al / np.pi, f / 1e9, ls, color=f"C{i}", alpha=0.4, label=f"{top}: αl")
            ax[1].plot(f / 1e9, zb.real, ls, color=f"C{i}", label=f"{top}: Re ZB")
            ax[1].plot(f / 1e9, zb.imag, ls, color=f"C{i}", alpha=0.5, label=f"{top}: Im ZB")
            ax[2].plot(f / 1e9, abcd.db(s[:, 1, 0]), ls, color=f"C{i}", label=f"{top}: S21")
            ax[3].plot(f / 1e9, abcd.db(s[:, 0, 0]), ls, color=f"C{i}", label=f"{top}: S11")
            if i == 0:
                self._ultimo = (f, s, v["zp"], top, v)
            passa = al < 1e-6
            bordas = f[np.flatnonzero(np.diff(passa.astype(int)))] / 1e9
            texto.append(f"[{top}] ZB em {f[0] / 1e9:g} GHz: {zb[0].real:.2f} {zb[0].imag:+.2f}j Ω")
            texto.append(f"  Bordas de banda: {', '.join(f'{b:.3f}' for b in bordas) or 'nenhuma'} GHz")
            passa_zb = zb[al < 1e-6].real
            if passa_zb.size:
                texto.append(f"  Re ZB na banda passante: {passa_zb.min():.2f} a {passa_zb.max():.2f} Ω")
            texto.append("  " + abcd.descrever_banda(f, s, v["s11max"], v["s21min"]))
            if top == "Treliça (X)":
                texto.append("  Rede balanceada (diferencial): requer linhas acopladas,\n"
                             "  duas camadas ou vias (Martín, pp. 434-437).")
        ax[0].set_xlabel("βl/π, αl/π")
        ax[0].set_ylabel("f [GHz]")
        ax[0].set_title("Diagrama de dispersão", fontsize=9)
        ax[1].set_ylabel("ZB [Ω]")
        ax[1].set_title("Impedância de Bloch", fontsize=9)
        ax[2].set_ylabel("|S21| [dB]")
        ax[2].set_title(f"{v['n']} células, portas {v['zp']:g} Ω", fontsize=9)
        ax[3].set_ylabel("|S11| [dB]")
        for a in ax:
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        for a in ax[1:]:
            a.set_xlabel("f [GHz]")
        ax[2].set_ylim(bottom=max(ax[2].get_ylim()[0], -80))
        ax[3].set_ylim(bottom=max(ax[3].get_ylim()[0], -60))
        self.grafico.desenhar()
        self.ver_resultados()
        self.saida.escrever("\n".join(texto))
        self.estado.publicar_resultado("celula_abcd", f, self._ultimo[1], v["zp"])

    def exportar_s2p(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s, zp, top, v = self._ultimo
            self.salvar_s2p(f"celula_{top}_{v['n']}cel", f, s, zp,
                            f"Celula {top}: Ls={v['ls']}nH Cs={v['cs']}pF Lp={v['lp']}nH Cp={v['cp']}pF N={v['n']}")
