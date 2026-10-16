"""Abas SIW (seção 6.5) e linha com capacitores negativos (não Foster, seção 6.4.2)."""
import numpy as np

from ..core import abcd
from ..core import naofoster as nf
from ..core import siw
from ..core.microstrip import C0, sintese
from ..geometria import construtores
from .comum import AbaBase, AreaTexto, Formulario, protegido
from .previa import PreviaEstrutura


class AbaSIW(AbaBase):
    TITULO = "Guia integrado ao substrato (SIW)"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("fc", "Corte do TE10", 6.0, "GHz"),
            ("d", "Diâmetro das vias", 0.8, "mm"),
            ("s", "Passo das vias (0 = 2d)", 1.4, "mm"),
            ("comp", "Comprimento do guia", 30.0, "mm"),
            ("wc", "Cunha: largura no guia (0 = 0,4·W)", 0.0, "mm"),
            ("lc", "Cunha: comprimento (0 = λg/4 em 1,5·fc)", 0.0, "mm"),
            ("lalim", "Alimentação", 5.0, "mm"),
        ], titulo="SIW (eqs. 6.31-6.34)")
        self.form.pack(fill="x")
        self.botoes(("Calcular", self.calcular), ("Enviar p/ exportação", self.publicar))
        self.saida = AreaTexto(self.painel, 20)
        self.saida.pack(fill="both", expand=True)
        self.observar_previa(self.form)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    def _gerar(self):
        v = self.form.get()
        sub = self.estado.substrato
        s = v["s"] * 1e-3 if v["s"] > 0 else None
        p = siw.projeto(v["fc"] * 1e9, sub.er, v["d"] * 1e-3, s)
        w50 = sintese(50, sub)
        wc = v["wc"] if v["wc"] > 0 else max(0.4 * p.w * 1e3, w50 * 1.01)
        lc = v["lc"] if v["lc"] > 0 else p.lambda_g(1.5 * p.fc10) / 4 * 1e3
        info = {"tipo": "SIW", "fc10_GHz": v["fc"], "W_mm": f"{p.w * 1e3:.4f}", "a_mm": f"{p.a * 1e3:.4f}",
                "d_mm": v["d"], "s_mm": f"{p.s * 1e3:.3f}", "comprimento_mm": v["comp"],
                "cunha_mm": f"{wc:.3f} x {lc:.3f}"}
        lay = construtores.linha_siw(p, v["comp"], sub, w50, wc, lc, v["lalim"], info=info)
        return v, sub, p, w50, wc, lc, lay

    def atualizar_previa(self):
        try:
            *_, lay = self._gerar()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        self.previa.mostrar(lay)

    @protegido
    def calcular(self):
        v, sub, p, w50, wc, lc, lay = self._gerar()
        self.previa.mostrar(lay)
        f = np.linspace(0.5 * p.fc10, 2.4 * p.fc10, 2000)
        comp = v["comp"] * 1e-3
        s = abcd.abcd_para_s(siw.abcd_siw(f, p, comp, sub.tand), abs(p.z_te(np.array([1.5 * p.fc10]))[0]))
        acima = f > p.fc10
        ax = self.grafico.eixos(1, 3)
        ax[0].plot(p.beta(f).real[acima] * p.a / np.pi, f[acima] / 1e9, label="TE10")
        k = 2 * np.pi * f * np.sqrt(p.er) / C0
        ax[0].plot(k * p.a / np.pi, f / 1e9, ":", color="gray", label="TEM no dielétrico")
        ax[0].axhline(p.fc20 / 1e9, color="C3", ls="--", lw=0.8, label="corte TE20")
        ax[0].set_xlabel("β·a/π")
        ax[0].set_ylabel("f [GHz]")
        ax[0].set_title("Dispersão do TE10", fontsize=9)
        ax[1].plot(f[acima] / 1e9, p.lambda_g(f[acima]) * 1e3, label="λg")
        ax[1].set_ylim(0, 6 * p.a * 1e3)
        ax2 = ax[1].twinx()
        ax2.plot(f[acima] / 1e9, np.abs(p.z_te(f[acima])), color="C1")
        ax2.set_ylabel("Z_TE [Ω]", color="C1")
        ax2.set_ylim(0, 3 * ETA(p))
        ax[1].set_xlabel("f [GHz]")
        ax[1].set_ylabel("λg [mm]")
        ax[1].set_title("Comprimento de onda e impedância de onda", fontsize=9)
        ax[2].plot(f / 1e9, abcd.db(s[:, 1, 0]), label=f"|S21|, {v['comp']:g} mm, tan δ = {sub.tand:g}")
        ax[2].axvspan(p.fc10 / 1e9, p.fc20 / 1e9, color="#e3f0e3", label="banda monomodo")
        ax[2].set_ylim(-40, 1)
        ax[2].set_xlabel("f [GHz]")
        ax[2].set_title("Guia (sem transições nem perdas de condutor)", fontsize=9)
        for a in ax:
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        fm = 1.5 * p.fc10
        av = siw.verificacoes(p)
        self.saida.escrever("\n".join([
            f"{sub.nome}: εr = {sub.er:g}, h = {sub.h_mm:g} mm",
            f"Largura efetiva a = c/(2fc√εr) = {p.a * 1e3:.3f} mm (eq. 6.34)",
            f"Distância entre fileiras W = {p.w * 1e3:.3f} mm (eq. 6.31)",
            f"Vias: d = {p.d * 1e3:g} mm (máx. {p.d_max * 1e3:.2f}, eq. 6.32), s = {p.s * 1e3:.3f} mm "
            f"(máx. 2d = {p.s_max * 1e3:.2f}, eq. 6.33)",
            f"Banda monomodo: {p.fc10 / 1e9:.3f} a {p.fc20 / 1e9:.3f} GHz (TE10 a TE20)",
            f"Em 1,5·fc = {fm / 1e9:.2f} GHz: λg = {p.lambda_g(fm) * 1e3:.2f} mm, "
            f"Z_TE = {abs(p.z_te(np.array([fm]))[0]):.1f} Ω, αd = {p.alfa_dieletrico(fm, sub.tand) * 8.686 / 100:.4f} dB/cm",
            "", f"Transição: linha de 50 Ω (W = {w50:.3f} mm) → cunha {wc:.3f} mm × {lc:.2f} mm",
            "(dimensões iniciais; otimizar em EM, livro p. 465)"] + [f"⚠ {a}" for a in av]))

    @protegido
    def publicar(self):
        *_, lay = self._gerar()
        self.previa.mostrar(lay)
        self.estado.publicar_layout(lay)
        self.estado.msg(f"Layout '{lay.nome}' enviado para Exportar (Ctrl+E).")


def ETA(p):
    """Impedância intrínseca do dielétrico (escala do gráfico de Z_TE)."""
    return 376.73 / np.sqrt(p.er)


class AbaNaoFoster(AbaBase):
    TITULO = "Linha com capacitores negativos (não Foster)"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("z", "Z da linha hospedeira", 31.5, "Ω"),
            ("c", "C' da linha hospedeira", 106.0, "pF/m"),
            ("cn", "|Cn| (capacitor negativo)", 6.1, "pF"),
            ("passo", "Espaçamento", 0.1, "m"),
            ("n", "Nº de células", 10, ""),
            ("fmax", "f máx", 200.0, "MHz"),
        ], titulo="Linha ENZ (Fig. 6.32)")
        self.form.pack(fill="x")
        self.botoes(("Calcular", self.calcular), ("Ex. Fig. 6.32 (45 pF/m)", self.exemplo))
        self.saida = AreaTexto(self.painel, 18)
        self.saida.pack(fill="both", expand=True)

    def exemplo(self):
        """Livro, Fig. 6.32: hospedeira 31,5 Ω, 106 pF/m → 45 pF/m carregada."""
        self.form.set(dict(z=31.5, c=106, cn=6.1, passo=0.1, n=10, fmax=200))
        self.calcular()

    @protegido
    def calcular(self):
        v = self.form.get()
        lin = nf.LinhaENZ(v["z"], v["c"] * 1e-12, v["cn"] * 1e-12, v["passo"])
        if not lin.estavel:
            raise ValueError(f"C' efetiva = {lin.c_ef * 1e12:.1f} pF/m ≤ 0: estrutura instável (livro, p. 462).")
        f = np.linspace(v["fmax"] * 1e6 / 2000, v["fmax"] * 1e6, 2000)
        s = abcd.abcd_para_s(lin.abcd(f, v["n"]), lin.z_ef)
        comp = v["n"] * v["passo"]
        tg = abcd.atraso_grupo(f, s[:, 1, 0])
        fase = -np.unwrap(np.angle(s[:, 1, 0]))
        vp = 2 * np.pi * f * comp / fase
        ax = self.grafico.eixos(1, 2)
        ax[0].plot(f / 1e6, comp / tg / C0, label="vg/c")
        ax[0].plot(f / 1e6, vp / C0, "--", label="vp/c")
        ax[0].axhline(1, color="gray", ls=":", lw=0.8)
        ax[0].set_ylim(0, 2.5)
        ax[0].set_xlabel("f [MHz]")
        ax[0].set_title("Velocidades (cf. Fig. 6.32b)", fontsize=9)
        ax[1].plot(f / 1e6, abcd.db(s[:, 0, 0]), label="|S11|")
        ax[1].plot(f / 1e6, abcd.db(s[:, 1, 0]), label="|S21|")
        ax[1].set_ylim(-60, 2)
        ax[1].set_xlabel("f [MHz]")
        ax[1].set_title(f"{v['n']} células, portas {lin.z_ef:.1f} Ω", fontsize=9)
        for a in ax:
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        self.saida.escrever("\n".join([
            f"L' = Z²·C' = {lin.l_por_m * 1e9:.1f} nH/m (hospedeira)",
            f"C' efetiva = {lin.c_ef * 1e12:.1f} pF/m (estável: > 0)",
            f"Z efetiva = {lin.z_ef:.1f} Ω   ε equivalente = c²L'C' = {lin.eps_ef:.3f}",
            f"vp = vg = {lin.v_sobre_c:.2f}·c (sem dispersão no limite contínuo)", "",
            "Livro (Fig. 6.32): 31,5 Ω e 106 pF/m carregada para 45 pF/m",
            "com -60 pF dá 50 Ω e 'ε equivalente 0,3'. Pelas mesmas",
            "relações a conta dá 48,3 Ω e ε = 0,43 (registrado).", "",
            "Capacitores negativos ideais (NIC sem limite de banda);",
            "na prática a banda é limitada pelo conversor de",
            "impedância negativa e a estabilidade exige C' > 0."]))
