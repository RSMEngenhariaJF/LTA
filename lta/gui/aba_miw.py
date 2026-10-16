"""Aba Linhas de atraso magnetoindutivas e eletroindutivas (seção 6.2)."""
import numpy as np

from ..core import abcd
from ..core import miw
from ..core.microstrip import C0
from .comum import AbaBase, AreaTexto, Formulario, protegido

TIPOS = ("MIW coplanar (M < 0, regressiva)", "MIW axial (M > 0, progressiva)", "EIW com CSRRs (regressiva)")


class AbaMIW(AbaBase):
    TITULO = "Linhas de atraso MIW/EIW"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("tipo", "Cadeia", TIPOS, ""),
            ("f0", "f0 (ressonância)", 2.45, "GHz"),
            ("k", "|M|/L ou CM/Cc", 0.0256, ""),
            ("l", "Período", 3.75, "mm"),
            ("L", "L do ressonador (circuito)", 10.0, "nH"),
            ("tau", "Largura do pulso do leitor", 15.0, "ns"),
            ("bits", "Código (bits)", "11", ""),
            ("npb", "Ressonadores por bit (0 = eq. 6.10)", 0, ""),
            ("b", "Refletor: susceptância normalizada", 1.0, ""),
        ], titulo="Cadeia de ressonadores acoplados")
        self.form.pack(fill="x")
        self.botoes(("Calcular", self.calcular), ("Ex. etiqueta 2,45 GHz (p. 429)", self.exemplo))
        self.saida = AreaTexto(self.painel, 20)
        self.saida.pack(fill="both", expand=True)

    def exemplo(self):
        """Livro, p. 429: f0 = 2,45 GHz, l = 3,75 mm, banda de 5,13% (|M|/L = 0,0256)."""
        self.form.set(dict(tipo=TIPOS[0], f0=2.45, k=0.0256, l=3.75, L=10, tau=15, bits="11", npb=0, b=1.0))
        self.calcular()

    @protegido
    def calcular(self):
        v = self.form.get()
        if not 0 < v["k"] < 0.5 or v["f0"] <= 0 or v["l"] <= 0 or v["L"] <= 0 or v["tau"] <= 0:
            raise ValueError("Use 0 < acoplamento < 0,5 e valores positivos.")
        bits = "".join(ch for ch in v["bits"] if ch in "01")
        if not bits:
            raise ValueError("Código deve conter 0 e 1.")
        f0, l, tau = v["f0"] * 1e9, v["l"] * 1e-3, v["tau"] * 1e-9
        kappa = {TIPOS[0]: -v["k"], TIPOS[1]: v["k"], TIPOS[2]: miw.fator_acoplamento(cm_sobre_cc=v["k"])}[v["tipo"]]
        fa, fb, fbw = miw.bordas(f0, kappa)
        proj = miw.projeto_etiqueta(f0, kappa, l, tau, n_bits=1)
        npb = v["npb"] if v["npb"] > 0 else proj["n_ressonadores"]
        lr = v["L"] * 1e-9
        cr = 1 / ((2 * np.pi * f0) ** 2 * lr)
        m = kappa * lr
        bl = np.linspace(1e-3, np.pi - 1e-3, 400)
        fb_l = miw.frequencia(bl, f0, kappa)
        vg = miw.velocidade_grupo(bl, f0, kappa, l)
        zb = abs(abcd.impedancia_bloch(miw.abcd_miw(np.array([f0]), lr, cr, m))[0][0])
        f = np.linspace(1e6, 2 * f0, 60000)
        s11 = miw.s11_etiqueta(f, lr, cr, m, npb, bits, v["b"], zb)
        t, env = miw.resposta_pulso(f, s11, f0, tau)
        t_bit = npb * proj["atraso_por_ressonador"]
        t_max = (len(bits) + 1.5) * t_bit + 2 * tau
        sel = t < t_max
        ax = self.grafico.eixos(1, 3)
        ax[0].plot(bl / np.pi, fb_l / 1e9)
        ax[0].set_xlabel("βl/π")
        ax[0].set_ylabel("f [GHz]")
        ax[0].set_title(f"Dispersão (eq. {'6.4' if v['tipo'] == TIPOS[2] else '6.1'})", fontsize=9)
        ax[1].plot(fb_l / 1e9, C0 / np.abs(vg))
        ax[1].set_yscale("log")
        ax[1].set_xlabel("f [GHz]")
        ax[1].set_ylabel("c/|vg|")
        ax[1].set_title("Fator de onda lenta (eq. 6.8)", fontsize=9)
        ax[2].plot(t[sel] * 1e9, env[sel], label="eco (envelope)")
        for k in range(1, len(bits) + 1):
            ax[2].axvline(k * t_bit * 1e9, color="gray", ls=":", lw=0.8)
        ax[2].set_xlabel("t [ns]")
        ax[2].set_title(f"Resposta da etiqueta '{bits}' a um pulso de {v['tau']:g} ns", fontsize=9)
        for a in ax:
            a.grid(True, alpha=0.3)
        self.grafico.desenhar()
        self.saida.escrever("\n".join([
            v["tipo"],
            f"Banda: {fa / 1e9:.4f} a {fb / 1e9:.4f} GHz (FBW = {100 * fbw:.2f}%, eqs. 6.3/6.6)",
            f"vg em f0 = {abs(l * kappa * 2 * np.pi * f0):.3g} m/s → c/vg = {proj['fator_lento']:.0f}",
            f"Atraso de ida e volta por ressonador: {proj['atraso_por_ressonador'] * 1e9:.3f} ns (eq. 6.9)", "",
            f"Etiqueta: {npb} ressonadores por bit (1 bit = Δτ/τ, eq. 6.10), {len(bits)} bits → "
            f"{npb * len(bits)} ressonadores, {npb * len(bits) * v['l']:.1f} mm",
            f"Ecos previstos em múltiplos de {t_bit * 1e9:.2f} ns (linhas pontilhadas).",
            f"Portas em ZB(f0) = {zb:.2f} Ω (modelo de circuito)", "",
            "O atraso cresce perto das bordas da banda: pulsos",
            "curtos (banda larga) chegam mais tarde e com ecos",
            "espúrios. Compromisso bits × banda × comprimento",
            "× perdas (livro, p. 429). Perdas não incluídas."]))
