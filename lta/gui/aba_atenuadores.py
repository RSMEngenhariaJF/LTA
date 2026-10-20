"""Aba Atenuadores resistivos T e Pi (Steer, seção 7.3.2)."""
import numpy as np

from ..core import passivos as pv
from .comum import AbaBase, AreaTexto, Formulario, protegido


class AbaAtenuadores(AbaBase):
    TITULO = "Atenuadores resistivos (pads T e Pi)"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("top", "Topologia", pv.TOPOLOGIAS, ""),
            ("z01", "Z01 (porta 1)", 50.0, "Ω"),
            ("z02", "Z02 (porta 2)", 50.0, "Ω"),
            ("k", "Atenuação", 10.0, "dB"),
            ("serie", "Série comercial", tuple(pv.SERIES), ""),
            ("pin", "Potência de entrada", 1.0, "W"),
        ], titulo="Pad resistivo")
        self.form.pack(fill="x")
        self.botoes(("Calcular", self.calcular), ("Pad 50 → 75 Ω mínimo", self.ex_75))
        self.saida = AreaTexto(self.painel, 22)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Steer, eqs. 7.2-7.14: pads casados nas duas portas (Γ = 0 em\n"
            "Z01 e em Z02). Com Z01 ≠ Z02 há uma atenuação mínima (eq. 7.14):\n"
            "50 → 75 Ω exige pelo menos 5,72 dB.\n\n"
            "Mostra os resistores ideais, as versões balanceadas, os valores\n"
            "E24/E96 mais próximos com o erro de atenuação e casamento\n"
            "resultante, e a potência dissipada em cada resistor.\n"
            "Resistores concentrados servem até alguns GHz; acima disso use\n"
            "filme resistivo distribuído (linha com perdas).")

    def ex_75(self):
        kmin = 10 * np.log10(pv.k_minimo(50, 75))
        self.form.set(dict(top=pv.TOPOLOGIAS[0], z01=50, z02=75, k=f"{kmin + 1e-6:.6f}"))
        self.calcular()

    @protegido
    def calcular(self):
        v = self.form.get()
        a = pv.projeto(v["top"], v["z01"], v["z02"], v["k"])
        rc, k_c, s11_c, s22_c = pv.efeito_arredondamento(a, v["serie"])
        p1, p2, p3, pout = a.potencias(v["pin"])
        kmin_db = 10 * np.log10(pv.k_minimo(v["z01"], v["z02"]))
        ks = np.linspace(max(kmin_db, 0.0) + 0.05, 40, 300)
        rs, ps = [], []
        for k in ks:
            ak = pv.projeto(v["top"], v["z01"], v["z02"], k)
            rs.append(ak.r)
            ps.append(ak.potencias(1.0)[:3])
        rs, ps = np.array(rs), np.array(ps)
        ax = self.grafico.eixos(1, 2)
        for i in range(3):
            ax[0].semilogy(ks, np.where(np.isfinite(rs[:, i]) & (rs[:, i] > 0), rs[:, i], np.nan), label=f"R{i + 1}")
            ax[0].plot([v["k"]], [a.r[i]] if np.isfinite(a.r[i]) and a.r[i] > 0 else [np.nan], "o", color=f"C{i}")
            ax[1].plot(ks, 100 * ps[:, i], label=f"R{i + 1}")
        ax[0].set_ylabel("R [Ω]")
        ax[0].set_title(f"Pad {v['top']}: {v['z01']:g} → {v['z02']:g} Ω", fontsize=9)
        ax[1].set_ylabel("Potência dissipada [% de Pentrada]")
        ax[1].set_title("Dissipação por resistor", fontsize=9)
        for x in ax:
            x.axvline(v["k"], color="gray", ls=":", lw=0.8)
            if kmin_db > 0:
                x.axvline(kmin_db, color="C3", ls="--", lw=0.8, label=f"Kmín = {kmin_db:.2f} dB")
            x.set_xlabel("Atenuação [dB]")
            x.grid(True, alpha=0.3, which="both")
            x.legend(fontsize=7)
        self.grafico.desenhar()
        linhas = a.texto + ["", f"Valores {v['serie']}: " + ", ".join(pv._fmt(x) for x in rc),
                            f"  atenuação {k_c:.3f} dB (erro {k_c - v['k']:+.3f} dB), |S11| = {s11_c:.1f} dB, "
                            f"|S22| = {s22_c:.1f} dB",
                            "", f"Com {v['pin']:g} W na entrada: R1 {p1:.4g} W, R2 {p2:.4g} W, R3 {p3:.4g} W; "
                                f"saída {pout:.4g} W",
                            "Escolha a potência nominal dos resistores com folga (≥ 2×) e verifique a tensão máxima."]
        self.saida.escrever("\n".join(linhas))
