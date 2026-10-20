"""Aba Tolerâncias, potência e orçamento de erros (Garg, Bahl e Bianchi, seções 2.4.2-2.4.6 e 4.5.3)."""
import numpy as np

from ..core import linhas_acopladas as la
from ..core import robustez as rb
from ..core.descontinuidades import texto_limites
from ..core.microstrip import sintese
from .comum import AbaBase, AreaTexto, Formulario, protegido

LINHAS = ("Microstrip", "Linhas acopladas (W, S)")
DEFINIR = ("Pela impedância Z0", "Pela largura W")


class AbaRobustez(AbaBase):
    TITULO = "Tolerâncias, potência e orçamento de erros"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("linha", "Estrutura", LINHAS, ""),
            ("definir", "Microstrip definida", DEFINIR, ""),
            ("z0", "Z0", 50.0, "Ω"),
            ("w", "Largura W", 1.8, "mm"),
            ("s", "Separação S (acopladas)", 0.3, "mm"),
            ("dw", "Tolerância ±ΔW (e ±ΔS)", 0.025, "mm"),
            ("dh", "Tolerância ±Δh", 0.025, "mm"),
            ("der", "Tolerância ±Δεr", 0.05, ""),
            ("dz", "Especificação: |ΔZ0/Z0| máx.", 2.0, "%"),
            ("f", "Frequência de análise", 10.0, "GHz"),
            ("material", "Material (K e rigidez)", tuple(rb.MATERIAIS), ""),
            ("k", "Condutividade térmica K (0 = do material)", 0.0, "W/m/°C"),
            ("tmax", "Temperatura máxima", 100.0, "°C"),
            ("tamb", "Temperatura ambiente", 25.0, "°C"),
            ("vswr_con", "VSWR do conector (orçamento)", 1.10, ""),
            ("w2", "Degrau para W2/W1 (0 = sem)", 0.5, ""),
        ], titulo="Projeto robusto")
        self.form.pack(fill="x")
        self.botoes(("Calcular", self.calcular))
        self.saida = AreaTexto(self.painel, 24)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "• Tolerâncias: sensibilidades S = (B/A)·∂A/∂B e pior caso\n"
            "  (eqs. 2.97-2.105; acopladas, 4.120-4.127), mais Monte Carlo.\n"
            "• Potência média: aquecimento com o modelo de placas paralelas\n"
            "  (eqs. 2.84-2.95); de pico: ruptura do dielétrico (eq. 2.96).\n"
            "  Conectores (SMA) e bordas da trilha costumam limitar antes.\n"
            "• Orçamento de erros como na Tabela 2.5 (espessura, dispersão,\n"
            "  tolerâncias, degrau, conector).")

    @protegido
    def calcular(self):
        v = self.form.get()
        sub = self.estado.substrato
        f = v["f"] * 1e9
        w = float(sintese(v["z0"], sub)) if v["definir"] == DEFINIR[0] else v["w"]
        k_mat, e_rup = rb.MATERIAIS[v["material"]]
        k = v["k"] if v["k"] > 0 else k_mat
        linhas = [f"Substrato: {sub.nome} (εr = {sub.er:g}, h = {sub.h_mm:g} mm, t = {sub.t_um:g} µm)", ""]
        ax = self.grafico.eixos(2, 2)
        if v["linha"] == LINHAS[1]:
            r = rb.pior_caso_acopladas(w, v["s"], sub, v["dw"], v["dw"], v["dh"], v["der"])
            linhas += [f"Linhas acopladas W = {w:.4f} mm, S = {v['s']:.4f} mm ({la.MODELOS[0]}):",
                       f"  Z0e = {r['ze']:.2f} Ω ± {100 * r['dze']:.2f}%, Z0o = {r['zo']:.2f} Ω ± {100 * r['dzo']:.2f}%",
                       f"  Z0 = √(Z0e·Z0o) ± {100 * r['dz0']:.2f}% (eq. 4.125)",
                       f"  C = {r['c_db']:.2f} dB ± {r['dc_db']:.3f} dB (pior caso)", ""]
        sens = rb.sensibilidades(w, sub)
        pc = rb.pior_caso(w, sub, v["dw"], v["dh"], v["der"])
        linhas += [f"Microstrip W = {w:.4f} mm:",
                   "  Sensibilidades (S de Z0, S de εeff): " + "; ".join(f"{k_}: {a:+.3f}, {b:+.3f}"
                                                                       for k_, (a, b) in sens.items()),
                   f"  Pior caso: ΔZ0/Z0 = ±{100 * pc.dz_rel:.2f}% (VSWR {pc.vswr:.3f}), "
                   f"Δεeff/εeff = ±{100 * pc.de_rel:.2f}%"]
        dw_req = rb.precisao_w(w, sub, v["dz"] / 100, 1.0, v["dh"], v["der"])
        linhas.append(f"  Para |ΔZ0/Z0| ≤ {v['dz']:g}%: " + (f"ΔW ≤ ±{dw_req * 1e3:.1f} µm (eqs. 2.100-2.101)"
                                                              if dw_req else "impossível só com ΔW: h e εr já excedem"))
        z, e = rb.monte_carlo(w, sub, v["dw"], v["dh"], v["der"], n=1500)
        linhas.append(f"  Monte Carlo (1500, uniformes): Z0 = {z.mean():.2f} ± {z.std():.2f} Ω "
                      f"[{z.min():.2f}; {z.max():.2f}], εeff = {e.mean():.4f} ± {e.std():.4f}")
        ax[0].hist(z, bins=40, color="C0", alpha=0.8)
        ax[0].set_xlabel("Z0 [Ω]")
        ax[0].set_title("Monte Carlo de Z0", fontsize=9)
        ax[1].hist(e, bins=40, color="C1", alpha=0.8)
        ax[1].set_xlabel("εeff")
        ax[1].set_title("Monte Carlo de εeff", fontsize=9)
        fs = np.linspace(0.1e9, max(2 * f, 2e9), 200)
        pt = rb.potencia(w, sub, fs, k, v["tmax"], v["tamb"], e_rup)
        ax[2].semilogy(fs / 1e9, pt.p_media)
        ax[2].set_xlabel("f [GHz]")
        ax[2].set_ylabel("P média máx. [W]")
        ax[2].set_title(f"Potência média (K = {k:g} W/m/°C)", fontsize=9)
        ax[3].plot(fs / 1e9, pt.dt_por_w)
        ax[3].set_xlabel("f [GHz]")
        ax[3].set_ylabel("ΔT por watt [°C/W]")
        ax[3].set_title("Aquecimento da trilha", fontsize=9)
        for a in ax:
            a.grid(True, alpha=0.3)
        self.grafico.desenhar()
        p0 = rb.potencia(w, sub, f, k, v["tmax"], v["tamb"], e_rup)
        linhas += ["", f"Potência em {v['f']:g} GHz (Tmax = {v['tmax']:g} °C, Tamb = {v['tamb']:g} °C):",
                   f"  ΔT = {p0.dt_por_w[0]:.4g} °C/W → P média ≤ {p0.p_media[0]:.4g} W (largura térmica "
                   f"We = {p0.we:.3f} mm)"]
        if np.isfinite(p0.p_pico):
            linhas.append(f"  P de pico ≤ {p0.p_pico:.4g} W pela ruptura do dielétrico ({e_rup:g} kV/cm × h); "
                          "o ar nas bordas (30 kV/cm) e o conector limitam bem antes")
        linhas.append("  Resistências térmicas de solda/fixação e o conector não estão incluídas.")
        tab, total = rb.orcamento(w, sub, f, v["dw"], v["dh"], v["der"], v["vswr_con"] or None,
                                  w * v["w2"] if v["w2"] > 0 else None)
        linhas += ["", f"Orçamento de erros (Tabela 2.5) em {v['f']:g} GHz:"]
        for nome, vs, de in tab:
            linhas.append(f"  {nome:48s} VSWR {vs:.3f}" + (f"   Δεeff {de:+.4f}" if de is not None else ""))
        linhas.append(f"  Pior caso combinado (produto dos VSWR): {total:.3f}")
        linhas += ["", "Limites do substrato:"] + texto_limites(sub, w, f)
        self.saida.escrever("\n".join(linhas))
