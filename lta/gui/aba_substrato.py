"""Aba Substrato: parâmetros do dielétrico e calculadora microstrip."""
import numpy as np

from ..core.microstrip import Substrato, analise, sintese
from .comum import AbaBase, AreaTexto, Formulario, protegido

# Valores nominais para referência rápida. Conferir sempre o datasheet do lote
# e a frequência de projeto antes de usar (Dk de projeto pode diferir do nominal).
PREDEFINIDOS = {
    "Rogers RO3010 (livro)": Substrato(10.2, 1.27, 35, 0.0022, "Rogers RO3010"),
    "Arlon er=2.4 h=0.675 (livro, Fig. 2.33)": Substrato(2.4, 0.675, 35, 0.0015, "Arlon er2.4"),
    "Rogers RO4003C 0.813 mm": Substrato(3.55, 0.813, 35, 0.0027, "Rogers RO4003C"),
    "FR-4 1.6 mm": Substrato(4.4, 1.6, 35, 0.02, "FR-4"),
}


class AbaSubstrato(AbaBase):
    TITULO = "Substrato"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        s = estado.substrato
        self.form = Formulario(self.painel, [
            ("predef", "Predefinido", tuple(PREDEFINIDOS), ""),
            ("nome", "Nome", s.nome, ""),
            ("er", "εr", s.er, ""),
            ("h_mm", "Espessura h", s.h_mm, "mm"),
            ("t_um", "Cobre t", s.t_um, "µm"),
            ("tand", "tan δ", s.tand, ""),
        ], titulo="Substrato (usado por todas as abas)")
        self.form.pack(fill="x")
        self.form.ao_mudar("predef", self._carregar_predef)
        self.botoes(("Aplicar", self.aplicar))

        self.calc = Formulario(self.painel, [
            ("z0", "Z0 desejada", 50.0, "Ω"),
            ("w", "Largura W", 1.2, "mm"),
        ], titulo="Calculadora microstrip (Hammerstad-Jensen)")
        self.calc.pack(fill="x", pady=(10, 0))
        self.botoes(("Z0 → W", self.z_para_w), ("W → Z0", self.w_para_z))
        self.saida = AreaTexto(self.painel, 8)
        self.saida.pack(fill="both", expand=True)
        estado.ouvir(lambda ev: ev == "substrato" and self._sincronizar())
        self.aplicar()

    def _carregar_predef(self):
        s = PREDEFINIDOS[self.form.get()["predef"]]
        self.form.set(dict(nome=s.nome, er=s.er, h_mm=s.h_mm, t_um=s.t_um, tand=s.tand))

    def _sincronizar(self):
        s = self.estado.substrato
        self.form.set(dict(nome=s.nome, er=s.er, h_mm=s.h_mm, t_um=s.t_um, tand=s.tand))
        self._plotar()

    @protegido
    def aplicar(self):
        v = self.form.get()
        if v["er"] < 1 or v["h_mm"] <= 0:
            raise ValueError("εr deve ser ≥ 1 e h > 0.")
        self.estado.definir_substrato(Substrato(v["er"], v["h_mm"], v["t_um"], v["tand"], v["nome"]))
        self.estado.msg(f"Substrato aplicado: {v['nome']} (εr={v['er']:g}, h={v['h_mm']:g} mm)")

    @protegido
    def z_para_w(self):
        z0 = self.calc.get()["z0"]
        w = sintese(z0, self.estado.substrato)
        _, ee = analise(w, self.estado.substrato)
        self.calc.set({"w": f"{w:.4f}"})
        self._relatorio(w, z0, ee)

    @protegido
    def w_para_z(self):
        w = self.calc.get()["w"]
        z0, ee = analise(w, self.estado.substrato)
        self.calc.set({"z0": f"{float(z0):.3f}"})
        self._relatorio(w, float(z0), float(ee))

    def _relatorio(self, w, z0, ee):
        lg = 299.792458 / np.sqrt(ee)  # mm a 1 GHz
        self.saida.escrever(
            f"W = {w:.4f} mm   (W/h = {w / self.estado.substrato.h_mm:.3f})\n"
            f"Z0 = {z0:.3f} Ω\nεeff = {ee:.4f}\n"
            f"λg a 1 GHz = {lg:.2f} mm\n\n"
            "Modelo quase-estático, t = 0, sem dispersão.\n"
            "Fonte: Hammerstad & Jensen, IEEE MTT-S 1980.")

    def _plotar(self):
        s = self.estado.substrato
        u = np.logspace(-1.3, 1.3, 200)
        z, ee = analise(u * s.h_mm, s)
        ax1, ax2 = self.grafico.eixos(2, 1, sharex=True)
        ax1.semilogx(u, z)
        ax1.set_ylabel("Z0 [Ω]")
        ax1.grid(True, which="both", alpha=0.3)
        ax1.set_title(f"Microstrip em {s.nome}: εr = {s.er:g}, h = {s.h_mm:g} mm", fontsize=10)
        ax2.semilogx(u, ee, color="C1")
        ax2.set_ylabel("εeff")
        ax2.set_xlabel("W/h")
        ax2.grid(True, which="both", alpha=0.3)
        self.grafico.desenhar()
