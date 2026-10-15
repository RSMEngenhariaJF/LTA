"""Aba Floquet: harmônicos espaciais e radiação leaky-wave (seção 2.2)."""
import numpy as np

from ..core import abcd
from ..core import floquet as fl
from ..core import modos_acoplados as cm
from ..core.microstrip import C0, analise, parametros_linha, sintese
from .aba_substrato import PREDEFINIDOS
from .comum import AbaBase, AreaTexto, Formulario, OpcoesModelo, protegido

ESTRUTURAS = ("EBG senoidal (largura)", "EBG quadrada (largura)", "Linha com carga capacitiva",
              "Linha uniforme (Fig. 2.1)")


class AbaFloquet(AbaBase):
    TITULO = "Floquet / harmônicos espaciais"

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("tipo", "Célula", ESTRUTURAS, ""),
            ("periodo", "Período l", 10.0, "mm"),
            ("k", "|K| (EBG)", 0.64, "1/cm"),
            ("zref", "Z referência (EBG)", 50.0, "Ω"),
            ("zo", "Zo da linha (carga)", 100.0, "Ω"),
            ("cls", "Cls (carga)", 0.5, "pF"),
            ("nmax", "Harmônicos ±n", 2, ""),
            ("fmax", "f máx", 25.0, "GHz"),
            ("npts", "Pontos", 2000, ""),
        ], titulo="Célula unitária periódica")
        self.form.pack(fill="x")
        self.modelo_fis = OpcoesModelo(self.painel)
        self.modelo_fis.pack(fill="x", pady=(8, 0))
        self.botoes(("Calcular", self.calcular), ("Ex. EBG Fig. 2.11", self.ex_ebg),
                    ("Ex. Fig. 2.1 (linha)", self.ex_linha), ("Ex. carga capacitiva", self.ex_carga))
        self.saida = AreaTexto(self.painel, 18)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Harmônicos espaciais βn = β + 2πn/l (eq. 2.6).\n"
            "Onda rápida |βn| < k0: o harmônico radia (leaky wave)\n"
            "com cos θ = βn/k0 (eq. 2.11), θ medido a partir do\n"
            "eixo da linha (θ = 90°: broadside).\n\n"
            "βn < 0: radiação para trás (backward, espectral);\n"
            "βn > 0: para a frente (forward, não espectral) (eq. 2.13).\n\n"
            "A região sombreada do diagrama é a de onda rápida\n"
            "(entre as retas ±k0·l, Fig. 2.2).")

    def _ro3010(self):
        self.estado.definir_substrato(PREDEFINIDOS["Rogers RO3010 (livro)"])

    def ex_ebg(self):
        self._ro3010()
        self.form.set(dict(tipo=ESTRUTURAS[0], periodo=10, k=0.64, zref=50, nmax=2, fmax=25, npts=2500))
        self.calcular()

    def ex_linha(self):
        self._ro3010()
        self.form.set(dict(tipo=ESTRUTURAS[3], periodo=10, zref=50, nmax=2, fmax=25, npts=1000))
        self.calcular()

    def ex_carga(self):
        self._ro3010()
        self.form.set(dict(tipo=ESTRUTURAS[2], periodo=5, zo=100, cls=0.3, nmax=2, fmax=25, npts=2500))
        self.calcular()

    def _celula(self, f, v):
        """ABCD da célula por frequência e εeff de referência (linha não perturbada)."""
        sub = self.estado.substrato
        perdas, dispersivo = self.modelo_fis.get()
        l = v["periodo"] * 1e-3
        if v["tipo"] in ESTRUTURAS[:2]:
            tipo = "Senoidal" if v["tipo"] == ESTRUTURAS[0] else "Quadrada"
            p = cm.perfil(tipo, v["k"] * 100, l, 1, v["zref"], 128)
            z0, ee, al = cm.parametros_linha(f, p, sub, perdas, dispersivo)
            return cm.abcd_perfil(f, p.z, z0, ee, al), cm.eeff_medio(p.z, ee[0])
        z_linha = v["zo"] if v["tipo"] == ESTRUTURAS[2] else v["zref"]
        w = sintese(z_linha, sub)
        z0, ee, al = (x[:, 0] for x in parametros_linha(w, sub, f, perdas, dispersivo))
        theta = (2 * np.pi * f * np.sqrt(ee) / C0 - 1j * al) * l
        if v["tipo"] == ESTRUTURAS[3]:
            return abcd.linha(z0, theta), float(analise(w, sub)[1])
        if v["cls"] <= 0:
            raise ValueError("Cls deve ser positiva.")
        meia = abcd.linha(z0, theta / 2)
        return abcd.cascata(meia, abcd.shunt(2j * np.pi * f * v["cls"] * 1e-12), meia), float(analise(w, sub)[1])

    @protegido
    def calcular(self):
        v = self.form.get()
        if v["periodo"] <= 0 or v["fmax"] <= 0 or v["npts"] < 50:
            raise ValueError("Período e f máx devem ser positivos; use ao menos 50 pontos.")
        if not 0 <= v["nmax"] <= 5:
            raise ValueError("Número de harmônicos entre 0 e 5.")
        l = v["periodo"] * 1e-3
        f = np.linspace(v["fmax"] * 1e9 / v["npts"], v["fmax"] * 1e9, v["npts"])
        m, ee_ref = self._celula(f, v)
        al, bl = abcd.dispersao(m)
        ext = fl.desdobrar(bl)
        harm = fl.harmonicos(f, ext, l, v["nmax"])
        k0l = 2 * np.pi * f / C0 * l
        fg = f / 1e9

        ax = self.grafico.eixos(2, 2)
        lim = 2 * v["nmax"] + 1.5
        ax[0].fill_betweenx(fg, -k0l / np.pi, k0l / np.pi, color="#cfe3f5", label="Onda rápida (radia)")
        for h in harm:
            cor = "C0" if h.n == 0 else f"C{abs(h.n) % 9 + 1}"
            larg = 1.8 if h.n == 0 else 1.0
            ax[0].plot(h.beta_l / np.pi, fg, color=cor, lw=larg, label=f"n = {h.n:+d}" if h.n else "n = 0 (fund.)")
            ax[0].plot(-h.beta_l / np.pi, fg, color=cor, lw=0.7, ls=":")
            rad = np.where(h.rapido, h.beta_l / np.pi, np.nan)
            ax[0].plot(rad, fg, color="k", lw=3, alpha=0.35)
        for x in (-1, 1):
            ax[0].axvline(x, color="gray", ls=":", lw=0.8)
        ax[0].set_xlim(-lim, lim)
        ax[0].set_ylim(0, fg[-1])
        ax[0].set_xlabel("βn·l/π")
        ax[0].set_ylabel("f [GHz]")
        ax[0].set_title("Diagrama de Floquet, zona estendida (Fig. 2.2)", fontsize=9)
        ax[0].legend(fontsize=6, loc="upper left", ncol=2)

        ax[1].plot(fg, al * 20 / np.log(10), label="αl [dB/célula]")
        ax[1].set_ylabel("Atenuação [dB/célula]")
        ax2 = ax[1].twinx()
        ax2.plot(fg, ext / np.pi, color="C1", lw=1, label="βl/π estendido")
        ax2.set_ylabel("βl/π (fundamental)", color="C1")
        ax[1].set_title("Bandas rejeitadas (eq. 2.21)", fontsize=9)

        for h in harm:
            if h.rapido.any():
                ax[2].plot(fg, h.theta_graus, label=f"n = {h.n:+d}")
        ax[2].axhline(90, color="gray", ls=":", lw=0.8)
        ax[2].set_ylim(0, 180)
        ax[2].set_ylabel("θ [°] (cos θ = βn/k0)")
        ax[2].set_title("Ângulo do feixe leaky-wave (eq. 2.11)", fontsize=9)

        vg = fl.velocidade_grupo(f, ext, l)
        ax[3].plot(fg, np.where(al < 1e-3, vg, np.nan), label="vg/c (todos os n, eq. 2.7)")
        for h in harm:
            if h.n in (0, -1):
                ax[3].plot(fg, np.clip(h.vp_c, -3, 3), ls="--", label=f"vp/c, n = {h.n:+d} (eq. 2.8)")
        ax[3].axhline(1, color="gray", ls=":", lw=0.8)
        ax[3].axhline(-1, color="gray", ls=":", lw=0.8)
        ax[3].set_ylim(-3, 3)
        ax[3].set_title("Velocidades de fase e de grupo", fontsize=9)
        for a in ax[1:]:
            a.set_xlabel("f [GHz]")
            a.legend(fontsize=7)
        for a in ax:
            a.grid(True, alpha=0.3)
        self.grafico.desenhar()
        self.saida.escrever(self._texto(v, f, al, ext, harm, ee_ref))

    def _texto(self, v, f, al, ext, harm, ee_ref):
        l = v["periodo"] * 1e-3
        linhas = [f"{v['tipo']}: l = {v['periodo']:g} mm, εeff (sem perturbação) = {ee_ref:.3f}",
                  self.modelo_fis.descrever(self.estado.substrato), "",
                  "Bandas rejeitadas (αl > 0,05 dB/célula):"]
        rej = al * 20 / np.log(10) > 0.05
        i, achou = 0, False
        while i < len(f):
            if rej[i]:
                j = i
                while j + 1 < len(f) and rej[j + 1]:
                    j += 1
                k = i + int(np.argmax(al[i:j + 1]))
                linhas.append(f"  {f[i] / 1e9:.3f} a {f[j] / 1e9:.3f} GHz, βl ≈ {ext[k] / np.pi:.0f}π, "
                              f"máx {al[k] * 20 / np.log(10):.2f} dB/célula")
                achou = True
                i = j + 1
            else:
                i += 1
        if not achou:
            linhas.append("  nenhuma na faixa")
        linhas += ["", f"Bragg (eq. 2.9), 1ª banda: f = c/(2l√εeff) = "
                       f"{cm.frequencia_bragg(l, ee_ref) / 1e9:.3f} GHz",
                   f"Linha pouco perturbada, n = -1: radiação a partir de "
                   f"{fl.inicio_radiacao(l, ee_ref) / 1e9:.3f} GHz (backfire),",
                   f"  broadside em {fl.frequencia_broadside(l, ee_ref) / 1e9:.3f} GHz", "",
                   "Radiação calculada (|βn| < k0):"]
        algum = False
        for h in harm:
            for fa, fb, ta, tb in fl.faixas_radiacao(f, h):
                sentido = "para trás → para a frente" if ta > 90 > tb else ("para trás" if ta > 90 else "para a frente")
                linhas.append(f"  n = {h.n:+d}: {fa / 1e9:.3f} a {fb / 1e9:.3f} GHz, θ de {ta:.0f}° a {tb:.0f}° "
                              f"({sentido})")
                algum = True
        if not algum:
            linhas.append("  nenhum harmônico rápido na faixa")
        linhas += ["", "Perto de βn = 0 (broadside) há banda rejeitada",
                   "(n-ésima banda de Bragg), que interrompe o feixe.",
                   "A radiação não entra no cálculo de α (só perdas",
                   "ôhmicas e dielétricas, se ligadas)."]
        return "\n".join(linhas)
