"""Teoria de modos acoplados (CMT) para linhas EBG não uniformes.

Referência: F. Martín (2015), seção 2.4 (eqs. 2.47-2.105).

Unidades internas SI (m, Hz, 1/m). Convenções:
    da+/dz = -j*beta*a+ + K*a-
    da-/dz =  K*a+ + j*beta*a-
    K(z)   = -(1/2Z0) dZ0/dz                       (eq. 2.50)
    S11 = a-(0)/a+(0), S21 = a+(L)/a+(0) com a-(L) = 0   (eq. 2.56)
"""
from dataclasses import dataclass, field

import numpy as np

from . import abcd
from .microstrip import C0, Substrato, analise, sintese
from .microstrip import parametros_linha as _parametros_microstrip

TIPOS = ("Senoidal", "Quadrada")


@dataclass
class Perfil:
    tipo: str
    z: np.ndarray            # m, posições (nós)
    z0: np.ndarray           # ohm, impedância característica nos nós
    periodo: float           # m
    n_periodos: int
    harmonicos: dict = field(default_factory=dict)  # n -> |Kn| [1/m]

    @property
    def comprimento(self):
        return self.periodo * self.n_periodos


def perfil(tipo, k_amp, periodo, n_periodos, z_ref=50.0, pts_por_periodo=64):
    """Perfil Z0(z) para coeficiente de acoplamento senoidal ou quadrado.

    tipo "Senoidal": K(z) = -2|K1| cos(2*pi*z/l)  -> eq. 2.102 (k_amp = |K1|)
    tipo "Quadrada": K(z) = -/+K em quartos de período -> eq. 2.104 (k_amp = K)
    """
    n_pts = pts_por_periodo * n_periodos + 1
    z = np.linspace(0, periodo * n_periodos, n_pts)
    zl = np.mod(z, periodo)
    zl[-1] = periodo if n_periodos else 0  # último nó fecha o período
    if tipo == "Senoidal":
        z0 = z_ref * np.exp(2 * k_amp * periodo / np.pi * np.sin(2 * np.pi * zl / periodo))
        harm = {1: k_amp}
    elif tipo == "Quadrada":
        q = periodo / 4
        z0 = np.where(zl <= q, z_ref * np.exp(2 * k_amp * zl),
                      np.where(zl <= 3 * q, z_ref * np.exp(-2 * k_amp * (zl - 2 * q)),
                               z_ref * np.exp(2 * k_amp * (zl - periodo))))
        harm = {n: 2 * k_amp / (n * np.pi) for n in (1, 3, 5, 7)}  # eq. 2.105
    else:
        raise ValueError(f"Tipo de perfil desconhecido: {tipo}")
    return Perfil(tipo, z, z0, periodo, n_periodos, harm)


def coeficiente_acoplamento(p: Perfil):
    """K(z) nos pontos médios dos segmentos (eq. 2.50, por diferenças finitas)."""
    return -0.5 * np.diff(np.log(p.z0)) / np.diff(p.z)


def larguras(p: Perfil, sub: Substrato):
    """W(z) [mm] e eps_eff(z) do perfil."""
    w = sintese(p.z0, sub)
    _, ee = analise(w, sub)
    return w, ee


def eeff_medio(z, eeff):
    """sqrt(eeff_medio) = (1/L) * integral(sqrt(eeff(z)) dz)  (eq. 2.57)."""
    return (np.trapezoid(np.sqrt(eeff), z) / (z[-1] - z[0])) ** 2


def parametros_linha(f, p: Perfil, sub: Substrato, perdas=False, dispersivo=False):
    """(Z0, εeff, α) do perfil por frequência, forma (Nf, Nnós), para os solvers.

    A largura vem da síntese quase-estática (como no layout); com dispersão, Z0 e
    εeff variam com f. Sem as opções, equivale a (p.z0, εeff, 0) repetidos em f.
    """
    w, _ = larguras(p, sub)
    z0, ee, alfa = _parametros_microstrip(w, sub, f, perdas, dispersivo)
    if not dispersivo:
        z0 = np.broadcast_to(p.z0, ee.shape).copy()   # exatamente o perfil projetado
    return z0, ee, alfa


def _por_frequencia(f, *arrays):
    """Converte arrays (Nnós,) ou (Nf, Nnós) para (Nf, Nnós)."""
    return [np.broadcast_to(np.asarray(a, dtype=float), (f.size, np.shape(a)[-1])) for a in arrays]


def resolver_cmt(f, z, z0, eeff, alfa=None):
    """Solução numérica das equações acopladas (2.48) -> (S11, S21).

    Discretiza em segmentos com γ = α + jβ e K constantes; a matriz de transição
    de cada segmento é exp(M dz) = cosh(g dz) I + sinh(g dz)/g M, com
    M = [[-γ, K], [K, γ]] e g^2 = K^2 + γ^2 (sem perdas: g^2 = K^2 - β^2).
    z0, eeff e alfa [Np/m] podem ser (Nnós,) ou (Nf, Nnós) (dispersão e perdas).
    """
    f = np.asarray(f, dtype=float)
    if alfa is None:
        alfa = np.zeros(np.shape(eeff)[-1])
    z0, eeff, alfa = _por_frequencia(f, z0, eeff, alfa)
    dz = np.diff(z)
    k = -0.5 * np.diff(np.log(z0), axis=1) / dz
    ee_m = 0.5 * (eeff[:, 1:] + eeff[:, :-1])
    al_m = 0.5 * (alfa[:, 1:] + alfa[:, :-1])
    t = np.broadcast_to(np.eye(2, dtype=complex), (f.size, 2, 2)).copy()
    p = np.empty((f.size, 2, 2), dtype=complex)
    for i, dzi in enumerate(dz):
        ki = k[:, i]
        gama = al_m[:, i] + 2j * np.pi * f * np.sqrt(ee_m[:, i]) / C0
        g = np.sqrt(ki**2 + gama**2 + 0j)
        ch = np.cosh(g * dzi)
        sh_g = np.where(np.abs(g) * dzi < 1e-12, dzi, np.sinh(g * dzi) / np.where(g == 0, 1, g))
        p[:, 0, 0] = ch - gama * sh_g
        p[:, 0, 1] = ki * sh_g
        p[:, 1, 0] = ki * sh_g
        p[:, 1, 1] = ch + gama * sh_g
        t = p @ t
    s11 = -t[:, 1, 0] / t[:, 1, 1]
    s21 = t[:, 0, 0] + t[:, 0, 1] * s11
    return s11, s21


def resolver_abcd(f, z, z0, eeff, z_ref=50.0, alfa=None):
    """Verificação independente: cascata de trechos uniformes curtos (θ = (β - jα)·dz)."""
    f = np.asarray(f, dtype=float)
    if alfa is None:
        alfa = np.zeros(np.shape(eeff)[-1])
    z0, eeff, alfa = _por_frequencia(f, z0, eeff, alfa)
    dz = np.diff(z)
    zm = np.sqrt(z0[:, 1:] * z0[:, :-1])
    ee_m = 0.5 * (eeff[:, 1:] + eeff[:, :-1])
    al_m = 0.5 * (alfa[:, 1:] + alfa[:, :-1])
    t = np.broadcast_to(np.eye(2, dtype=complex), (f.size, 2, 2)).copy()
    for i, dzi in enumerate(dz):
        theta = (2 * np.pi * f * np.sqrt(ee_m[:, i]) / C0 - 1j * al_m[:, i]) * dzi
        t = t @ abcd.linha(zm[:, i], theta)
    s = abcd.abcd_para_s(t, z_ref)
    return s[:, 0, 0], s[:, 1, 0]


# ---------------------------------------------------------------------------
# Soluções analíticas aproximadas (vizinhança da n-ésima banda rejeitada)
# ---------------------------------------------------------------------------
@dataclass
class BandaAnalitica:
    n: int
    k: float          # |Kn| [1/m]
    f_max: float      # Hz (eq. 2.84)
    s11_max: float    # linear (eq. 2.85)
    s21_min: float    # linear (eq. 2.86)
    banda: float      # Hz (eq. 2.89)


def bandas_analiticas(p: Perfil, eeff_med):
    out = []
    comp = p.comprimento
    for n, k in sorted(p.harmonicos.items()):
        f_max = n * C0 / (2 * p.periodo * np.sqrt(eeff_med))
        bw = C0 * k / (np.pi * np.sqrt(eeff_med)) * np.sqrt(1 + (np.pi / (k * comp)) ** 2)
        out.append(BandaAnalitica(n, k, f_max, np.tanh(k * comp), 1 / np.cosh(k * comp), bw))
    return out


def s_analitico(f, n, k, periodo, comprimento, eeff_med):
    """S11 e S21 aproximados perto da banda n (eqs. 2.76a,b)."""
    beta = 2 * np.pi * np.asarray(f) * np.sqrt(eeff_med) / C0
    db = beta - n * np.pi / periodo
    g = np.sqrt(k**2 - db**2 + 0j)
    g = np.where(np.abs(g) < 1e-12, 1e-12, g)
    den = -db * np.sinh(g * comprimento) + 1j * g * np.cosh(g * comprimento)
    s11 = k * np.sinh(g * comprimento) / den
    s21 = 1j * g * np.exp(-1j * n * np.pi * comprimento / periodo) / den
    return s11, s21


def frequencia_bragg(periodo, eeff, n=1):
    """Condição de Bragg f = n*c/(2*l*sqrt(eeff)) (eqs. 2.9 / 2.84)."""
    return n * C0 / (2 * periodo * np.sqrt(eeff))
