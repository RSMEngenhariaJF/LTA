"""Amplificador com transistor: estabilidade, ganhos e casamento conjugado simultâneo.

Referência: M. Steer, Microwave and RF Design: A Systems Approach, SciTech, 2010, seções 11.5-11.7:
  - fator de Rollett k e |Δ| (eqs. 11.62 e seguintes), μ de Edwards-Sinsky, círculos de estabilidade;
  - ganhos GT (eq. 11.93), GTU,max (11.101), GMA (11.100), GMS (11.102) e U (11.103);
  - casamento conjugado simultâneo (ΓS, ΓL) quando k > 1 e |Δ| < 1.
Os casamentos de entrada e saída são projetados com casamento.py (Z de ΓS e ΓL).
"""
from dataclasses import dataclass

import numpy as np

# Fig. 11-8 (Steer): pHEMT FPD6836P70 (RFMD), VDS = 5 V, ID = 55 mA. f [GHz], |S11| ∠, |S21| ∠, |S12| ∠, |S22| ∠
PHEMT_STEER = np.array([
    [0.5, 0.976, -20.9, 11.395, 161.5, 0.011, 78.3, 0.635, -11.5],
    [1, 0.925, -41.3, 10.729, 145.1, 0.021, 67.8, 0.614, -22.2],
    [2, 0.796, -78.2, 8.842, 116.7, 0.034, 51.4, 0.553, -37.9],
    [3, 0.694, -106.8, 7.180, 94.5, 0.041, 40.4, 0.506, -48.9],
    [4, 0.614, -127.3, 6.002, 76.7, 0.044, 33.9, 0.475, -57.7],
    [5, 0.555, -147.0, 5.249, 60.3, 0.048, 28.4, 0.453, -66.4],
    [6, 0.511, -170.2, 4.729, 43.7, 0.052, 23.3, 0.438, -76.0],
    [7, 0.493, 163.9, 4.261, 26.8, 0.057, 14.0, 0.391, -87.6],
    [8, 0.486, 140.4, 3.784, 11.2, 0.057, 6.4, 0.340, -99.1],
    [9, 0.473, 122.5, 3.448, -2.4, 0.059, 5.2, 0.332, -109.6],
    [10, 0.488, 103.4, 3.339, -17.3, 0.073, 0.9, 0.355, -124.8],
    [11, 0.539, 79.8, 3.166, -35.0, 0.086, -10.1, 0.349, -145.6],
    [12, 0.626, 60.8, 2.877, -51.9, 0.095, -21.4, 0.307, -169.6],
    [13, 0.685, 47.6, 2.604, -68.2, 0.100, -32.5, 0.295, 165.3],
    [14, 0.724, 36.2, 2.392, -83.8, 0.106, -43.3, 0.312, 142.7],
    [15, 0.787, 20.9, 2.225, -99.7, 0.109, -55.1, 0.320, 125.4],
    [16, 0.818, 5.2, 2.067, -116.6, 0.112, -68.4, 0.340, 103.9],
    [17, 0.831, -9.6, 1.855, -134.4, 0.108, -83.5, 0.373, 76.1],
    [18, 0.852, -19.5, 1.603, -148.6, 0.103, -94.2, 0.406, 54.7],
    [19, 0.815, -20.5, 1.440, -159.3, 0.102, -103.0, 0.449, 43.1],
    [20, 0.780, -26.8, 1.382, -171.2, 0.106, -113.5, 0.460, 37.9],
    [21, 0.779, -46.8, 1.333, 171.2, 0.109, -130.7, 0.438, 31.4],
    [22, 0.786, -62.1, 1.195, 152.0, 0.110, -148.4, 0.417, 6.0],
    [23, 0.774, -70.1, 1.073, 137.2, 0.108, -162.4, 0.428, -16.5],
    [24, 0.744, -81.7, 1.025, 123.5, 0.112, -175.2, 0.433, -29.0],
    [25, 0.704, -90.9, 1.061, 107.3, 0.132, 170.0, 0.396, -46.5],
    [26, 0.677, -111.1, 1.065, 85.8, 0.148, 147.8, 0.298, -71.0],
])


def s_de_tabela(tab):
    """(f [Hz], S (Nf, 2, 2)) de uma tabela [f GHz, |S11|, ∠S11, |S21|, ∠S21, |S12|, ∠S12, |S22|, ∠S22]."""
    c = lambda m, a: m * np.exp(1j * np.radians(a))          # noqa: E731
    s = np.empty((len(tab), 2, 2), dtype=complex)
    s[:, 0, 0] = c(tab[:, 1], tab[:, 2])
    s[:, 1, 0] = c(tab[:, 3], tab[:, 4])
    s[:, 0, 1] = c(tab[:, 5], tab[:, 6])
    s[:, 1, 1] = c(tab[:, 7], tab[:, 8])
    return tab[:, 0] * 1e9, s


@dataclass
class Estabilidade:
    k: np.ndarray
    delta: np.ndarray       # |Δ|
    mu: np.ndarray          # μ (entrada, Edwards-Sinsky)
    mu_l: np.ndarray        # μ' (saída)
    incondicional: np.ndarray


def estabilidade(s):
    """k = (1 - |S11|² - |S22|² + |Δ|²)/(2|S12S21|), Δ = S11S22 - S12S21; μ = (1 - |S11|²)/(|S22 - ΔS11*| + |S12S21|)."""
    s11, s12, s21, s22 = s[:, 0, 0], s[:, 0, 1], s[:, 1, 0], s[:, 1, 1]
    d = s11 * s22 - s12 * s21
    k = (1 - abs(s11) ** 2 - abs(s22) ** 2 + abs(d) ** 2) / (2 * abs(s12 * s21))
    mu = (1 - abs(s11) ** 2) / (abs(s22 - d * np.conj(s11)) + abs(s12 * s21))
    mu_l = (1 - abs(s22) ** 2) / (abs(s11 - d * np.conj(s22)) + abs(s12 * s21))
    return Estabilidade(k, abs(d), mu, mu_l, (k > 1) & (abs(d) < 1))


def circulos(s):
    """Círculos de estabilidade ((centro_S, raio_S), (centro_L, raio_L)): Γ em que |Γ_OUT| ou |Γ_IN| = 1."""
    s11, s12, s21, s22 = s[:, 0, 0], s[:, 0, 1], s[:, 1, 0], s[:, 1, 1]
    d = s11 * s22 - s12 * s21
    cl = np.conj(s22 - d * np.conj(s11)) / (abs(s22) ** 2 - abs(d) ** 2)
    rl = abs(s12 * s21 / (abs(s22) ** 2 - abs(d) ** 2))
    cs = np.conj(s11 - d * np.conj(s22)) / (abs(s11) ** 2 - abs(d) ** 2)
    rs = abs(s12 * s21 / (abs(s11) ** 2 - abs(d) ** 2))
    return (cs, rs), (cl, rl)


def ganhos(s):
    """Dicionário de ganhos lineares: GTU,max, GMS, GMA (onde k ≥ 1; senão NaN) e U."""
    s11, s12, s21, s22 = s[:, 0, 0], s[:, 0, 1], s[:, 1, 0], s[:, 1, 1]
    est = estabilidade(s)
    gtu = abs(s21) ** 2 / ((1 - abs(s11) ** 2) * (1 - abs(s22) ** 2))
    gms = abs(s21 / s12)
    with np.errstate(invalid="ignore"):
        gma = np.where(est.k >= 1, gms * (est.k - np.sqrt(np.maximum(est.k**2 - 1, 0))), np.nan)
    r = s21 / s12
    u = 0.5 * abs(r - 1) ** 2 / (est.k * abs(r) - r.real)
    return dict(GTU_max=gtu, GMS=gms, GMA=gma, U=u, S21=abs(s21) ** 2)


def conjugado_simultaneo(s):
    """(ΓS, ΓL) do casamento conjugado simultâneo (válido onde k > 1 e |Δ| < 1; senão NaN)."""
    s11, s12, s21, s22 = s[:, 0, 0], s[:, 0, 1], s[:, 1, 0], s[:, 1, 1]
    d = s11 * s22 - s12 * s21
    b1 = 1 + abs(s11) ** 2 - abs(s22) ** 2 - abs(d) ** 2
    b2 = 1 + abs(s22) ** 2 - abs(s11) ** 2 - abs(d) ** 2
    c1 = s11 - d * np.conj(s22)
    c2 = s22 - d * np.conj(s11)
    with np.errstate(invalid="ignore"):
        gs = (b1 - np.sign(b1) * np.sqrt(b1**2 - 4 * abs(c1) ** 2 + 0j)) / (2 * c1)
        gl = (b2 - np.sign(b2) * np.sqrt(b2**2 - 4 * abs(c2) ** 2 + 0j)) / (2 * c2)
    ok = estabilidade(s).incondicional
    return np.where(ok, gs, np.nan), np.where(ok, gl, np.nan)


def gt(s, gs, gl):
    """Ganho de transdução (eq. 11.93) para ΓS e ΓL (escalares ou por frequência)."""
    s11, s12, s21, s22 = s[..., 0, 0], s[..., 0, 1], s[..., 1, 0], s[..., 1, 1]
    num = (1 - abs(gs) ** 2) * abs(s21) ** 2 * (1 - abs(gl) ** 2)
    return num / abs((1 - s11 * gs) * (1 - s22 * gl) - s12 * s21 * gs * gl) ** 2


def gamas_entrada_saida(s, gs, gl):
    """(ΓIN, ΓOUT) do transistor com as terminações ΓL e ΓS (eqs. 11.34-11.35)."""
    s11, s12, s21, s22 = s[..., 0, 0], s[..., 0, 1], s[..., 1, 0], s[..., 1, 1]
    return s11 + s12 * s21 * gl / (1 - s22 * gl), s22 + s12 * s21 * gs / (1 - s11 * gs)


def z_de_gama(g, z0=50.0):
    return z0 * (1 + g) / (1 - g)


def interpolar(f, s, f0):
    """S em f0 por interpolação linear em módulo e fase (desenrolada) de cada termo."""
    out = np.empty((2, 2), dtype=complex)
    for i in range(2):
        for j in range(2):
            mag = np.interp(f0, f, abs(s[:, i, j]))
            fase = np.interp(f0, f, np.unwrap(np.angle(s[:, i, j])))
            out[i, j] = mag * np.exp(1j * fase)
    return out
