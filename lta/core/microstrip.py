"""Análise e síntese de linhas microstrip: quase-estática, dispersão e perdas.

Modelo quase-estático: Hammerstad & Jensen, "Accurate Models for Microstrip
Computer-Aided Design", IEEE MTT-S Int. Microwave Symp. Dig., 1980, pp. 407-409
(espessura de metal nula, aproximação de primeira ordem também adotada no livro
de F. Martín, cap. 2).

Opcionais (funções separadas, não alteram analise/sintese):
  - dispersão de εeff: Kirschning & Jansen, Electron. Lett. 18(6), 1982, pp. 272-273;
    Z0(f) pela fórmula de Bianco et al. (1978), com a definição potência-corrente.
  - perdas: dielétrico, eq. 1.68 de Martín (2015); condutor, eq. 1.90 (W >= 2h) e
    Pucel/Gupta (Microstrip Lines and Slotlines, 1979) para tiras estreitas, com
    profundidade pelicular (eq. 1.79), espessura finita e rugosidade de Hammerstad.
"""
from dataclasses import dataclass

import numpy as np

C0 = 299_792_458.0  # velocidade da luz no vácuo [m/s]
ETA0 = 376.730313668  # impedância do vácuo [ohm]
MU0 = 4e-7 * np.pi    # permeabilidade do vácuo [H/m]
NP_DB = 20 / np.log(10)   # 1 Np = 8,686 dB


@dataclass
class Substrato:
    er: float = 10.2        # permissividade relativa
    h_mm: float = 1.27      # espessura do dielétrico [mm]
    t_um: float = 35.0      # espessura do cobre [um] (perdas de condutor e exportação 3D)
    tand: float = 0.0023    # tangente de perdas (perdas dielétricas e exportação 3D)
    nome: str = "Rogers RO3010"
    sigma: float = 5.8e7    # condutividade do metal [S/m] (cobre)
    rug_um: float = 0.0     # rugosidade RMS do cobre [um] (0 = superfície lisa)
    # True: Z0 e εeff incluem a espessura t do metal (Bahl e Garg; Garg, Bahl e Bianchi,
    # Microstrip Lines and Slotlines, 4ª ed., eqs. 2.124-2.127). Padrão False (t = 0), como no livro do Martín.
    incluir_t: bool = False


def _eeff_hj(u, er):
    a = (1 + np.log((u**4 + (u / 52) ** 2) / (u**4 + 0.432)) / 49
         + np.log(1 + (u / 18.1) ** 3) / 18.7)
    b = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053
    return (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u) ** (-a * b)


def _z_ar_hj(u):
    f = 6 + (2 * np.pi - 6) * np.exp(-((30.666 / u) ** 0.7528))
    return ETA0 / (2 * np.pi) * np.log(f / u + np.sqrt(1 + (2 / u) ** 2))


def largura_efetiva_t(w_mm, h_mm, t_mm):
    """We/h com espessura t (eq. 2.125): W/h + 1,25 t/(πh)·(1 + ln(4πW/t) ou ln(2h/t))."""
    u = np.asarray(w_mm, dtype=float) / h_mm
    if t_mm <= 0:
        return u
    return u + 1.25 * t_mm / (np.pi * h_mm) * np.where(u <= 1 / (2 * np.pi),
                                                       1 + np.log(4 * np.pi * np.asarray(w_mm) / t_mm),
                                                       1 + np.log(2 * h_mm / t_mm))


def analise(w_mm, sub: Substrato):
    """Retorna (Z0 [ohm], eps_eff) para largura(s) w_mm (Hammerstad-Jensen).

    Com sub.incluir_t, a impedância em ar usa a largura efetiva We (eq. 2.125) e
    εeff perde C = (εr - 1)/4,6·(t/h)/√(W/h) (eqs. 2.126-2.127).
    """
    u = np.asarray(w_mm, dtype=float) / sub.h_mm
    er = sub.er
    eeff = _eeff_hj(u, er)
    if sub.incluir_t and sub.t_um > 0:
        t = sub.t_um * 1e-3
        eeff = eeff - (er - 1) / 4.6 * (t / sub.h_mm) / np.sqrt(u)
        u = largura_efetiva_t(w_mm, sub.h_mm, t)
    return _z_ar_hj(u) / np.sqrt(eeff), eeff


def sintese(z0, sub: Substrato):
    """Largura [mm] que produz a impedância z0 (escalar ou array).

    Inversão numérica por bisseção em log(W/h); Z0 é monotônica em W.
    """
    z0 = np.atleast_1d(np.asarray(z0, dtype=float))
    lo = np.full_like(z0, np.log(1e-3))
    hi = np.full_like(z0, np.log(1e3))
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        zm, _ = analise(np.exp(mid) * sub.h_mm, sub)
        maior = zm > z0          # largura pequena demais -> aumentar
        lo = np.where(maior, mid, lo)
        hi = np.where(maior, hi, mid)
    w = np.exp(0.5 * (lo + hi)) * sub.h_mm
    return w if w.size > 1 else float(w[0])


# ---------------------------------------------------------------------------
# Dispersão em frequência
# ---------------------------------------------------------------------------
def _grade(w_mm, f):
    """Larguras e frequências em grade (Nf, Nw)."""
    w = np.atleast_1d(np.asarray(w_mm, dtype=float))[None, :]
    f = np.atleast_1d(np.asarray(f, dtype=float))[:, None]
    return w, f


def eeff_dispersivo(w_mm, sub: Substrato, f):
    """εeff(f) de Kirschning-Jansen, forma (Nf, Nw). Validade: εr <= 20, f·h <= 25 GHz·mm."""
    w, f = _grade(w_mm, f)
    u = w / sub.h_mm
    er = sub.er
    _, ee0 = analise(w, sub)
    fn = f * 1e-9 * sub.h_mm                       # GHz·mm
    p1 = 0.27488 + (0.6315 + 0.525 / (1 + 0.0157 * fn) ** 20) * u - 0.065683 * np.exp(-8.7513 * u)
    p2 = 0.33622 * (1 - np.exp(-0.03442 * er))
    p3 = 0.0363 * np.exp(-4.6 * u) * (1 - np.exp(-((fn / 38.7) ** 4.97)))
    p4 = 1 + 2.751 * (1 - np.exp(-((er / 15.916) ** 8)))
    p = p1 * p2 * ((0.1844 + p3 * p4) * fn) ** 1.5763
    return er - (er - ee0) / (1 + p)


def dispersao(w_mm, sub: Substrato, f):
    """(Z0(f), εeff(f)), forma (Nf, Nw).

    Z0(f) = Z0(0)·√(εeff(0)/εeff(f))·(εeff(f)-1)/(εeff(0)-1)  (Bianco et al., 1978).
    """
    w, _ = _grade(w_mm, f)
    z00, ee0 = analise(w, sub)
    ee = eeff_dispersivo(w_mm, sub, f)
    if sub.er <= 1:
        return np.broadcast_to(z00, ee.shape).copy(), ee
    return z00 * np.sqrt(ee0 / ee) * (ee - 1) / (ee0 - 1), ee


# ---------------------------------------------------------------------------
# Perdas
# ---------------------------------------------------------------------------
def profundidade_pelicular(f, sigma=5.8e7):
    """δp = √(2/(ωμσ)) [m] (eq. 1.79). Cobre: 2,06 µm a 1 GHz."""
    return np.sqrt(1 / (np.pi * np.asarray(f, dtype=float) * MU0 * sigma))


def resistencia_superficial(f, sub: Substrato):
    """Rs [Ω/□] = 1/(σδp) (eq. 1.83), corrigida para espessura finita e rugosidade.

    - Espessura: Rs = 1/(σ·δ·(1 - e^(-t/δ))) tende a 1/(σt) quando t << δ (corrente uniforme).
    - Rugosidade (Hammerstad): K = 1 + (2/π)·atan(1,4·(Δ/δ)²), Δ = rugosidade RMS.
    """
    d = profundidade_pelicular(np.maximum(f, 1.0), sub.sigma)
    t = max(sub.t_um, 1e-3) * 1e-6
    rs = 1 / (sub.sigma * d * (1 - np.exp(-t / d)))
    if sub.rug_um > 0:
        rs = rs * (1 + 2 / np.pi * np.arctan(1.4 * (sub.rug_um * 1e-6 / d) ** 2))
    return rs


def _largura_efetiva(w, h, t):
    """We com correção de espessura (Pucel; Gupta et al., 1979). Mesmas unidades de w."""
    return np.where(w / h >= 1 / (2 * np.pi),
                    w + t / np.pi * (1 + np.log(2 * h / t)),
                    w + t / np.pi * (1 + np.log(4 * np.pi * w / t)))


def alfa_condutor(w_mm, sub: Substrato, f, z0=None):
    """Atenuação por perdas no condutor αc [Np/m], forma (Nf, Nw).

    W/h >= 2: eq. 1.90 do livro. W/h < 2: expressões de Pucel (Gupta et al., 1979),
    citadas como ref. [12] no livro. z0: impedância a usar (padrão: quase-estática).
    """
    w, f = _grade(w_mm, f)
    h = sub.h_mm
    t = max(sub.t_um, 1e-3) * 1e-3                 # mm
    if z0 is None:
        z0, _ = analise(w, sub)
    rs = resistencia_superficial(f, sub)
    we = _largura_efetiva(w, h, t)
    u = w / h
    ue = we / h
    # tiras estreitas e médias (Pucel)
    fator_e = np.where(u <= 1 / (2 * np.pi),
                       np.log(4 * np.pi * w / t) + t / w,
                       np.log(2 * h / t) - t / h)
    estreita = (1 - (ue / 4) ** 2) * (1 + 1 / ue + fator_e / (np.pi * ue)) / (2 * np.pi)
    # tiras largas (eq. 1.90)
    a = ue + 2 / np.pi * np.log(2 * np.pi * np.e * (ue / 2 + 0.94))
    larga = (ue + ue / np.pi / (ue / 2 + 0.94)) / a**2 * (
        1 + 1 / ue + 1 / (np.pi * ue) * (np.log(2 * h / t + 1) - (1 - t / h) / (1 + t / (2 * h))))
    g = np.where(u >= 2, larga, estreita)
    return rs / (z0 * h * 1e-3) * g


def alfa_dieletrico(w_mm, sub: Substrato, f, eeff=None):
    """Atenuação por perdas no dielétrico αd [Np/m] (eq. 1.68), forma (Nf, Nw)."""
    w, f = _grade(w_mm, f)
    if eeff is None:
        _, eeff = analise(w, sub)
    er = sub.er
    q = er * (eeff - 1) / (eeff * (er - 1)) if er > 1 else np.ones_like(eeff)
    return np.pi * f / C0 * q * np.sqrt(eeff) * sub.tand


def parametros_linha(w_mm, sub: Substrato, f, perdas=False, dispersivo=False):
    """Z0, εeff e α = αc + αd [Np/m] por frequência e largura, forma (Nf, Nw).

    Sem as opções, reproduz o modelo quase-estático sem perdas (α = 0).
    """
    w, ff = _grade(w_mm, f)
    if dispersivo:
        z0, ee = dispersao(w_mm, sub, f)
    else:
        z0, ee = analise(w, sub)
        z0, ee = (np.broadcast_to(x, (ff.shape[0], w.shape[1])).copy() for x in (z0, ee))
    alfa = np.zeros_like(ee)
    if perdas:
        alfa = alfa_condutor(w_mm, sub, f, z0) + alfa_dieletrico(w_mm, sub, f, ee)
    return z0, ee, alfa
