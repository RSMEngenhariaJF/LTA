"""Elementos de sintonia e linhas metamateriais sintonizáveis.

Referência: F. Martín (2015), seções 5.2 e 5.3 (eqs. 5.1-5.12).
"""
import numpy as np

from . import abcd
from .crlh import CelulaCRLH

EPS0 = 8.8541878128e-12


# ---------------------------------------------------------------------------
# Dispositivos (5.2)
# ---------------------------------------------------------------------------
def varactor(v, cj0, phi=0.7, gama=0.5):
    """Capacitância de junção C = Cj0/(1 + V/φ)^γ (V reverso ≥ 0).

    γ = 1/2 para junção abrupta; o livro (nota 1, p. 361) indica largura de depleção
    proporcional à raiz m-ésima da tensão, 2 < m < 3, isto é, γ = 1/m entre 1/3 e 1/2.
    """
    v = np.asarray(v, dtype=float)
    if np.any(v < -phi):
        raise ValueError("Tensão reversa abaixo de -φ: o diodo conduziria.")
    return cj0 / (1 + v / phi) ** gama


def tensao_varactor(c, cj0, phi=0.7, gama=0.5):
    """Inversa de varactor(): tensão reversa que dá a capacitância c."""
    return phi * ((cj0 / np.asarray(c, dtype=float)) ** (1 / gama) - 1)


def hbv(v, cm, k, c0=0.0):
    """Varactor de heterobarreira (simétrico): C = C0 + CM·sech(kV) (livro, seção 5.4)."""
    return c0 + cm / np.cosh(k * np.asarray(v, dtype=float))


def mems_capacitancia(g, area):
    """C(g) = ε0·A/g (eq. 5.2, camada isolante desprezada)."""
    return EPS0 * area / np.asarray(g, dtype=float)


def mems_tensao(g, g0, k_mola, area):
    """Tensão que leva a ponte à altura g: V = √(2k·g²(g0 - g)/(ε0·A)) (eq. 5.4)."""
    g = np.asarray(g, dtype=float)
    return np.sqrt(2 * k_mola * g**2 * (g0 - g) / (EPS0 * area))


def mems_pull_in(g0, k_mola, area):
    """Tensão de colapso Vp = √(8k·g0³/(27ε0·A)), em g = 2g0/3 (eq. 5.5)."""
    return np.sqrt(8 * k_mola * g0**3 / (27 * EPS0 * area))


# ---------------------------------------------------------------------------
# Linha com SRR carregado por varactor (VLSRR), filtro rejeita-faixa (5.3.1.1)
# ---------------------------------------------------------------------------
def vlsrr_ceq(cl, cvar):
    """Capacitância do anel: CL (meia-borda esquerda) em série com o varactor (CR desprezada)."""
    return cl * cvar / (cl + cvar)


def vlsrr_tanque(m, lr, ceq, rs):
    """Tanque RLC refletido na linha (eqs. 5.6): Ls' = 4M²Ceqω0², Cs' = Lr/(4M²ω0²), R = 4M²ω0²/Rs.

    Q = ω0·R·Cs' = ω0·Lr/Rs (cresce com ω0: rejeição mais funda e estreita no topo da faixa).
    """
    w0 = 1 / np.sqrt(lr * ceq)
    return 4 * m**2 * ceq * w0**2, lr / (4 * m**2 * w0**2), 4 * m**2 * w0**2 / rs, w0


def vlsrr_abcd(f, l, c, m, lr, ceq, rs):
    """Célula π: C/2 - [L' = 2L - Ls' em série com o tanque RLC] - C/2 (Fig. 5.8a, transformada)."""
    ls, cs, r, _ = vlsrr_tanque(m, lr, ceq, rs)
    w = 2 * np.pi * np.asarray(f, dtype=float)
    tanque = 1 / (1 / r + 1j * w * cs + 1 / (1j * w * ls))
    zs = 1j * w * (2 * l - ls) + tanque
    y = 1j * w * c / 2
    return abcd.cascata(abcd.shunt(y), abcd.series(zs), abcd.shunt(y))


# ---------------------------------------------------------------------------
# CRLH sintonizável (5.3.2): defasador e antena leaky-wave
# ---------------------------------------------------------------------------
def crlh_sintonizada(c: CelulaCRLH, x):
    """Multiplica CL e LL pelo mesmo fator x (eqs. 5.9-5.11).

    Na linha balanceada (LR·CL = LL·CR), ZB fica constante em toda a sintonia; a
    transição passa a ω0/√x.
    """
    return CelulaCRLH(c.LR, c.CR, c.LL * x, c.CL * x, c.tipo)


def figura_merito(fases_graus, s21_db):
    """FoM = (φmáx - φmín) / máx(atenuação) [°/dB] sobre os estados de sintonia (eq. 5.12)."""
    perda = -np.min(s21_db)
    if perda <= 0:
        return np.inf
    return (np.max(fases_graus) - np.min(fases_graus)) / perda


def crlh_abcd_q(c: CelulaCRLH, f, q=np.inf):
    """ABCD da célula CRLH com perdas nos elementos sintonizáveis (CL e LL) de fator Q.

    Capacitor: Y = jωC(1 - j/Q); indutor: Z = jωL(1 - j/Q). Com Q = ∞, igual a c.abcd(f).
    """
    w = 2 * np.pi * np.asarray(f, dtype=float)
    k = 1 - 1j / q if np.isfinite(q) else 1.0
    zs = 1j * w * c.LR + 1 / (1j * w * c.CL * k)
    yp = 1j * w * c.CR + 1 / (1j * w * c.LL * k)
    return abcd.cascata(abcd.series(zs / 2), abcd.shunt(yp), abcd.series(zs / 2))
