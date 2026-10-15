"""Análise e síntese quase-estática de linhas microstrip.

Modelo: Hammerstad & Jensen, "Accurate Models for Microstrip Computer-Aided
Design", IEEE MTT-S Int. Microwave Symp. Dig., 1980, pp. 407-409.
Hipóteses: espessura de metal nula, sem dispersão em frequência
(aproximação de primeira ordem também adotada no livro de F. Martín, cap. 2).
"""
from dataclasses import dataclass

import numpy as np

C0 = 299_792_458.0  # velocidade da luz no vácuo [m/s]
ETA0 = 376.730313668  # impedância do vácuo [ohm]


@dataclass
class Substrato:
    er: float = 10.2        # permissividade relativa
    h_mm: float = 1.27      # espessura do dielétrico [mm]
    t_um: float = 35.0      # espessura do cobre [um] (usada só na exportação 3D)
    tand: float = 0.0023    # tangente de perdas (usada só na exportação 3D)
    nome: str = "Rogers RO3010"


def analise(w_mm, sub: Substrato):
    """Retorna (Z0 [ohm], eps_eff) para largura(s) w_mm."""
    u = np.asarray(w_mm, dtype=float) / sub.h_mm
    er = sub.er
    a = (1 + np.log((u**4 + (u / 52) ** 2) / (u**4 + 0.432)) / 49
         + np.log(1 + (u / 18.1) ** 3) / 18.7)
    b = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053
    eeff = (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u) ** (-a * b)
    f = 6 + (2 * np.pi - 6) * np.exp(-((30.666 / u) ** 0.7528))
    z01 = ETA0 / (2 * np.pi) * np.log(f / u + np.sqrt(1 + (2 / u) ** 2))
    return z01 / np.sqrt(eeff), eeff


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
