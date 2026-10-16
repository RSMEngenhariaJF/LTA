"""Linha carregada com capacitores negativos (elementos não Foster): linha ENZ superluminal.

Referência: F. Martín (2015), seção 6.4.2 (Fig. 6.32).

Linha hospedeira TEM (Z, C' por metro) com capacitores negativos -Cn a cada d metros:
C'_ef = C' - Cn/d. Estável enquanto C'_ef > 0 (livro, p. 462); com L' = Z²C' fixo,
ε_ef = c²·L'·C'_ef < 1 e vp = vg > c em banda larga (sem dispersão no limite contínuo).
"""
from dataclasses import dataclass

import numpy as np

from . import abcd
from .microstrip import C0


@dataclass
class LinhaENZ:
    z_host: float     # Ω
    c_host: float     # F/m
    c_neg: float      # F (módulo do capacitor negativo)
    passo: float      # m entre capacitores

    @property
    def l_por_m(self):
        return self.z_host**2 * self.c_host

    @property
    def c_ef(self):
        return self.c_host - self.c_neg / self.passo

    @property
    def estavel(self):
        return self.c_ef > 0

    @property
    def z_ef(self):
        return np.sqrt(self.l_por_m / self.c_ef) if self.estavel else np.nan

    @property
    def eps_ef(self):
        """c²·L'·C'_ef: permissividade equivalente relativa a uma linha TEM no ar."""
        return C0**2 * self.l_por_m * self.c_ef

    @property
    def v_sobre_c(self):
        return 1 / np.sqrt(self.eps_ef) if self.estavel else np.nan

    def abcd(self, f, n):
        """n células: trecho de linha (passo) + derivação de -Cn (ideal, independente de f)."""
        f = np.asarray(f, dtype=float)
        w = 2 * np.pi * f
        theta = w * self.passo * np.sqrt(self.l_por_m * self.c_host)
        meia = abcd.linha(self.z_host, theta / 2)
        cel = abcd.cascata(meia, abcd.shunt(-1j * w * self.c_neg), meia)
        return abcd.potencia(cel, n)


def capacitor_para_eps(z_host, c_host, passo, eps_alvo):
    """|Cn| que leva a ε_ef desejada: C'_ef = ε_alvo/(c²·L')."""
    l = z_host**2 * c_host
    return (c_host - eps_alvo / (C0**2 * l)) * passo
