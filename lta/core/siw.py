"""Guia de onda integrado ao substrato (SIW): projeto e modo TE10.

Referência: F. Martín (2015), seção 6.5 (eqs. 6.31-6.34), que resume as refs. [66-68].
O SIW se comporta como guia retangular preenchido com dielétrico de largura efetiva a;
só os modos TEn0 existem (as vias não sustentam corrente longitudinal).
"""
from dataclasses import dataclass

import numpy as np

from . import abcd
from .microstrip import C0, ETA0


@dataclass
class ProjetoSIW:
    a: float          # m, largura efetiva (eq. 6.34)
    w: float          # m, distância entre centros das fileiras de vias (eq. 6.31)
    d: float          # m, diâmetro das vias
    s: float          # m, passo das vias
    er: float
    fc10: float       # Hz
    d_max: float      # m (eq. 6.32)
    s_max: float      # m (eq. 6.33)

    @property
    def fc20(self):
        return 2 * self.fc10

    def beta(self, f):
        """β do TE10 [rad/m] (imaginário abaixo do corte)."""
        k = 2 * np.pi * np.asarray(f, dtype=float) * np.sqrt(self.er) / C0
        return np.sqrt(k**2 - (np.pi / self.a) ** 2 + 0j)

    def lambda_g(self, f):
        return 2 * np.pi / self.beta(f).real

    def z_te(self, f):
        """Impedância de onda do TE10: η/√εr / √(1 - (fc/f)²)."""
        return ETA0 / np.sqrt(self.er) / np.sqrt(1 - (self.fc10 / np.asarray(f, dtype=float)) ** 2 + 0j)

    def alfa_dieletrico(self, f, tand):
        """αd = k²·tanδ/(2β) [Np/m] do guia preenchido (acima do corte)."""
        k = 2 * np.pi * np.asarray(f, dtype=float) * np.sqrt(self.er) / C0
        return k**2 * tand / (2 * self.beta(f).real)


def projeto(fc10, er, d, s=None, n_superior=2):
    """Largura efetiva pelo corte (6.34), largura física (6.31) e limites das vias (6.32-6.33).

    n_superior: ordem do primeiro modo superior que limita a banda (TE20 → 2).
    """
    if fc10 <= 0 or er < 1 or d <= 0:
        raise ValueError("Use fc10 > 0, εr ≥ 1 e diâmetro de via > 0.")
    a = C0 / (2 * fc10 * np.sqrt(er))
    w = 0.5 * (a + np.sqrt((a + 0.54 * d) ** 2 - 0.4 * d**2)) + 0.27 * d
    d_max = 2 * a / (5 * np.sqrt(n_superior**2 - 1))
    s = 2 * d if s is None else s
    return ProjetoSIW(a, w, d, s, er, fc10, d_max, 2 * d)


def verificacoes(p: ProjetoSIW):
    """Avisos de projeto (eqs. 6.32-6.33 e regra de vazamento de Deslandes e Wu)."""
    av = []
    if p.d > p.d_max:
        av.append(f"Via de {p.d * 1e3:.2f} mm acima do limite da eq. 6.32 ({p.d_max * 1e3:.2f} mm).")
    if p.s > p.s_max:
        av.append(f"Passo {p.s * 1e3:.2f} mm > 2d = {p.s_max * 1e3:.2f} mm (eq. 6.33): vazamento de radiação.")
    if p.s <= p.d:
        av.append("Passo menor ou igual ao diâmetro: vias se sobrepõem.")
    return av


def abcd_siw(f, p: ProjetoSIW, comprimento, tand=0.0):
    """Trecho de SIW como linha equivalente do TE10 (Z = Z_TE), com perdas dielétricas opcionais."""
    f = np.asarray(f, dtype=float)
    beta = p.beta(f)
    alfa = p.alfa_dieletrico(f, tand) if tand else 0.0
    gama_l = (1j * beta + alfa) * comprimento
    return abcd.linha(p.z_te(f), -1j * gama_l)
