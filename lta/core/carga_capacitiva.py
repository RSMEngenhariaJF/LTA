"""Linha com carga capacitiva periódica (efeito slow-wave).

Referência: F. Martín (2015), seção 2.5.2, eqs. 2.116-2.126.
Modelo concentrado da célula: indutância L e capacitância C por seção da
linha não carregada, mais a capacitância de carga Cls em derivação.
"""
from dataclasses import dataclass

import numpy as np

from . import abcd
from .microstrip import C0, Substrato, analise, sintese


@dataclass
class ProjetoCarga:
    L: float      # H por célula
    C: float      # F por célula (linha não carregada)
    Cls: float    # F, capacitância de carga

    @property
    def zo(self):
        """Impedância da linha não carregada (eq. 2.123)."""
        return np.sqrt(self.L / self.C)

    @property
    def zb(self):
        """Impedância de Bloch em baixa frequência (eq. 2.118)."""
        return np.sqrt(self.L / (self.C + self.Cls))

    @property
    def fc(self):
        """Frequência de corte da 1a banda (eq. 2.117)."""
        return 1 / (np.pi * np.sqrt(self.L * (self.C + self.Cls)))

    @property
    def fs(self):
        """Início da 1a banda espúria (eq. 2.119)."""
        return 1 / (2 * np.sqrt(self.L * self.C))

    @property
    def fator_slow_wave(self):
        """v_p(não carregada) / v_p(carregada) em baixa frequência (eq. 2.124)."""
        return np.sqrt((self.C + self.Cls) / self.C)

    def fase_por_celula(self, f):
        """Defasagem por célula em baixa frequência [rad] (eq. 2.126)."""
        return 2 * np.pi * f * np.sqrt(self.L * (self.C + self.Cls))


def projeto_por_fc(fc, zb, fs):
    """Dados fc, ZB e fs, obtém L, C e Cls (inversão das eqs. 2.117-2.119)."""
    L = zb / (np.pi * fc)
    ct = 1 / (np.pi * fc * zb)
    C = 1 / (4 * fs**2 * L)
    return _valida(ProjetoCarga(L, C, ct - C))


def projeto_por_fase(fase_graus, f0, zb, fs):
    """Dados a fase por célula em f0, ZB e fs (procedimento do ex. da Fig. 2.34)."""
    raiz_lct = np.deg2rad(fase_graus) / (2 * np.pi * f0)
    L = zb * raiz_lct
    ct = raiz_lct / zb
    C = 1 / (4 * fs**2 * L)
    return _valida(ProjetoCarga(L, C, ct - C))


def _valida(p):
    if p.Cls <= 0:
        raise ValueError("Especificação inviável: Cls <= 0. Aumente fs ou reduza fc.")
    return p


def celula_distribuida(f, p: ProjetoCarga):
    """Meia linha + Cls em derivação + meia linha (simétrica). kl = w*sqrt(LC)."""
    kl = 2 * np.pi * np.asarray(f) * np.sqrt(p.L * p.C)
    meia = abcd.linha(p.zo, kl / 2)
    return abcd.cascata(meia, abcd.shunt(1j * 2 * np.pi * np.asarray(f) * p.Cls), meia)


def celula_concentrada(f, p: ProjetoCarga):
    """Célula T: L/2 - (C+Cls) - L/2 (Fig. 2.27b)."""
    return abcd.celula(f, "T", ls=p.L, cp=p.C + p.Cls)


@dataclass
class GeometriaCarga:
    w_linha: float    # mm, linha de alta impedância (Zo)
    eeff_linha: float
    periodo: float    # mm
    w_patch: float    # mm
    l_patch: float    # mm
    eeff_patch: float
    w_porta: float    # mm, linha de alimentação


def geometria(p: ProjetoCarga, sub: Substrato, z_patch=20.0, z_porta=50.0):
    """Estimativa inicial das dimensões físicas em microstrip.

    - Período: sqrt(LC) = l*sqrt(eeff)/c da linha de alta impedância.
    - Patch: trecho largo de impedância z_patch, cuja capacitância excedente
      em relação à linha substituída vale Cls:
        Cls ~ l_p*(sqrt(eeff_p)/Zp - sqrt(eeff)/Zo)/c
    É ponto de partida para otimização em simulador EM (como feito no livro).
    """
    w_l = sintese(p.zo, sub)
    _, ee_l = analise(w_l, sub)
    periodo = C0 * np.sqrt(p.L * p.C) / np.sqrt(ee_l) * 1e3
    w_p = sintese(z_patch, sub)
    _, ee_p = analise(w_p, sub)
    c_por_m = np.sqrt(ee_p) / (z_patch * C0) - np.sqrt(ee_l) / (p.zo * C0)
    if c_por_m <= 0:
        raise ValueError("A impedância do patch deve ser menor que Zo da linha.")
    l_patch = p.Cls / c_por_m * 1e3
    return GeometriaCarga(float(w_l), float(ee_l), float(periodo), float(w_p),
                          float(l_patch), float(ee_p), float(sintese(z_porta, sub)))
