"""Atenuadores resistivos T e Pi (pads), com impedâncias diferentes nos dois lados.

Referência: M. Steer, Microwave and RF Design: A Systems Approach, SciTech, 2010, seção 7.3.2:
  - T: eqs. 7.2-7.6; Pi: eqs. 7.7-7.11; K = Pentrada/Psaída (7.12); Kmín com Z01 ≠ Z02 (7.14).
Versões balanceadas: no T, R1 e R2 se dividem em duas metades (uma em cada fio); no Pi, R3 se
divide em duas metades. Também: valores comerciais E24/E96 mais próximos e o efeito do
arredondamento, e a potência dissipada em cada resistor.
"""
from dataclasses import dataclass, field

import numpy as np

from . import abcd

TOPOLOGIAS = ("T", "Pi")
SERIES = {"E24": (1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6,
                  6.2, 6.8, 7.5, 8.2, 9.1),
          "E96": tuple(round(10 ** (i / 96), 2) for i in range(96))}


def k_minimo(z01, z02):
    """Menor K (linear) de um pad casado nos dois lados (eq. 7.14, com r = Zmaior/Zmenor)."""
    r = max(z01, z02) / min(z01, z02)
    return float(2 * r - 1 + 2 * np.sqrt(r * (r - 1)))


@dataclass
class Atenuador:
    topologia: str
    z01: float
    z02: float
    k_db: float
    r: tuple                   # (R1, R2, R3) [Ω]
    texto: list = field(default_factory=list)

    def abcd(self, f, r=None):
        r1, r2, r3 = self.r if r is None else r
        n = np.shape(np.asarray(f, dtype=float))
        um = np.ones(n)
        if self.topologia == "T":
            return abcd.cascata(abcd.series(r1 * um), abcd.shunt(um / r3), abcd.series(r2 * um))
        return abcd.cascata(abcd.shunt(um / r1), abcd.series(r3 * um), abcd.shunt(um / r2))

    def s(self, f, r=None):
        return abcd.abcd_para_s2(self.abcd(f, r), self.z01, self.z02)

    def potencias(self, p_entrada=1.0):
        """Potência dissipada em (R1, R2, R3) e entregue à carga, para p_entrada [W] casada."""
        r1, r2, r3 = self.r
        k = 10 ** (self.k_db / 10)
        p_saida = p_entrada / k
        v1, v2 = np.sqrt(p_entrada * self.z01), np.sqrt(p_saida * self.z02)
        if self.topologia == "T":
            p1, p2 = (v1 / self.z01) ** 2 * r1, (v2 / self.z02) ** 2 * r2
        else:
            p1, p2 = v1**2 / r1, v2**2 / r2
        return p1, p2, p_entrada - p_saida - p1 - p2, p_saida


def projeto(topologia, z01, z02, k_db):
    """Pad de K dB entre Z01 (porta 1) e Z02 (porta 2)."""
    if z01 <= 0 or z02 <= 0:
        raise ValueError("As impedâncias devem ser positivas.")
    k = 10 ** (k_db / 10)
    kmin = k_minimo(z01, z02)
    if k <= 1 or k < kmin * (1 - 1e-9):
        raise ValueError(f"Com Z01 = {z01:g} Ω e Z02 = {z02:g} Ω a atenuação mínima de um pad casado é "
                         f"{10 * np.log10(kmin):.3f} dB (eq. 7.14).")
    rk = np.sqrt(k)
    if topologia == "T":
        r3 = 2 * np.sqrt(k * z01 * z02) / (k - 1)
        r1 = (z01 * (k + 1) - 2 * np.sqrt(k * z01 * z02)) / (k - 1)
        r2 = (z02 * (k + 1) - 2 * np.sqrt(k * z01 * z02)) / (k - 1)
        eqs = "eqs. 7.2-7.4"
    elif topologia == "Pi":
        r3 = (k - 1) / 2 * np.sqrt(z01 * z02 / k)
        den1 = (k + 1) * np.sqrt(z02) - 2 * rk * np.sqrt(z01)
        den2 = (k + 1) * np.sqrt(z01) - 2 * rk * np.sqrt(z02)
        r1 = z01 * (k - 1) * np.sqrt(z02) / den1 if abs(den1) > 1e-12 else np.inf
        r2 = z02 * (k - 1) * np.sqrt(z01) / den2 if abs(den2) > 1e-12 else np.inf
        eqs = "eqs. 7.7-7.9"
    else:
        raise ValueError(f"Topologia desconhecida: {topologia}")
    a = Atenuador(topologia, float(z01), float(z02), float(k_db), (float(r1), float(r2), float(r3)))
    nomes = ("R1 (série, lado 1)", "R2 (série, lado 2)", "R3 (derivação)") if topologia == "T" else \
        ("R1 (derivação, lado 1)", "R2 (derivação, lado 2)", "R3 (série)")
    a.texto = [f"Pad {topologia} de {k_db:g} dB (K = {k:.4f}), Z01 = {z01:g} Ω, Z02 = {z02:g} Ω ({eqs}); "
               f"Kmín = {10 * np.log10(kmin):.3f} dB"]
    a.texto += [f"  {n} = {_fmt(x)}" for n, x in zip(nomes, a.r)]
    if topologia == "T":
        a.texto.append(f"  Balanceado: R1/2 = {_fmt(r1 / 2)} e R2/2 = {_fmt(r2 / 2)} em cada fio, R3 = {_fmt(r3)}")
    else:
        a.texto.append(f"  Balanceado: R3/2 = {_fmt(r3 / 2)} em cada fio, R1 = {_fmt(r1)}, R2 = {_fmt(r2)}")
    return a


def _fmt(r):
    if not np.isfinite(r):
        return "∞ (aberto)"
    return f"{r / 1e3:.4g} kΩ" if r >= 1e3 else f"{r:.4g} Ω"


def comercial(r, serie="E24"):
    """Valor da série (E24/E96) mais próximo de r (escala logarítmica); r ≤ 0 ou ∞ ficam como estão."""
    if not np.isfinite(r) or r <= 0:
        return r
    dec = 10 ** np.floor(np.log10(r))
    cand = np.array([v * d for d in (dec / 10, dec, dec * 10) for v in SERIES[serie]])
    return float(cand[np.argmin(np.abs(np.log(cand / r)))])


def efeito_arredondamento(a: Atenuador, serie="E24"):
    """(resistores comerciais, atenuação [dB], |S11| [dB], |S22| [dB]) com os valores da série."""
    rc = tuple(comercial(x, serie) for x in a.r)
    s = a.s(np.array([1.0]), rc)[0]
    return rc, -20 * np.log10(abs(s[1, 0])), 20 * np.log10(abs(s[0, 0]) + 1e-15), 20 * np.log10(abs(s[1, 1]) + 1e-15)
