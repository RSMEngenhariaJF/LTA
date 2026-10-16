"""Linhas planares convencionais: microstrip, stripline, CPW, CPW com terra (GCPW) e CPS.

Referências (modelos quase-estáticos, metal de espessura nula):
  - microstrip: Hammerstad e Jensen (1980), já em microstrip.py, com dispersão e perdas;
  - stripline centrada: Cohn (1954), solução exata por transformação conforme;
  - CPW em substrato de espessura h, GCPW e CPS: R. N. Simons, "Coplanar Waveguide
    Circuits, Components, and Systems", Wiley, 2001 (integrais elípticas).
Livro: F. Martín (2015), seções 1.2-1.5 (linhas planares, linhas carregadas, perdas).
"""
from dataclasses import dataclass

import numpy as np

from . import abcd
from .csrr import elipk
from .microstrip import C0, ETA0, Substrato, alfa_condutor, analise

TIPOS = ("Microstrip", "Stripline", "CPW", "CPW com terra (GCPW)", "Coplanar strips (CPS)")


def _razao(k):
    """K(k)/K(k')."""
    return elipk(k) / elipk(np.sqrt(1 - k * k))


def _sinh_razao(a, b):
    """sinh(a)/sinh(b) sem estouro para argumentos grandes (a < b)."""
    return np.exp(a - b) * (-np.expm1(-2 * a)) / (-np.expm1(-2 * b))


# ---------------------------------------------------------------------------
# Análise (dimensões em mm)
# ---------------------------------------------------------------------------
def stripline(w, b, er):
    """Stripline centrada (Cohn): Z0 = η0/(4√εr)·K(k)/K(k'), k = sech(πW/2b). Retorna (Z0, εeff)."""
    x = min(np.pi * w / (2 * b), 700.0)
    k = 1 / np.cosh(x)
    return ETA0 / (4 * np.sqrt(er)) * _razao(k), er


def cpw(s, g, h, er):
    """CPW sobre substrato de espessura h, sem terra inferior (Simons, eqs. 2.29-2.32)."""
    k0 = s / (s + 2 * g)
    k1 = _sinh_razao(np.pi * s / (4 * h), np.pi * (s + 2 * g) / (4 * h))
    eeff = 1 + (er - 1) / 2 * _razao(k1) / _razao(k0)
    return 30 * np.pi / np.sqrt(eeff) / _razao(k0), eeff


def gcpw(s, g, h, er):
    """CPW com plano de terra inferior (Simons, eqs. 2.38-2.40)."""
    k0 = s / (s + 2 * g)
    k3 = np.tanh(np.pi * s / (4 * h)) / np.tanh(np.pi * (s + 2 * g) / (4 * h))
    r0, r3 = _razao(k0), _razao(k3)
    eeff = (1 + er * r3 / r0) / (1 + r3 / r0)
    return 60 * np.pi / np.sqrt(eeff) / (r0 + r3), eeff


def cps(s, w, h, er):
    """Coplanar strips: duas fitas de largura w separadas por s (Simons, eqs. 7.1-7.4)."""
    k = s / (s + 2 * w)
    k1 = _sinh_razao(np.pi * s / (4 * h), np.pi * (s + 2 * w) / (4 * h))
    eeff = 1 + (er - 1) / 2 * _razao(k1) / _razao(k)
    return 120 * np.pi / np.sqrt(eeff) * _razao(k), eeff


@dataclass
class Geometria:
    """Dimensões [mm]: w = largura (microstrip/stripline) ou fita central S (CPW) ou fitas (CPS);
    g = fenda (CPW/GCPW) ou separação entre fitas (CPS); b = distância entre planos (stripline)."""
    tipo: str
    w: float
    g: float = 0.0
    b: float = 0.0


def analisar(geo: Geometria, sub: Substrato):
    """(Z0, εeff) da geometria."""
    if geo.w <= 0:
        raise ValueError("Largura deve ser positiva.")
    if geo.tipo == TIPOS[0]:
        z, e = analise(geo.w, sub)
        return float(z), float(e)
    if geo.tipo == TIPOS[1]:
        if geo.b <= 0:
            raise ValueError("Informe a distância entre os planos de terra (b).")
        return stripline(geo.w, geo.b, sub.er)
    if geo.g <= 0:
        raise ValueError("Informe a fenda/separação (G).")
    if geo.tipo == TIPOS[2]:
        return cpw(geo.w, geo.g, sub.h_mm, sub.er)
    if geo.tipo == TIPOS[3]:
        return gcpw(geo.w, geo.g, sub.h_mm, sub.er)
    if geo.tipo == TIPOS[4]:
        return cps(geo.g, geo.w, sub.h_mm, sub.er)
    raise ValueError(f"Tipo desconhecido: {geo.tipo}")


def sintetizar(tipo, z0, sub: Substrato, g=0.0, b=0.0):
    """Largura (W, S ou largura das fitas) [mm] para Z0, com G ou b fixos (bisseção em log)."""
    def z_de(w):
        return analisar(Geometria(tipo, w, g, b), sub)[0]
    lo, hi = np.log(1e-3), np.log(200.0)
    zlo, zhi = z_de(np.exp(lo)), z_de(np.exp(hi))
    cresce = zhi > zlo                                  # em todas as linhas Z cai com a largura
    if not min(zlo, zhi) < z0 < max(zlo, zhi):
        raise ValueError(f"Z0 = {z0:g} Ω fora da faixa realizável ({min(zlo, zhi):.1f} a {max(zlo, zhi):.1f} Ω) "
                         "com esses parâmetros.")
    for _ in range(80):
        m = 0.5 * (lo + hi)
        if (z_de(np.exp(m)) > z0) != cresce:
            lo = m
        else:
            hi = m
    return float(np.exp(0.5 * (lo + hi)))


# ---------------------------------------------------------------------------
# Propagação e elementos
# ---------------------------------------------------------------------------
def alfa_dieletrico(f, er, eeff, tand):
    """αd = π f/c · εr(εeff - 1)/(√εeff(εr - 1)) · tanδ [Np/m] (eq. 1.68; vale para linhas quase-TEM)."""
    q = er * (eeff - 1) / (eeff * (er - 1)) if er > 1 else np.sqrt(eeff)
    return np.pi * np.asarray(f, dtype=float) / C0 * q * np.sqrt(eeff) * tand


def gama(f, geo: Geometria, sub: Substrato, perdas=True):
    """γ = α + jβ [1/m]; perdas de condutor só na microstrip (eq. 1.90/Pucel)."""
    z0, ee = analisar(geo, sub)
    f = np.asarray(f, dtype=float)
    beta = 2 * np.pi * f * np.sqrt(ee) / C0
    alfa = 0.0
    if perdas:
        alfa = alfa_dieletrico(f, sub.er, ee, sub.tand)
        if geo.tipo == TIPOS[0]:
            alfa = alfa + alfa_condutor(geo.w, sub, f)[:, 0]
    return alfa + 1j * beta, z0


def abcd_trecho(f, geo, sub, comprimento_mm, perdas=True):
    g, z0 = gama(f, geo, sub, perdas)
    return abcd.linha(z0, -1j * g * comprimento_mm * 1e-3)


def abcd_stub(f, geo, sub, comprimento_mm, curto=False, perdas=True):
    """Stub em derivação: Yin = tanh(γl)/Z0 (aberto) ou 1/(Z0·tanh(γl)) (curto)."""
    g, z0 = gama(f, geo, sub, perdas)
    t = np.tanh(g * comprimento_mm * 1e-3)
    y = 1 / (z0 * t) if curto else t / z0
    return abcd.shunt(y)


def comprimento_eletrico(geo, sub, f, graus):
    """Comprimento físico [mm] para `graus` de fase em f."""
    _, ee = analisar(geo, sub)
    return np.radians(graus) * C0 / (2 * np.pi * f * np.sqrt(ee)) * 1e3
