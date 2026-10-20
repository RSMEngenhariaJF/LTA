"""Linhas planares convencionais: microstrip, stripline, CPW, CPW com terra (GCPW), CPS e slotline.

Referências (modelos quase-estáticos):
  - microstrip: Hammerstad e Jensen (1980), já em microstrip.py, com dispersão, perdas e
    espessura opcional;
  - stripline centrada: Cohn (1954), solução exata por transformação conforme;
  - CPW em substrato de espessura h, GCPW e CPS: R. N. Simons, "Coplanar Waveguide
    Circuits, Components, and Systems", Wiley, 2001 (integrais elípticas);
  - espessura do metal em CPW/CPS (eqs. 7.110-7.114), perda no condutor das coplanares
    (eq. 7.126) e slotline (Janaswamy e Schaubert; Garg e Gupta, eqs. 5.36-5.43): R. Garg,
    I. Bahl e M. Bianchi, Microstrip Lines and Slotlines, 4ª ed., 2024, caps. 5 e 7.
Livro: F. Martín (2015), seções 1.2-1.5 (linhas planares, linhas carregadas, perdas).
"""
from dataclasses import dataclass

import numpy as np

from . import abcd
from .csrr import elipk
from .microstrip import C0, ETA0, Substrato, alfa_condutor, analise, resistencia_superficial

TIPOS = ("Microstrip", "Stripline", "CPW", "CPW com terra (GCPW)", "Coplanar strips (CPS)", "Slotline")
EPS0 = 8.8541878128e-12


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


def _delta_t(t, w_ref, b, cpw=True):
    """Δ das larguras efetivas com metal de espessura t (eqs. 7.111 e 7.113; b = largura total)."""
    if t <= 0:
        return 0.0
    b = min(b, 500 * t)                                 # validade: 0,002 < t/b ≤ 0,2
    if cpw:
        return t / np.pi * (4.089 + (0.9536 + 3.864e-3 * b / t) * np.log(4 * np.pi * w_ref / t))
    return t / np.pi * (5.188 + (0.9959 + 8.0808e-4 * b / t) * np.log(4 * np.pi * w_ref / t))


def _com_espessura(z0, eeff, ca_sup0, ca_sup_t):
    """Troca a capacitância em ar do semiplano superior (t = 0) pela com espessura e recalcula (Z0, εeff)."""
    ca0 = 1 / (C0 * z0 * np.sqrt(eeff))
    cd = (eeff - 1) * ca0
    ca = ca0 - ca_sup0 + ca_sup_t
    return 1 / (C0 * np.sqrt(ca * (ca + cd))), (ca + cd) / ca


def cpw(s, g, h, er, t=0.0):
    """CPW sobre substrato de espessura h, sem terra inferior (Simons, eqs. 2.29-2.32).

    t > 0 [mm]: larguras efetivas Se = S + Δ e We = G - Δ no semiplano superior (eqs. 7.110-7.114).
    """
    k0 = s / (s + 2 * g)
    k1 = _sinh_razao(np.pi * s / (4 * h), np.pi * (s + 2 * g) / (4 * h))
    eeff = 1 + (er - 1) / 2 * _razao(k1) / _razao(k0)
    z = 30 * np.pi / np.sqrt(eeff) / _razao(k0)
    if t > 0:
        d = _delta_t(t, g, s + 2 * g)
        d = min(d, 0.5 * g)                             # fora da faixa de ajuste da eq. 7.111
        ke = (s + d) / (s + d + 2 * (g - d))
        z, eeff = _com_espessura(z, eeff, 2 * EPS0 * _razao(k0), 2 * EPS0 * _razao(ke))
    return z, eeff


def gcpw(s, g, h, er, t=0.0):
    """CPW com plano de terra inferior (Simons, eqs. 2.38-2.40); t como em cpw()."""
    k0 = s / (s + 2 * g)
    k3 = np.tanh(np.pi * s / (4 * h)) / np.tanh(np.pi * (s + 2 * g) / (4 * h))
    r0, r3 = _razao(k0), _razao(k3)
    eeff = (1 + er * r3 / r0) / (1 + r3 / r0)
    z = 60 * np.pi / np.sqrt(eeff) / (r0 + r3)
    if t > 0:
        d = _delta_t(t, g, s + 2 * g)
        d = min(d, 0.5 * g)                             # fora da faixa de ajuste da eq. 7.111
        ke = (s + d) / (s + d + 2 * (g - d))
        z, eeff = _com_espessura(z, eeff, 2 * EPS0 * r0, 2 * EPS0 * _razao(ke))
    return z, eeff


def cps(s, w, h, er, t=0.0):
    """Coplanar strips: duas fitas de largura w separadas por s (Simons, eqs. 7.1-7.4).

    t > 0: We = W + Δ e Se = S - Δ no semiplano superior (eqs. 7.112-7.113).
    """
    k = s / (s + 2 * w)
    k1 = _sinh_razao(np.pi * s / (4 * h), np.pi * (s + 2 * w) / (4 * h))
    eeff = 1 + (er - 1) / 2 * _razao(k1) / _razao(k)
    z = 120 * np.pi / np.sqrt(eeff) * _razao(k)
    if t > 0:
        d = _delta_t(t, w, s + 2 * w, cpw=False)
        d = min(d, 0.5 * s)                             # fora da faixa de ajuste da eq. 7.113
        ke = (s - d) / (s - d + 2 * (w + d))
        z, eeff = _com_espessura(z, eeff, EPS0 / (2 * _razao(k)), EPS0 / (2 * _razao(ke)))
    return z, eeff


def slotline(w, h, er, f):
    """Slotline de fenda w [mm] em substrato h [mm] na frequência f [Hz] (dispersiva).

    Retorna (Z0s [Ω], εeff = (λ0/λs)²). Janaswamy e Schaubert (eqs. 5.40-5.43) para
    2,22 ≤ εr ≤ 9,8, W/λ0 ≤ 0,075 e 0,006 ≤ h/λ0 ≤ 0,06; Garg e Gupta (eqs. 5.36-5.39) para
    9,8 < εr ≤ 20 e 0,02 ≤ W/h ≤ 1. Fora disso o valor é extrapolado.
    """
    f = np.asarray(f, dtype=float)
    lam0 = C0 / f * 1e3                                      # mm
    u, hl, wl = w / h, h / lam0, w / lam0
    ln, lg = np.log, np.log10
    if er <= 9.8:
        if er <= 3.8:
            ls = (1.045 - 0.365 * ln(er) + 6.3 * u * er**0.945 / (238.64 + 100 * u)
                  - (0.148 - 8.81 * (er + 0.95) / (100 * er)) * ln(hl))
            z = (60 + 3.69 * np.sin((er - 2.22) * np.pi / 2.36) + 133.5 * ln(10 * er) * np.sqrt(wl)
                 + 2.81 * (1 - 0.011 * er * (4.48 - ln(er))) * u * ln(100 * hl)
                 + 131.1 * (1.028 - ln(er)) * np.sqrt(hl)
                 + 12.48 * (1 + 0.18 * ln(er)) * u / np.sqrt(er - 2.06 + 0.85 * u**2))
        else:
            ls = (0.9217 - 0.277 * ln(er) + 0.0322 * u * np.sqrt(er / (u + 0.435))
                  - 0.01 * ln(hl) * (4.6 - 3.65 / (er**2 * np.sqrt(wl) * (9.06 - 100 * wl))))
            z = (73.6 - 2.15 * er + (638.9 - 31.37 * er) * wl**0.6
                 + (36.23 * np.sqrt(er**2 + 41) - 225) * u / (u + 0.876 * er - 2)
                 + 0.51 * (er + 2.12) * u * ln(100 * hl) - 0.753 * er * hl / np.sqrt(wl))
    else:
        h100 = hl * 100
        if u < 0.2:
            ls = 0.923 - 0.448 * lg(er) + 0.2 * u - (0.29 * u + 0.047) * lg(h100)
            z = (72.62 - 35.29 * lg(er) + 50 * (u - 0.02) * (u - 0.1) / u + lg(u * 100) * (44.28 - 19.58 * lg(er))
                 - (0.32 * lg(er) - 0.11 + u * (1.07 * lg(er) + 1.44)) * (11.4 - 6.07 * lg(er) - h100) ** 2)
        else:
            ls = 0.987 - 0.483 * lg(er) + u * (0.111 - 0.0022 * er) - (0.121 + 0.094 * u - 0.0032 * er) * lg(h100)
            z = (113.19 - 53.55 * lg(er) + 1.25 * u * (114.59 - 51.88 * lg(er)) + 20 * (u - 0.2) * (1 - u)
                 - (0.15 + 0.23 * lg(er) + u * (-0.79 + 2.07 * lg(er)))
                 * (10.25 - 5 * lg(er) + u * (2.1 - 1.42 * lg(er)) - h100) ** 2)
    return z, 1 / ls**2


@dataclass
class Geometria:
    """Dimensões [mm]: w = largura (microstrip/stripline), fita central S (CPW), fitas (CPS) ou
    fenda (slotline); g = fenda (CPW/GCPW) ou separação entre fitas (CPS); b = distância entre
    planos (stripline)."""
    tipo: str
    w: float
    g: float = 0.0
    b: float = 0.0


def _t_mm(sub):
    return sub.t_um * 1e-3 if sub.incluir_t else 0.0


def analisar(geo: Geometria, sub: Substrato, f=None):
    """(Z0, εeff) da geometria; a slotline exige a frequência f [Hz] (escalar)."""
    if geo.w <= 0:
        raise ValueError("Largura deve ser positiva.")
    if geo.tipo == TIPOS[0]:
        z, e = analise(geo.w, sub)
        return float(z), float(e)
    if geo.tipo == TIPOS[1]:
        if geo.b <= 0:
            raise ValueError("Informe a distância entre os planos de terra (b).")
        return stripline(geo.w, geo.b, sub.er)
    if geo.tipo == TIPOS[5]:
        if f is None:
            raise ValueError("A slotline é dispersiva: informe a frequência.")
        z, e = slotline(geo.w, sub.h_mm, sub.er, float(np.ravel(f)[0]))
        return float(z), float(e)
    if geo.g <= 0:
        raise ValueError("Informe a fenda/separação (G).")
    t = _t_mm(sub)
    if geo.tipo == TIPOS[2]:
        return cpw(geo.w, geo.g, sub.h_mm, sub.er, t)
    if geo.tipo == TIPOS[3]:
        return gcpw(geo.w, geo.g, sub.h_mm, sub.er, t)
    if geo.tipo == TIPOS[4]:
        return cps(geo.g, geo.w, sub.h_mm, sub.er, t)
    raise ValueError(f"Tipo desconhecido: {geo.tipo}")


def sintetizar(tipo, z0, sub: Substrato, g=0.0, b=0.0, f=None):
    """Largura (W, S, fitas da CPS ou fenda da slotline) [mm] para Z0, com G ou b fixos (bisseção em log)."""
    def z_de(w):
        return analisar(Geometria(tipo, w, g, b), sub, f)[0]
    lo, hi = np.log(1e-3), np.log(200.0)
    if tipo == TIPOS[5]:
        if not f:
            raise ValueError("A slotline é dispersiva: informe a frequência.")
        lo, hi = np.log(0.02 * sub.h_mm), np.log(min(2 * sub.h_mm, 0.075 * C0 / f * 1e3))
    zlo, zhi = z_de(np.exp(lo)), z_de(np.exp(hi))
    cresce = zhi > zlo                                  # na slotline Z cresce com a fenda
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


def alfa_condutor_coplanar(geo: Geometria, sub: Substrato, f):
    """Perda no condutor [Np/m] de CPW/GCPW (a = S/2, b = S/2 + G) ou CPS (a = S/2, b = S/2 + W),
    eq. 7.126 (válida para t > 3δ e t << a, b - a)."""
    t = max(sub.t_um, 1e-3) * 1e-6
    if geo.tipo == TIPOS[4]:
        a, b = geo.g / 2 * 1e-3, (geo.g / 2 + geo.w) * 1e-3
    else:
        a, b = geo.w / 2 * 1e-3, (geo.w / 2 + geo.g) * 1e-3
    k = a / b
    kk = elipk(k) * elipk(np.sqrt(1 - k * k))
    _, ee = analisar(geo, sub)
    rs = resistencia_superficial(np.asarray(f, dtype=float), sub)
    termo = (1 / a * (np.pi + np.log(8 * np.pi * a * (1 - k) / (t * (1 + k))))
             + 1 / b * (np.pi + np.log(8 * np.pi * b * (1 - k) / (t * (1 + k)))))
    return rs * np.sqrt(ee) / (480 * np.pi * kk * (1 - k * k)) * termo


def gama(f, geo: Geometria, sub: Substrato, perdas=True):
    """(γ = α + jβ [1/m], Z0); perdas de condutor na microstrip (eq. 1.90/Pucel) e nas coplanares
    (eq. 7.126); a slotline varia Z0 e εeff com f e só tem perda dielétrica."""
    f = np.asarray(f, dtype=float)
    if geo.tipo == TIPOS[5]:
        z0, ee = slotline(geo.w, sub.h_mm, sub.er, f)
    else:
        z0, ee = analisar(geo, sub)
    beta = 2 * np.pi * f * np.sqrt(ee) / C0
    alfa = 0.0
    if perdas:
        alfa = alfa_dieletrico(f, sub.er, ee, sub.tand)
        if geo.tipo == TIPOS[0]:
            alfa = alfa + alfa_condutor(geo.w, sub, f)[:, 0]
        elif geo.tipo in TIPOS[2:5]:
            alfa = alfa + alfa_condutor_coplanar(geo, sub, f)
    return alfa + 1j * beta, z0


def abcd_trecho(f, geo, sub, comprimento_mm, perdas=True):
    g, z0 = gama(f, geo, sub, perdas)
    return abcd.linha(z0, -1j * g * comprimento_mm * 1e-3)


def abcd_stub(f, geo, sub, comprimento_mm, curto=False, perdas=True, raio_via=None):
    """Stub em derivação: Yin = tanh(γl)/Z0 (aberto) ou 1/(Z0·tanh(γl)) (curto).

    raio_via [mm]: o curto é a via real, Zt = R + jωL (Goldfarb e Pucel, eqs. 3.61-3.62).
    """
    g, z0 = gama(f, geo, sub, perdas)
    t = np.tanh(g * comprimento_mm * 1e-3)
    if curto and raio_via:
        from .descontinuidades import via
        lv, rv = via(raio_via, sub, f)
        zt = rv + 1j * 2 * np.pi * np.asarray(f, dtype=float) * lv
        return abcd.shunt((z0 + zt * t) / (z0 * (zt + z0 * t)))
    y = 1 / (z0 * t) if curto else t / z0
    return abcd.shunt(y)


def encurtamento_via(raio_mm, geo, sub, f0):
    """Quanto encurtar [mm] um stub em curto para compensar a indutância da via em f0."""
    from .descontinuidades import via
    lv, _ = via(raio_mm, sub, f0)
    z0, ee = analisar(geo, sub, f0)
    return float(np.arctan(2 * np.pi * f0 * lv / z0) * C0 / (2 * np.pi * f0 * np.sqrt(ee)) * 1e3)


def comprimento_eletrico(geo, sub, f, graus):
    """Comprimento físico [mm] para `graus` de fase em f."""
    _, ee = analisar(geo, sub, f)
    return np.radians(graus) * C0 / (2 * np.pi * f * np.sqrt(ee)) * 1e3
