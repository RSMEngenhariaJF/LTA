"""Correções físicas da microstrip: descontinuidades, limites do substrato e fórmulas do Pozar.

Referências:
  - extremidade aberta: M. Kirschning, R. H. Jansen e N. H. L. Koster, "Accurate model for
    open end effect of microstrip lines", Electron. Lett. 17(3), 1981, pp. 123-125;
  - gap série: P. Garg e I. J. Bahl, "Microstrip discontinuities", Int. J. Electronics 45(1),
    1978, pp. 81-87 (modelo π: Cg em série, Cp em derivação; 0,5 ≤ W/h ≤ 2, 2,5 ≤ εr ≤ 15);
  - ondas de superfície e modos superiores: D. M. Pozar, Microwave Engineering, 4ª ed., Wiley,
    2012, seções 3.6 e 3.8 (eq. 3.174 e Ex. 3.4);
  - fórmulas fechadas da microstrip: Pozar, eqs. 3.195-3.199 (referência de comparação);
  - degrau, curva, chanfro, via, entalhe, fT e Q de radiação: R. Garg, I. Bahl e M. Bianchi,
    Microstrip Lines and Slotlines, 4ª ed., Artech House, 2024, seções 2.4.6, 3.4 e 3.8.
"""
import numpy as np

from .microstrip import C0, Substrato, analise, resistencia_superficial


# ---------------------------------------------------------------------------
# Extremidade aberta
# ---------------------------------------------------------------------------
def extensao_ponta_aberta(w_mm, sub: Substrato):
    """Δl [mm] da extremidade aberta (Kirschning-Jansen-Koster): a ponta irradia como se a
    linha fosse Δl mais longa. Validade: 0,01 ≤ W/h ≤ 100, εr ≤ 128 (erro < 2,5%)."""
    u = np.asarray(w_mm, dtype=float) / sub.h_mm
    _, ee = analise(w_mm, sub)
    er = sub.er
    x1 = 0.434907 * (ee**0.81 + 0.26) / (ee**0.81 - 0.189) * (u**0.8544 + 0.236) / (u**0.8544 + 0.87)
    x2 = 1 + u**0.371 / (2.358 * er + 1)
    x3 = 1 + 0.5274 * np.arctan(0.084 * u ** (1.9413 / x2)) / ee**0.9236
    x4 = 1 + 0.0377 * np.arctan(0.067 * u**1.456) * (6 - 5 * np.exp(0.036 * (1 - er)))
    x5 = 1 - 0.218 * np.exp(-7.5 * u)
    dl = sub.h_mm * x1 * x3 * x5 / x4
    return float(dl) if np.ndim(dl) == 0 else dl


def capacitancia_ponta(w_mm, sub: Substrato):
    """Capacitância equivalente da ponta aberta [F]: C = Δl·√εeff/(c·Z0)."""
    z0, ee = analise(w_mm, sub)
    return float(extensao_ponta_aberta(w_mm, sub) * 1e-3 * np.sqrt(ee) / (C0 * z0))


# ---------------------------------------------------------------------------
# Gap série (modelo π)
# ---------------------------------------------------------------------------
def gap(w_mm, s_mm, sub: Substrato):
    """(Cg, Cp) [F] do gap simétrico de largura s numa linha de largura W (Garg-Bahl).

    Cg é a capacitância série e Cp a capacitância de cada lado para o terra. Fora da faixa
    de validade (0,1 ≤ s/W ≤ 1; 0,5 ≤ W/h ≤ 2) o valor é extrapolado: verificar em EM.
    """
    if w_mm <= 0 or s_mm <= 0:
        raise ValueError("Largura e gap devem ser positivos.")
    u = w_mm / sub.h_mm
    r = s_mm / w_mm
    lg = np.log10(u)
    mo = u * (0.619 * lg - 0.3853)
    ko = 4.26 - 1.453 * lg
    if r <= 0.3:
        me, ke = 0.8675, 2.043 * u**0.12
    else:
        me, ke = 1.565 / u**0.16 - 1, 1.97 - 0.03 / u
    w_m = w_mm * 1e-3
    c_impar = (sub.er / 9.6) ** 0.8 * r**mo * np.exp(ko) * 1e-12 * w_m
    c_par = 12 * (sub.er / 9.6) ** 0.9 * r**me * np.exp(ke) * 1e-12 * w_m
    cp = c_par / 2
    cg = (c_impar - c_par / 2) / 2
    return float(max(cg, 0.0)), float(cp)


def gap_para_capacitancia(w_mm, cg_alvo, sub: Substrato):
    """Largura do gap s [mm] que dá a capacitância série cg_alvo [F] (bisseção em log s)."""
    lo, hi = np.log(0.005 * w_mm), np.log(5.0 * w_mm)
    c_lo, c_hi = gap(w_mm, np.exp(lo), sub)[0], gap(w_mm, np.exp(hi), sub)[0]
    if not c_hi < cg_alvo < c_lo:
        raise ValueError(f"Cg = {cg_alvo * 1e12:.4g} pF fora do alcance de um gap simples nesta linha "
                         f"({c_hi * 1e12:.3g} a {c_lo * 1e12:.3g} pF): use interdigital ou capacitor SMD.")
    for _ in range(80):
        m = 0.5 * (lo + hi)
        if gap(w_mm, np.exp(m), sub)[0] > cg_alvo:
            lo = m
        else:
            hi = m
    return float(np.exp(0.5 * (lo + hi)))


# ---------------------------------------------------------------------------
# Limites do substrato (modos indesejados)
# ---------------------------------------------------------------------------
def frequencia_ft(sub: Substrato):
    """fT [Hz] de acoplamento forte entre o modo quase-TEM e a onda de superfície TM0
    (Garg, Bahl e Bianchi, eq. 2.140): fT = 150/(πh)·√(2/(εr - 1))·atan(εr) GHz, h em mm."""
    return 150 / (np.pi * sub.h_mm) * np.sqrt(2 / (sub.er - 1)) * np.arctan(sub.er) * 1e9


def q_radiacao(z0, sub: Substrato, f):
    """Q de radiação de um ressonador λ/2 em aberto (eq. 2.142): Qr = 3εr·Z0·λ0²/(32·η·h²)."""
    lam0 = C0 / np.asarray(f, dtype=float)
    return 3 * sub.er * z0 * lam0**2 / (32 * 120 * np.pi * (sub.h_mm * 1e-3) ** 2)


# ---------------------------------------------------------------------------
# Degrau de largura, curva, chanfro, via e entalhe (Garg, Bahl e Bianchi, cap. 3)
# ---------------------------------------------------------------------------
def degrau(w1_mm, w2_mm, sub: Substrato):
    """Δl [mm] do degrau de W1 (larga) para W2 (estreita): a linha larga fica Δl mais longa e a
    estreita Δl mais curta, Δl = Δloc(W1)·(1 - W2/W1) (Hammerstad e Bekkadal, seção 3.4.3)."""
    wl, we = max(w1_mm, w2_mm), min(w1_mm, w2_mm)
    return float(extensao_ponta_aberta(wl, sub) * (1 - we / wl))


def curva(w_mm, sub: Substrato):
    """(Cb [F], Lb [H]) da curva de 90° sem compensação (eqs. 3.57-3.58; 2,5 ≤ εr ≤ 15,
    0,1 ≤ W/h ≤ 5 para C e 0,5 ≤ W/h ≤ 2 para L). Lb é o excesso em relação às linhas de centro
    (pode ser negativo)."""
    u = w_mm / sub.h_mm
    er = sub.er
    w_m = w_mm * 1e-3
    if u < 1:
        cb = ((14 * er + 12.5) * u - (1.83 * er - 2.25)) / np.sqrt(u) + 0.02 * er / u
    else:
        cb = (9.5 * er + 1.25) * u + 5.2 * er + 7
    lb = 100 * (4 * np.sqrt(u) - 4.21)
    return float(cb * 1e-12 * w_m), float(lb * 1e-9 * w_m)


def chanfro_otimo(w_mm, sub: Substrato):
    """Chanfro ótimo M [%] da curva de 90° (eq. 3.115, W/h ≥ 0,25, εr ≤ 25): M = 52 + 65·exp(-1,35W/h).
    M = X/d·100, com d a diagonal do quadrado do canto e X o trecho cortado a partir do canto externo."""
    return float(52 + 65 * np.exp(-1.35 * w_mm / sub.h_mm))


def encurtamento_chanfro(w_mm, sub: Substrato, f):
    """Δb [mm]: a curva chanfrada é eletricamente Δb mais curta que as linhas de centro
    (eqs. 3.116-3.117): Δb/D = 0,16·(2 - (f/fp)²), D = η0·h/(Z0·√εeff), fp = 0,4·Z0/h GHz."""
    z0, ee = analise(w_mm, sub)
    d = 120 * np.pi * sub.h_mm / (float(z0) * np.sqrt(float(ee)))
    fp = 0.4 * float(z0) / sub.h_mm * 1e9
    return float(d * 0.16 * (2 - (np.asarray(f, dtype=float) / fp) ** 2))


def via(raio_mm, sub: Substrato, f):
    """(L [H], R [Ω]) de uma via de aterramento cilíndrica maciça (Goldfarb e Pucel, eqs. 3.61-3.62).

    L(pH) = 0,2·[h·ln((h + √(r² + h²))/r) + 1,5·(r - √(r² + h²))], h e r em µm;
    R = Rdc·√(1 + π f µ0 σ r²), Rdc = h/(σπr²).
    """
    h, r = sub.h_mm * 1e3, raio_mm * 1e3
    s = np.sqrt(r**2 + h**2)
    l_via = 0.2 * (h * np.log((h + s) / r) + 1.5 * (r - s)) * 1e-12
    rdc = sub.h_mm * 1e-3 / (sub.sigma * np.pi * (raio_mm * 1e-3) ** 2)
    r_via = rdc * np.sqrt(1 + np.pi * np.asarray(f, dtype=float) * 4e-7 * np.pi * sub.sigma * (raio_mm * 1e-3) ** 2)
    return float(l_via), r_via


def entalhe(w_mm, b_mm, sub: Substrato):
    """L [H] de um entalhe transversal de profundidade b numa linha W (eq. 3.59, 0 ≤ b/W ≤ 0,9):
    L/h (nH/mm) = 2·[1 - (Z0/Z0')·√(εeff/εeff')]², com Z0', εeff' da largura W - b."""
    if not 0 < b_mm < 0.9 * w_mm:
        raise ValueError("Profundidade do entalhe deve estar entre 0 e 0,9·W.")
    z, ee = analise(w_mm, sub)
    z2, ee2 = analise(w_mm - b_mm, sub)
    return float(2 * (1 - float(z) / float(z2) * np.sqrt(float(ee) / float(ee2))) ** 2 * sub.h_mm * 1e-9)


def limites_substrato(sub: Substrato, w_mm=None):
    """Frequências [Hz] a partir das quais surgem modos indesejados.

    - TE1 e TM1 de onda de superfície no dielétrico aterrado (Pozar, eq. 3.174 e Ex. 3.4);
      o TM0 não tem corte, mas só acopla fortemente à microstrip perto do TE1;
    - modo transversal da microstrip de largura W: fc ≈ c/(√εr·(2W + 0,8h)) (Pozar, seção 3.8).
    """
    h = sub.h_mm * 1e-3
    out = {}
    if sub.er > 1:
        out["acoplamento com TM0 (fT, Garg-Bahl eq. 2.140)"] = frequencia_ft(sub)
        out["TE1 (onda de superfície)"] = C0 / (4 * h * np.sqrt(sub.er - 1))
        out["TM1 (onda de superfície)"] = C0 / (2 * h * np.sqrt(sub.er - 1))
    if w_mm:
        out[f"modo transversal (W = {w_mm:.3g} mm)"] = C0 / (np.sqrt(sub.er) * (2 * w_mm * 1e-3 + 0.8 * h))
    return out


def texto_limites(sub: Substrato, w_mm=None, f_max=None):
    """Linhas de texto com os limites e um aviso se f_max [Hz] os ultrapassar."""
    lim = limites_substrato(sub, w_mm)
    linhas = [f"  {k}: {v / 1e9:.2f} GHz" for k, v in lim.items()]
    if f_max is not None:
        menor = min(lim.values()) if lim else np.inf
        if f_max > 0.8 * menor:
            linhas.append(f"  ⚠ f máx = {f_max / 1e9:g} GHz está perto/acima de {menor / 1e9:.2f} GHz: "
                          "modos superiores e radiação invalidam o modelo quase-TEM.")
    return linhas


# ---------------------------------------------------------------------------
# Fórmulas fechadas do Pozar (comparação)
# ---------------------------------------------------------------------------
def pozar_analise(w_mm, sub: Substrato):
    """(Z0, εe) pelas eqs. 3.195-3.196 do Pozar (precisão ~1%, t = 0)."""
    u = np.asarray(w_mm, dtype=float) / sub.h_mm
    ee = (sub.er + 1) / 2 + (sub.er - 1) / 2 / np.sqrt(1 + 12 / u)
    z = np.where(u <= 1, 60 / np.sqrt(ee) * np.log(8 / u + u / 4),
                 120 * np.pi / (np.sqrt(ee) * (u + 1.393 + 0.667 * np.log(u + 1.444))))
    return z, ee


def pozar_sintese(z0, sub: Substrato):
    """W [mm] pela eq. 3.197 do Pozar."""
    er = sub.er
    a = z0 / 60 * np.sqrt((er + 1) / 2) + (er - 1) / (er + 1) * (0.23 + 0.11 / er)
    u = 8 * np.exp(a) / (np.exp(2 * a) - 2)
    if u > 2:
        b = 377 * np.pi / (2 * z0 * np.sqrt(er))
        u = 2 / np.pi * (b - 1 - np.log(2 * b - 1) + (er - 1) / (2 * er) * (np.log(b - 1) + 0.39 - 0.61 / er))
    return float(u * sub.h_mm)


def pozar_alfa_condutor(w_mm, sub: Substrato, f):
    """αc ≈ Rs/(Z0·W) [Np/m] (eq. 3.199, aproximação para tiras largas)."""
    z0, _ = pozar_analise(w_mm, sub)
    return resistencia_superficial(np.asarray(f, dtype=float), sub) / (z0 * w_mm * 1e-3)
