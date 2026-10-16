"""Projeto de linhas artificiais dual-band e multibanda (impedância e dispersão).

Referência: F. Martín (2015), seção 4.2.2 (eqs. 4.20-4.37).

Condição geral de uma célula T sem perdas com braços Za = jXa e derivação Zp = jXp
para fase φ (comprimento elétrico, < 0 na banda LH) e impedância de Bloch Z
(inversão das eqs. 2.30 e 2.33):
    Xa = Z·tan(φ/2),   Xp = -Z/sen(φ)
Para φ = -90° / +90° isso dá as eqs. 4.28 (Xa = ∓Z, Xp = ±Z).
"""
from dataclasses import dataclass

import numpy as np

from . import abcd
from .crlh import CelulaCRLH
from .csrr import CircuitoCSRR


def condicoes(fase, z):
    """(Xa, Xp) que dão fase [rad] e impedância z na célula T."""
    if np.isclose(np.sin(fase), 0):
        raise ValueError("Fases múltiplas de 180° não definem uma célula T única.")
    return z * np.tan(fase / 2), -z / np.sin(fase)


def _raizes(g, xs, tol_rel=1e-6):
    """Raízes de g por bisseção nas trocas de sinal de uma varredura, descartando polos."""
    gs = np.array([g(x) for x in xs])
    out = []
    for i in range(len(xs) - 1):
        if not (np.isfinite(gs[i]) and np.isfinite(gs[i + 1])) or np.sign(gs[i]) == np.sign(gs[i + 1]):
            continue
        a, b = xs[i], xs[i + 1]
        ga = gs[i]
        for _ in range(100):
            m = 0.5 * (a + b)
            gm = g(m)
            if np.sign(gm) == np.sign(ga):
                a, ga = m, gm
            else:
                b = m
        escala = max(abs(gs[i]), abs(gs[i + 1]), 1e-300)
        if abs(g(0.5 * (a + b))) < tol_rel * escala * 1e3 or abs(g(0.5 * (a + b))) < 1e-9 * escala:
            out.append(0.5 * (a + b))
    return out


# ---------------------------------------------------------------------------
# Linha convencional + linha PLH (eqs. 4.20-4.27)
# ---------------------------------------------------------------------------
@dataclass
class LinhaRHPLH:
    atraso: float     # l/vp [s] da linha convencional (eq. 4.21)
    CL: float         # F (eq. 4.22)
    LL: float         # H (eq. 4.23)
    n: int
    z0: float

    def fase(self, f):
        """φ = ω·l/vp - N/(ω√(LL·CL)) (eq. 4.20, PLH homogênea)."""
        w = 2 * np.pi * np.asarray(f, dtype=float)
        return w * self.atraso - self.n / (w * np.sqrt(self.LL * self.CL))


def rh_plh(f1, fase1, f2, fase2, z0, n=1):
    """Linha convencional em cascata com N células PLH (LL, CL), eqs. 4.21-4.23 e 4.24."""
    if not 0 < f1 < f2:
        raise ValueError("É preciso 0 < f1 < f2.")
    w1, w2 = 2 * np.pi * f1, 2 * np.pi * f2
    cond = fase2 * w1 - fase1 * w2
    if cond <= 0:
        raise ValueError("Condição 4.24 violada (φ2·ω1 - φ1·ω2 ≤ 0): LL e CL seriam negativos. "
                         f"Com fases positivas, exige φ2/φ1 > f2/f1 = {f2 / f1:.3f}.")
    atraso = (fase2 * w2 - fase1 * w1) / (w2**2 - w1**2)
    b = w1 * w2 * cond / (w2**2 - w1**2)               # N/√(LL·CL)
    return LinhaRHPLH(atraso, n / (z0 * b), z0 * n / b, n, z0)


# ---------------------------------------------------------------------------
# CRLH canônica (eqs. 4.28-4.32), fases e impedâncias quaisquer
# ---------------------------------------------------------------------------
def crlh_dual(f1, fase1, z1, f2, fase2, z2):
    """Célula CRLH (Fig. 3.20) com fase/impedância especificadas em f1 e f2.

    Braço: X = ωLR/2 - 1/(2ωCL); derivação: B = ωCR - 1/(ωLL) = -1/Xp.
    Com φ1 = -90°, φ2 = +90° reproduz as eqs. 4.29; com z1 = z2, a célula é
    balanceada e ωs = ωp = √(ω1ω2) (eqs. 4.31-4.32).
    """
    w = 2 * np.pi * np.array([f1, f2])
    xa, xp = zip(condicoes(fase1, z1), condicoes(fase2, z2))
    lr2, icl2 = np.linalg.solve(np.column_stack([w, -1 / w]), xa)      # LR/2 e 1/(2CL)
    cr, ill = np.linalg.solve(np.column_stack([w, -1 / w]), -1 / np.array(xp))
    if min(lr2, icl2, cr, ill) <= 0:
        raise ValueError("Especificação sem solução com elementos positivos (fases devem crescer com f).")
    return CelulaCRLH(2 * lr2, cr, 1 / ill, 1 / (2 * icl2))


# ---------------------------------------------------------------------------
# CSRR + gap (eqs. 4.33-4.34), com zero de transmissão
# ---------------------------------------------------------------------------
def csrr_dual(f1, fase1, z1, f2, fase2, z2, fz):
    """Célula CSRR + gap (Fig. 3.38b) com fase/impedância em f1, f2 e zero de transmissão fz.

    Ramo série (L, Cg) sai direto (eqs. 4.33a-b). Ramo em derivação: C, Lc, Cc com
    Xp(ω1), Xp(ω2) e ωz² = 1/(Lc(C + Cc)) (eqs. 4.34 generalizadas, como na nota 14
    do livro). Para cada x = Lc·Cc as duas primeiras são lineares em (1/C, Lc); a
    terceira dá x por busca de raiz. Retorna a lista de soluções físicas.
    """
    if not 0 < fz < f1 < f2:
        raise ValueError("É preciso 0 < fz < f1 < f2 (zero de transmissão abaixo da banda LH).")
    w1, w2, wz = 2 * np.pi * np.array([f1, f2, fz])
    (xa1, xp1), (xa2, xp2) = condicoes(fase1, z1), condicoes(fase2, z2)
    l2, icg2 = np.linalg.solve([[w1, -1 / w1], [w2, -1 / w2]], [xa1, xa2])
    if l2 <= 0 or icg2 <= 0:
        raise ValueError("Ramo série sem solução positiva.")

    def lineares(x):
        m = [[-1 / w1, w1 / (1 - w1**2 * x)], [-1 / w2, w2 / (1 - w2**2 * x)]]
        return np.linalg.solve(m, [xp1, xp2])         # (1/C, Lc)

    def g(x):
        with np.errstate(all="ignore"):
            u, lc = lineares(x)
            return (lc / u + x) * wz**2 - 1

    xs = np.logspace(np.log10(1 / (w2 * 50) ** 2), np.log10(1 / (wz / 50) ** 2), 3000)
    sols = []
    for x in _raizes(g, xs):
        u, lc = lineares(x)
        if u > 0 and lc > 0 and x > 0:
            sols.append(CircuitoCSRR(2 * l2, 1 / (2 * icg2), 1 / u, lc, x / lc))
    if not sols:
        raise ValueError("Sem solução física para o ramo em derivação (C, Lc, Cc > 0). Com Z1 = Z2 e "
                         "φ = ∓90°, C → ∞ (livro, p. 256): use trechos de linha ou outra fase.")
    return sols


def compensar_linhas(fase1, fase2, f1, f2, fase_linha_f1):
    """Fases que a célula deve ter quando há dois trechos de linha (Z = Zo) nas pontas.

    φcélula = φ - 2φL(f), com φL proporcional a f (linha não dispersiva). Livro, p. 259:
    φL(f1) = 15° leva φ1 = -90° → -120° e φ2 = +90° → +30°.
    """
    fl1 = np.radians(fase_linha_f1)
    return fase1 - 2 * fl1, fase2 - 2 * fl1 * f2 / f1


# ---------------------------------------------------------------------------
# CRLH estendida de ordem 4 (E-CRLH, Fig. 3.25c), quad-band (eq. 4.37)
# ---------------------------------------------------------------------------
@dataclass
class CelulaECRLH:
    """Braço série: Lhs/2 + 2Chs em série com (2Lhp ‖ Chp/2).
    Derivação: (Lvp ‖ Cvp) em paralelo com (Lvs + Cvs). Topologia conferida com os
    valores publicados (livro, p. 267)."""
    Lhs: float
    Chs: float
    Lhp: float
    Chp: float
    Lvs: float
    Cvs: float
    Lvp: float
    Cvp: float

    def ramos(self, f):
        w = 2 * np.pi * np.asarray(f, dtype=float)
        xa = w * self.Lhs / 2 - 1 / (w * 2 * self.Chs) + w * 2 * self.Lhp / (1 - w**2 * self.Lhp * self.Chp)
        b = w * self.Cvp - 1 / (w * self.Lvp) + w * self.Cvs / (1 - w**2 * self.Lvs * self.Cvs)
        return xa, b

    def abcd(self, f):
        xa, b = self.ramos(f)
        return abcd.cascata(abcd.series(1j * xa), abcd.shunt(1j * b), abcd.series(1j * xa))

    def params(self):
        return {k: getattr(self, k) for k in ("Lhs", "Chs", "Lhp", "Chp", "Lvs", "Cvs", "Lvp", "Cvp")}


def _foster4(w, alvo, sinal):
    """Resolve y·(a·ω - c/ω) + d·ω/(1 - ω²·r) = alvo nas 4 frequências (3 lineares + r).

    sinal=+1: reatância do braço (a = Lhs/2, c = 1/(2Chs), d = 2Lhp, r = Lhp·Chp);
    mesma forma para a susceptância da derivação (a = Cvp, c = 1/Lvp, d = Cvs, r = Lvs·Cvs).
    """
    def lineares(r):
        m = np.column_stack([w, -1 / w, w / (1 - w**2 * r)])
        return np.linalg.solve(m[:3], alvo[:3])

    def g(r):
        with np.errstate(all="ignore"):
            a, c, d = lineares(r)
            return (a * w[3] - c / w[3] + d * w[3] / (1 - w[3] ** 2 * r)) / abs(alvo[3]) - np.sign(alvo[3])

    rs = np.logspace(np.log10(1 / (w[-1] * 20) ** 2), np.log10(1 / (w[0] / 20) ** 2), 4000)
    sols = []
    for r in _raizes(g, rs):
        a, c, d = lineares(r)
        if min(a, c, d) > 0:
            sols.append((a, c, d, r))
    return sols


def ecrlh_quad(freqs, z):
    """E-CRLH com φ = -90° nas frequências ímpares e +90° nas pares (eqs. 4.37)."""
    f = np.sort(np.asarray(freqs, dtype=float))
    if f.size != 4 or f[0] <= 0:
        raise ValueError("Informe 4 frequências positivas.")
    w = 2 * np.pi * f
    n = np.arange(1, 5)
    xa = (-1.0) ** n * z                  # eq. 4.37a
    bp = -1 / ((-1.0) ** (n + 1) * z)     # B = -1/Xp, eq. 4.37b
    serie, deriv = _foster4(w, xa, 1), _foster4(w, bp, 1)
    if not serie or not deriv:
        raise ValueError("Sem solução E-CRLH com elementos positivos para essas frequências.")
    out = []
    for a, c, d, r in serie:
        for ap, cp, dp, rp in deriv:
            lhp = d / 2
            out.append(CelulaECRLH(2 * a, 1 / (2 * c), lhp, r / lhp, rp / dp, dp, 1 / cp, ap))
    return out
