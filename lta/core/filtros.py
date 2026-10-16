"""Filtros a partir do protótipo passa-baixas: passa-faixa, dual-band e células híbridas com CSRR.

Referência: F. Martín (2015), seção 4.2.3 (eqs. 4.38-4.50); protótipos de Matthaei,
Young e Jones (fórmulas de g para Butterworth e Chebyshev).
"""
import numpy as np

from . import abcd
from .multibanda import _raizes


# ---------------------------------------------------------------------------
# Protótipo passa-baixas
# ---------------------------------------------------------------------------
def prototipo(tipo, n, ripple_db=0.1):
    """Coeficientes g0..g(n+1) do protótipo passa-baixas (Ωc = 1 rad/s)."""
    if n < 1:
        raise ValueError("Ordem deve ser ≥ 1.")
    k = np.arange(1, n + 1)
    if tipo == "Butterworth":
        return np.concatenate([[1.0], 2 * np.sin((2 * k - 1) * np.pi / (2 * n)), [1.0]])
    if tipo != "Chebyshev":
        raise ValueError(f"Aproximação desconhecida: {tipo}")
    if ripple_db <= 0:
        raise ValueError("Ripple deve ser positivo.")
    beta = np.log(1 / np.tanh(ripple_db / 17.37))
    gama = np.sinh(beta / (2 * n))
    a = np.sin((2 * k - 1) * np.pi / (2 * n))
    b = gama**2 + np.sin(k * np.pi / n) ** 2
    g = [1.0, 2 * a[0] / gama]
    for i in range(1, n):
        g.append(4 * a[i - 1] * a[i] / (b[i - 1] * g[-1]))
    g.append(1.0 if n % 2 else 1 / np.tanh(beta / 4) ** 2)
    return np.array(g)


def g_periodico(n=3):
    """g único de um protótipo com n elementos iguais e corte de 3 dB em Ω = 1 (livro, p. 283: 1,521 para n = 3)."""
    def s21_db(g):
        m = np.eye(2, dtype=complex)[None]
        for i in range(n):
            el = abcd.shunt(np.array([1j * g])) if i % 2 == 0 else abcd.series(np.array([1j * g]))
            m = m @ el
        return 20 * np.log10(abs(abcd.abcd_para_s(m, 1.0)[0, 1, 0]))
    return _raizes(lambda g: s21_db(g) + 3.0103, np.linspace(0.3, 5, 400))[0]


# ---------------------------------------------------------------------------
# Passa-faixa com inversores de admitância e ressonadores em derivação (Fig. 4.28b)
# ---------------------------------------------------------------------------
def inversor(k):
    """Inversor ideal de impedância K (independente da frequência)."""
    return lambda f: np.broadcast_to(np.array([[0, 1j * k], [1j / k, 0]], dtype=complex),
                                     np.shape(f) + (2, 2))


def ressonadores(g, f0, fbw, z0=50.0):
    """(Ceq, Leq) dos tanques LC (eqs. 4.40-4.41) e banda de 3 dB de cada um (eq. 4.42)."""
    w0 = 2 * np.pi * f0
    ceq = g[1:-1] / (fbw * w0 * z0)
    leq = 1 / (w0**2 * ceq)
    delta = 2 * fbw * w0 / g[1:-1]
    return ceq, leq, delta


def resposta_inversores(f, g, y_ressonadores, z0=50.0):
    """S do filtro: inversores Z0 entre portas e ressonadores; o último ajusta g(n+1) (Chebyshev par)."""
    f = np.asarray(f, dtype=float)
    m = inversor(z0)(f)
    for y in y_ressonadores:
        m = m @ abcd.shunt(y(f)) @ inversor(z0)(f)
    m = m.copy()
    if not np.isclose(g[-1], 1.0):
        # carga g(n+1) ≠ 1 (Chebyshev de ordem par): o último inversor passa a K = Z0·√g(n+1)
        # (dois inversores Z0 em cascata se cancelam)
        m = m @ inversor(z0)(f) @ inversor(z0 * np.sqrt(g[-1]))(f)
    return abcd.abcd_para_s(m, z0)


def tanque(c, l):
    return lambda f: 1j * 2 * np.pi * np.asarray(f) * c + 1 / (1j * 2 * np.pi * np.asarray(f) * l)


# ---------------------------------------------------------------------------
# Transformação dual-band (eqs. 4.49-4.50)
# ---------------------------------------------------------------------------
def mapa_dual(f, f1, f2):
    """x(ω) = (ω/ω0 - ω0/ω)/w0, ω0 = √(ω1ω2), w0 = (ω2 - ω1)/ω0 (eqs. 4.50): x(ω1) = -1, x(ω2) = +1."""
    w0n = (f2 - f1) / np.sqrt(f1 * f2)
    f0 = np.sqrt(f1 * f2)
    return (np.asarray(f) / f0 - f0 / np.asarray(f)) / w0n


def ressonador_dual(ceq, leq, f1, f2):
    """Tanque (Ceq ‖ Leq) centrado em ω0 = √(ω1ω2) → ressonador dual-band (eq. 4.50a em cada elemento).

    Ceq vira C' ‖ L'; Leq vira L'' + C'' em série: resultado (L' ‖ C') ‖ (L'' + C''), a mesma
    derivação da E-CRLH (Fig. 3.25c). Retorna (Cp, Lp, Ls, Cs).
    """
    w0 = 2 * np.pi * np.sqrt(f1 * f2)
    w0n = (f2 - f1) / np.sqrt(f1 * f2)
    cp = ceq / w0n
    lp = w0n / (w0**2 * ceq)
    ls = leq / w0n
    cs = w0n / (w0**2 * leq)
    return cp, lp, ls, cs


def admitancia_dual(cp, lp, ls, cs):
    def y(f):
        w = 2 * np.pi * np.asarray(f, dtype=float)
        return 1j * w * cp + 1 / (1j * w * lp) + 1 / (1j * w * ls + 1 / (1j * w * cs))
    return y


# ---------------------------------------------------------------------------
# Célula híbrida com CSRR (Fig. 3.60b), eqs. 4.38-4.46
# ---------------------------------------------------------------------------
def celula_hibrida(f0, fbw, g, fz, z0=50.0):
    """Elementos da célula: Cg (eq. 4.46) e Ld, C, Lc, Cc (eqs. 3.86 e 4.43-4.45).

    Derivação: Ld ‖ [C em série com (Lc ‖ Cc)]. Condições: Zp(ω0) = jZ0, Zp(ω1) = jZ0/2,
    Zp(ω2) = ∞, ω1,2 = ω0 ∓ Δ/2 com Δ = 2·FBW·ω0/g (eq. 4.42), e zero em fz (eq. 3.86).
    A indutância série da linha é desprezada, como no livro.
    """
    w0 = 2 * np.pi * f0
    d = 2 * fbw * w0 / g
    w1, w2 = w0 - d / 2, w0 + d / 2
    y = 1 / (2 * np.pi * fz) ** 2                     # Lc·(C + Cc)
    if fz <= f0:
        raise ValueError("Use fz > f0 (zero de transmissão acima da banda; o livro sugere fz = 2f0).")

    def partes(x):                                    # x = Lc·Cc
        ldc = (1 - w2**2 * y) / (w2**2 * (1 - w2**2 * x))            # Yp(ω2) = 0
        den0 = 1 - w0**2 * y - w0**2 * ldc * (1 - w0**2 * x)
        return ldc, z0 * den0 / (w0 * (1 - w0**2 * y))               # X(ω0) = Z0 → Ld

    def resid(x):
        with np.errstate(all="ignore"):
            ldc, ld = partes(x)
            den = 1 - w1**2 * y - w1**2 * ldc * (1 - w1**2 * x)
            return w1 * ld * (1 - w1**2 * y) / den / (z0 / 2) - 1       # X(ω1) = Z0/2

    sols = []
    for x in _raizes(resid, np.logspace(np.log10(y) - 6, np.log10(y), 6000)):
        ldc, ld = partes(x)
        c = ldc / ld
        lc = (y - x) / c
        if min(ld, c, lc, x) > 0:
            sols.append(dict(Cg=1 / (2 * z0 * w0), Ld=ld, C=c, Lc=lc, Cc=x / lc))
    if not sols:
        raise ValueError("Sem solução física para a célula híbrida: altere fz ou a banda.")
    return sols[0]


def abcd_hibrida(f, p):
    """Célula T da Fig. 3.60b sem a indutância série: braços 2Cg; derivação Ld ‖ [C + (Lc ‖ Cc)]."""
    w = 2 * np.pi * np.asarray(f, dtype=float)
    zb = 1 / (1j * w * p["C"]) + 1j * w * p["Lc"] / (1 - w**2 * p["Lc"] * p["Cc"])
    yp = 1 / (1j * w * p["Ld"]) + 1 / zb
    braco = 1 / (2j * w * p["Cg"])
    return abcd.cascata(abcd.series(braco), abcd.shunt(yp), abcd.series(braco))


def banda_3db(f, s21):
    """(f_ini, f_fim, f_centro geométrico) da banda de 3 dB em torno do máximo de |S21|."""
    y = np.nan_to_num(abcd.db(s21), nan=-np.inf)        # pontos singulares do modelo (ex.: em fz)
    i = int(np.argmax(y))
    ok = y > y[i] - 3.0103
    a, b = i, i
    while a > 0 and ok[a - 1]:
        a -= 1
    while b < len(f) - 1 and ok[b + 1]:
        b += 1
    return f[a], f[b], np.sqrt(f[a] * f[b])
