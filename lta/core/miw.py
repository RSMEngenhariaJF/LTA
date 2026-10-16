"""Linhas de atraso magnetoindutivas (MIW) e eletroindutivas (EIW).

Referência: F. Martín (2015), seção 6.2 (eqs. 6.1-6.10).

MIW: cadeia de ressonadores (L, C) com indutância mútua M entre vizinhos (M < 0 no
arranjo coplanar, onda regressiva; M > 0 no axial). Análise de malhas:
    jωL + 1/(jωC) + 2jωM·cos βl = 0  ⇒  ω = ω0/√(1 + 2(M/L)·cos βl)   (eq. 6.1)
EIW (CSRRs com capacitância mútua CM): ω = ω0/√(1 - 2(CM/Cc)·cos βl) (eq. 6.4), dual.
"""
import numpy as np

from . import abcd
from .microstrip import C0


def fator_acoplamento(m_sobre_l=None, cm_sobre_cc=None):
    """κ da forma comum ω = ω0/√(1 + 2κ·cos βl): κ = M/L (MIW) ou κ = -CM/Cc (EIW)."""
    if m_sobre_l is not None:
        return m_sobre_l
    return -cm_sobre_cc


def frequencia(beta_l, f0, kappa):
    """f(βl) (eqs. 6.1 e 6.4)."""
    return f0 / np.sqrt(1 + 2 * kappa * np.cos(beta_l))


def bordas(f0, kappa):
    """Limites da banda: f0/√(1 ± 2|κ|) (eqs. 6.3 e 6.6). Retorna (f_inf, f_sup, FBW)."""
    if abs(kappa) >= 0.5:
        raise ValueError("|κ| deve ser < 0,5 (acoplamento fraco entre vizinhos).")
    a, b = f0 / np.sqrt(1 + 2 * abs(kappa)), f0 / np.sqrt(1 - 2 * abs(kappa))
    return a, b, (b - a) / f0


def velocidade_grupo(beta_l, f0, kappa, periodo):
    """vg = dω/dβ = l·κ·ω³·sen βl/ω0² (eq. 6.8); no centro (βl = π/2), vg = l·κ·ω0."""
    w0 = 2 * np.pi * f0
    w = 2 * np.pi * frequencia(beta_l, f0, kappa)
    return periodo * kappa * w**3 * np.sin(beta_l) / w0**2


def projeto_etiqueta(f0, kappa, periodo, largura_pulso, n_bits=None, n_ressonadores=None):
    """Etiqueta RFID por reflectometria (eqs. 6.7-6.10).

    Δτ = 2·n·l/vg(ω0) = 2n/(ω0·|κ|) (eq. 6.9) e nb = Δτ/τ (eq. 6.7). Informe n_bits
    (dá o nº de ressonadores) ou n_ressonadores (dá os bits).
    """
    w0 = 2 * np.pi * f0
    atraso_por_ressonador = 2 / (w0 * abs(kappa))
    if n_ressonadores is None:
        if not n_bits:
            raise ValueError("Informe o nº de bits ou de ressonadores.")
        n_ressonadores = int(np.ceil(n_bits * largura_pulso / atraso_por_ressonador))
    dt = n_ressonadores * atraso_por_ressonador
    return dict(n_ressonadores=n_ressonadores, atraso_total=dt, bits=dt / largura_pulso,
                comprimento=n_ressonadores * periodo, fator_lento=C0 / abs(periodo * kappa * w0),
                atraso_por_ressonador=atraso_por_ressonador)


# ---------------------------------------------------------------------------
# Circuito: célula T equivalente (Fig. 6.2b)
# ---------------------------------------------------------------------------
def abcd_miw(f, l_res, c_res, m):
    """Célula T: braços Z/2 + jωM e derivação -jωM, Z = jωL + 1/(jωC) (equivale à eq. 6.1)."""
    w = 2 * np.pi * np.asarray(f, dtype=float)
    z = 1j * w * l_res + 1 / (1j * w * c_res)
    return abcd.cascata(abcd.series(z / 2 + 1j * w * m), abcd.shunt(1 / (-1j * w * m)),
                        abcd.series(z / 2 + 1j * w * m))


def s11_etiqueta(f, l_res, c_res, m, n_por_bit, bits, refletor, z_ref):
    """S11 da linha de atraso com refletores ao fim de cada segmento de bit '1' e carga casada.

    refletor: susceptância normalizada b (admitância em derivação j·b/z_ref) que marca o bit.
    """
    cel = abcd_miw(f, l_res, c_res, m)
    seg = abcd.potencia(cel, n_por_bit)
    tot = np.broadcast_to(np.eye(2, dtype=complex), seg.shape).copy()
    y = np.full(np.shape(f), 1j * refletor / z_ref)
    for bit in bits:
        tot = tot @ seg
        if bit == "1":
            tot = tot @ abcd.shunt(y)
    return abcd.abcd_para_s(tot, z_ref)[:, 0, 0]


def resposta_pulso(f, s11, f0, largura):
    """Eco no tempo de um pulso gaussiano modulado em f0 com largura (FWHM do envelope) `largura`.

    f deve ser uniforme e começar perto de 0; retorna (t, envelope do eco).
    """
    sigma_t = largura / (2 * np.sqrt(2 * np.log(2)))
    espectro = np.exp(-((f - f0) ** 2) * (2 * np.pi * sigma_t) ** 2 / 2)
    df = f[1] - f[0]
    n = 1 << int(np.ceil(np.log2(8 * f.size)))
    x = np.zeros(n, dtype=complex)
    x[: f.size] = s11 * espectro
    sinal = np.fft.ifft(x) * n * df          # sinal analítico (só frequências positivas)
    t = np.arange(n) / (n * df)
    return t, np.abs(sinal) / np.abs(np.fft.ifft(np.pad(espectro, (0, n - f.size))) * n * df).max()
