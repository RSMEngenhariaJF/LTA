"""Linhas de transmissão não lineares (NLTL): sólitons, decomposição de pulsos e harmônicos.

Referência: F. Martín (2015), seção 5.4 (eqs. 5.13-5.25).

Modelo (Fig. 5.39b): seções com indutância série L e capacitância C(V) = Co + CD(V).
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class LinhaNL:
    L: float           # H por seção
    Co: float          # F por seção (linha)
    capacitancia: callable   # CD(V) [F] do dispositivo não linear

    def C(self, v):
        return self.Co + self.capacitancia(v)

    @property
    def z0_linha(self):
        """√(L/Co): impedância da linha sem os dispositivos (referência de projeto)."""
        return np.sqrt(self.L / self.Co)

    def z_pequeno_sinal(self, v=0.0):
        return np.sqrt(self.L / self.C(v))

    def f_bragg(self, v=0.0):
        """Corte da rede LC (pequeno sinal): 1/(π√(L·C(V)))."""
        return 1 / (np.pi * np.sqrt(self.L * self.C(v)))


# ---------------------------------------------------------------------------
# Sólitons (eqs. 5.20-5.25)
# ---------------------------------------------------------------------------
def _grade(v_max, n=4001):
    return np.linspace(0.0, v_max, n)


def funcao_F(linha: LinhaNL, v):
    """F(V) = ∫0^V C(V') dV' (eq. 5.21), por trapézios em uma grade de 0 a V."""
    vs = _grade(float(v))
    c = linha.C(vs)
    return np.trapezoid(c, vs)


def funcoes_FG(linha: LinhaNL, v_max, n=4001):
    """Grade V e as integrais F(V) (5.21) e G(V) = ∫0^V F (5.25)."""
    v = _grade(v_max, n)
    c = linha.C(v)
    dv = np.diff(v)
    f = np.concatenate([[0.0], np.cumsum(0.5 * (c[1:] + c[:-1]) * dv)])
    g = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * dv)])
    return v, f, g


def faixa_atraso(linha: LinhaNL, v_inf=1e3):
    """Atrasos por seção possíveis para sólitons: de √(L·C(∞)) a √(L·C(0)) (livro, p. 411)."""
    return np.sqrt(linha.L * linha.C(v_inf)), np.sqrt(linha.L * linha.C(0.0))


def atraso_soliton(linha: LinhaNL, v_max):
    """Atraso por seção do sóliton de amplitude Vmax: T = √(2L·G(Vmax))/Vmax (eq. 5.24)."""
    v, f, g = funcoes_FG(linha, v_max)
    return np.sqrt(2 * linha.L * g[-1]) / v_max


def amplitude_soliton(linha: LinhaNL, t, v_lim=200.0):
    """Inversa da eq. 5.24 (bisseção): Vmax do sóliton com atraso T por seção."""
    tmin, tmax = faixa_atraso(linha)
    if not tmin < t < tmax:
        raise ValueError(f"T fora da faixa de sólitons ({tmin * 1e12:.3f} a {tmax * 1e12:.3f} ps).")
    a, b = 1e-6, v_lim
    if atraso_soliton(linha, b) > t:
        raise ValueError("Amplitude necessária acima do limite de busca.")
    for _ in range(80):
        m = 0.5 * (a + b)
        if atraso_soliton(linha, m) > t:     # amplitudes maiores são mais rápidas (T menor)
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def perfil_soliton(linha: LinhaNL, t, n=2001):
    """Forma V(s) do sóliton com atraso T por seção, por quadratura da eq. 5.20.

    Multiplicando 5.20 por dV/ds e integrando: (T⁴/24L)·(dV/ds)² = G(V) - T²V²/(2L), logo
    s(V) = ∫V^Vmax dV / √((24L/T⁴)(G - T²V²/2L)). Retorna (s [s], V [V]) simétrico em s.
    """
    vmax = amplitude_soliton(linha, t)
    v, _, g = funcoes_FG(linha, vmax, 20001)
    h = g - t**2 * v**2 / (2 * linha.L)
    h = np.maximum(h, 1e-300)
    dsdv = 1 / np.sqrt(24 * linha.L / t**4 * h)
    # s(V) de Vmax para baixo; a singularidade integrável em Vmax é tratada pelo passo fino
    s = np.concatenate([np.cumsum((0.5 * (dsdv[1:] + dsdv[:-1]) * np.diff(v))[::-1])[::-1], [0.0]])
    sel = (v > vmax * 1e-3) & np.isfinite(s)
    v, s = v[sel], s[sel]
    idx = np.linspace(0, v.size - 1, min(n, v.size)).astype(int)
    v, s = v[idx], s[idx]
    return np.concatenate([-s, s[::-1]]), np.concatenate([v, v[::-1]])


# ---------------------------------------------------------------------------
# Simulação no tempo (eqs. 5.13-5.15), várias linhas de uma vez
# ---------------------------------------------------------------------------
def simular(linha: LinhaNL, n_secoes, fonte, t_final, dt, r_fonte=None, r_carga=None, nos=None):
    """Integra (RK4) as linhas de n_secoes (lista) com a mesma fonte.

    Nó i: C(Vi)·dVi/dt = I(i) - I(i+1); indutor I(i): L·dI/dt = V(i-1) - V(i), com V(-1) a
    tensão da fonte vs(t) atrás de r_fonte (0 = fonte ideal). Carga r_carga no último nó.
    Retorna (t, saídas (n_linhas, n_t), nós (n_t, len(nos)) da 1ª linha).
    """
    ns = np.atleast_1d(np.asarray(n_secoes, dtype=int))
    nmax = int(ns.max())
    rf = linha.z0_linha if r_fonte is None else r_fonte
    rc = linha.z0_linha if r_carga is None else r_carga
    b = ns.size
    ativo = np.arange(nmax)[None, :] < ns[:, None]           # (b, nmax)
    ultimo = np.zeros((b, nmax), dtype=bool)
    ultimo[np.arange(b), ns - 1] = True
    v = np.zeros((b, nmax))
    i = np.zeros((b, nmax))                                   # i[:, k]: indutor que chega ao nó k
    passos = int(np.ceil(t_final / dt))
    t = np.arange(passos) * dt
    saida = np.empty((b, passos))
    reg = None if nos is None else np.empty((passos, len(nos)))
    vfonte = np.zeros(b)

    def deriv(tt, v, i):
        vs = fonte(tt)
        if rf > 0:
            v_in = vs - rf * i[:, 0]                          # tensão após a resistência da fonte
        else:
            v_in = np.full(b, vs)
        vprev = np.concatenate([v_in[:, None], v[:, :-1]], axis=1)
        di = np.where(ativo, (vprev - v) / linha.L, 0.0)
        iprox = np.concatenate([i[:, 1:], np.zeros((b, 1))], axis=1)
        iprox = np.where(ultimo, v / rc, iprox)
        dv = np.where(ativo, (i - iprox) / linha.C(v), 0.0)
        return dv, di

    for n, tt in enumerate(t):
        a1, b1 = deriv(tt, v, i)
        a2, b2 = deriv(tt + dt / 2, v + dt / 2 * a1, i + dt / 2 * b1)
        a3, b3 = deriv(tt + dt / 2, v + dt / 2 * a2, i + dt / 2 * b2)
        a4, b4 = deriv(tt + dt, v + dt * a3, i + dt * b3)
        v = v + dt / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        i = i + dt / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
        saida[:, n] = v[np.arange(b), ns - 1]
        if reg is not None:
            reg[n] = v[0, nos]
    return t + dt, saida, reg


def eficiencia_harmonica(t, v_saida, f0, ordem, r_carga, p_ref, n_periodos=8):
    """P(ordem·f0) na carga / p_ref, pelos últimos n_periodos (regime)."""
    sel = t >= t[-1] - n_periodos / f0
    tt, vv = t[sel], np.atleast_2d(v_saida)[:, sel]
    amp = 2 * np.abs(np.mean(vv * np.exp(-2j * np.pi * ordem * f0 * tt), axis=1))
    return amp**2 / (2 * r_carga) / p_ref
