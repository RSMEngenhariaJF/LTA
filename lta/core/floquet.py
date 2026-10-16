"""Análise de Floquet: harmônicos espaciais e radiação leaky-wave.

Referência: F. Martín (2015), seção 2.2 (eqs. 2.1-2.14).

O método ABCD (eq. 2.21) dá βl na zona reduzida [0, π]. Aqui:
  - desdobra βl na zona estendida (fundamental contínuo, β crescente com f);
  - gera os harmônicos espaciais βn = β + 2πn/l (eq. 2.6);
  - identifica onda rápida |βn| < k0 (eq. 2.10), onde há radiação, e o ângulo
    do feixe cos θ = βn/k0 (eq. 2.11), medido a partir do eixo z da linha;
  - velocidades vp,n = ω/βn (eq. 2.8) e vg = dω/dβ (eq. 2.7, igual para todos).

O desdobramento supõe estrutura de mão direita (β cresce com f a partir de 0);
linhas de mão esquerda/CRLH (cap. 3) exigem outro critério.
"""
from dataclasses import dataclass

import numpy as np

from .microstrip import C0


def desdobrar(bl, tol=1e-3):
    """βl reduzido em [0, π] -> βl estendido contínuo e não decrescente, para f crescente.

    Os candidatos são 2πk ± βl; em cada frequência escolhe-se o mais próximo do
    valor anterior entre os que não o reduzem (tolerância tol para ruído numérico).
    Funciona com perdas (α > 0 em toda a faixa) e com bandas rejeitadas mais estreitas
    que o passo de frequência, que um critério baseado em α não detectaria.
    Nas bandas rejeitadas, βl fica preso na borda (0 ou π, a menos de perdas).
    """
    bl = np.asarray(bl, dtype=float)
    ext = np.empty_like(bl)
    ant = 0.0
    for i, b in enumerate(bl):
        k = np.floor(ant / (2 * np.pi))
        cand = np.array([2 * np.pi * (k + j) + s * b for j in (0, 1) for s in (1, -1)])
        cand = cand[cand >= ant - tol]
        ant = ext[i] = cand[np.argmin(cand - ant)] if cand.size else ant
    return ext


@dataclass
class Harmonico:
    n: int
    beta_l: np.ndarray        # βn·l [rad] (eq. 2.6)
    rapido: np.ndarray        # bool, |βn| < k0 (radia, eq. 2.10)
    theta_graus: np.ndarray   # ângulo do feixe (eq. 2.11), NaN onde não radia
    vp_c: np.ndarray          # vp,n / c (eq. 2.8)


def harmonicos(f, beta_l_ext, periodo, n_max=2):
    """Harmônicos espaciais n = -n_max..n_max a partir do fundamental estendido."""
    f = np.asarray(f, dtype=float)
    k0l = 2 * np.pi * f / C0 * periodo
    out = []
    for n in range(-n_max, n_max + 1):
        bn = beta_l_ext + 2 * np.pi * n
        rapido = np.abs(bn) < k0l
        with np.errstate(divide="ignore", invalid="ignore"):
            cos_t = np.where(k0l > 0, bn / k0l, np.nan)
            theta = np.where(rapido, np.degrees(np.arccos(np.clip(cos_t, -1, 1))), np.nan)
            vp = np.where(bn != 0, k0l / bn, np.inf)
        out.append(Harmonico(n, bn, rapido, theta, vp))
    return out


def velocidade_grupo(f, beta_l_ext, periodo):
    """vg/c = dω/dβ / c (eq. 2.7); NaN onde β não varia (bandas rejeitadas)."""
    w = 2 * np.pi * np.asarray(f, dtype=float)
    dbeta = np.gradient(beta_l_ext / periodo, w)
    with np.errstate(divide="ignore", invalid="ignore"):
        vg = np.where(np.abs(dbeta) > 1e-30, 1 / dbeta, np.nan)
    return vg / C0


def faixas_radiacao(f, harm: Harmonico):
    """Faixas contíguas em que o harmônico radia: [(f_ini, f_fim, θ_ini, θ_fim)]."""
    faixas = []
    r = harm.rapido
    i = 0
    while i < len(f):
        if r[i]:
            j = i
            while j + 1 < len(f) and r[j + 1]:
                j += 1
            faixas.append((f[i], f[j], harm.theta_graus[i], harm.theta_graus[j]))
            i = j + 1
        else:
            i += 1
    return faixas


def frequencia_broadside(periodo, eeff, n=-1):
    """Linha pouco perturbada: β_n = 0 em f = |n|·c/(l·√εeff) (feixe a 90°)."""
    return abs(n) * C0 / (periodo * np.sqrt(eeff))


def inicio_radiacao(periodo, eeff, n=-1):
    """Linha pouco perturbada: β_n = -k0 (backfire) em f = |n|·c/(l·(√εeff + 1))."""
    return abs(n) * C0 / (periodo * (np.sqrt(eeff) + 1))


# ---------------------------------------------------------------------------
# Antena leaky-wave (seção 4.2.4)
# ---------------------------------------------------------------------------
def celulas_lwa(alfa_k0, f, d, fracao=0.95):
    """Nº de células para irradiar `fracao` da potência: P_N = P0·e^(-2αNd) (eq. 4.51).

    alfa_k0: α/k0 (obtido por simulação EM de N células, como no livro); d: célula [m].
    Livro: α/k0 = 0,076 e 95% de potência irradiada levaram a N = 17.
    """
    if not 0 < fracao < 1 or alfa_k0 <= 0 or d <= 0:
        raise ValueError("Use 0 < fração < 1, α/k0 > 0 e célula > 0.")
    alfa = alfa_k0 * 2 * np.pi * f / C0
    return int(np.ceil(-np.log(1 - fracao) / (2 * alfa * d)))


def largura_feixe(comprimento, f, theta_graus):
    """Largura de feixe de meia potência [°] de uma fonte linear leaky-wave de comprimento L.

    Aproximação de abertura (Jackson e Oliner): Δθ ≈ λ0/(L·sen θ), θ a partir do eixo da
    linha (90° = broadside). Vale para α pequeno e feixe longe de endfire/backfire.
    """
    lam = C0 / np.asarray(f, dtype=float)
    s = np.sin(np.radians(np.asarray(theta_graus, dtype=float)))
    with np.errstate(divide="ignore"):
        return np.degrees(lam / (comprimento * s))
