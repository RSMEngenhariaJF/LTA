"""Aplicações de linhas EBG não uniformes: linha de atraso chirpada e ressonador de Bragg.

Referência: F. Martín (2015), seção 2.5.1 (eqs. 2.106 e 2.109-2.115).
Os perfis são de modulação da largura da trilha, analisados por CMT (seção 2.4).
"""
from dataclasses import dataclass

import numpy as np

from . import modos_acoplados as cm
from .microstrip import C0, Substrato, analise, parametros_linha, sintese

JANELAS_CDL = ("Gaussiana assimétrica (livro)", "Hamming", "Nenhuma")


def _eeff_ref(sub, z_ref):
    return float(analise(sintese(z_ref, sub), sub)[1])


# ---------------------------------------------------------------------------
# Linha de atraso chirpada (CDL), seção 2.5.1.5
# ---------------------------------------------------------------------------
@dataclass
class ProjetoCDL:
    f0: float          # Hz, frequência central (eq. 2.113)
    df: float          # Hz, banda de operação (eq. 2.112), com sinal de C
    psi: float         # s/Hz, inclinação do atraso de grupo (eq. 2.114)
    l0: float          # m, período no centro
    chirp: float       # 1/m², parâmetro C
    comprimento: float  # m, L = N·l0
    n_periodos: int
    eeff: float
    amplitude: float   # A

    def f_local(self, z):
        """Frequência refletida na posição z ∈ [-L/2, L/2] (eq. 2.111)."""
        return C0 / (2 * np.sqrt(self.eeff)) * (1 / self.l0 + self.chirp * np.asarray(z) / np.pi)


def projeto_cdl(f0, df, psi, sub: Substrato, z_ref=50.0, amplitude=0.4, n_periodos=None):
    """Síntese da CDL a partir de f0 [Hz], Δf [Hz] e ψ [s/Hz] (eq. 2.115).

    O sinal de ψ define o de C: ψ < 0 reflete primeiro as frequências altas (mais
    perdas), o que equaliza a perda de reflexão (livro, p. 101). L é arredondado para
    um número inteiro de períodos l0 (portas em Zref) e C é recalculado pela
    eq. 2.112 para manter Δf; ψ resultante muda na mesma proporção que L.
    """
    if f0 <= 0 or df <= 0 or psi == 0:
        raise ValueError("f0 e Δf devem ser positivos e ψ diferente de zero.")
    if df >= 2 * f0:
        raise ValueError("Δf deve ser menor que 2·f0 (a banda não pode incluir f ≤ 0).")
    ee = _eeff_ref(sub, z_ref)
    l0 = C0 / (2 * f0 * np.sqrt(ee))
    comp = C0 * abs(psi) * df / (2 * np.sqrt(ee))
    n = int(n_periodos) if n_periodos else max(1, int(round(comp / l0)))
    comp = n * l0
    chirp = np.sign(psi) * 2 * np.pi * np.sqrt(ee) * df / (C0 * comp)
    psi_real = 4 * np.pi * ee / (chirp * C0**2)
    return ProjetoCDL(f0, df, psi_real, l0, chirp, comp, n, ee, amplitude)


def janela_cdl(z, comp, tipo=JANELAS_CDL[0]):
    """W(z) para z ∈ [-L/2, L/2]. Livro (p. 101): exp(-4·((z - L/4)/L)²)."""
    z = np.asarray(z)
    if tipo == JANELAS_CDL[0]:
        return np.exp(-4 * ((z - comp / 4) / comp) ** 2)
    if tipo == JANELAS_CDL[1]:
        return 0.54 + 0.46 * np.cos(2 * np.pi * z / comp)
    return np.ones_like(z, dtype=float)


def perfil_cdl(p: ProjetoCDL, z_ref=50.0, janela=JANELAS_CDL[0], pts_por_periodo=32):
    """Z0(z) = Zref·exp(A·W(z)·sin(2πz/l0 + Cz² - CL²/4))  (eq. 2.109), z de 0 a L."""
    n_pts = pts_por_periodo * p.n_periodos + 1
    zc = np.linspace(-p.comprimento / 2, p.comprimento / 2, n_pts)
    fase = 2 * np.pi * zc / p.l0 + p.chirp * zc**2 - p.chirp * p.comprimento**2 / 4
    z0 = z_ref * np.exp(p.amplitude * janela_cdl(zc, p.comprimento, janela) * np.sin(fase))
    return cm.Perfil("CDL", zc + p.comprimento / 2, z0, p.l0, p.n_periodos)


def inclinacao_atraso(f, tau, f_ini, f_fim):
    """Ajuste linear de τ(f) na faixa: (ψ [s/Hz], τ no centro [s])."""
    sel = (f >= f_ini) & (f <= f_fim)
    if sel.sum() < 3:
        raise ValueError("Faixa insuficiente para ajustar o atraso de grupo.")
    a, b = np.polyfit(f[sel], tau[sel], 1)
    return a, a * 0.5 * (f_ini + f_fim) + b


# ---------------------------------------------------------------------------
# Ressonador de Bragg, seção 2.5.1.2
# ---------------------------------------------------------------------------
def comprimento_cavidade(phi1, phi2, beta, n=None):
    """Lm da eq. 2.106: φ1 + φ2 - 2βLm = 2πn. Retorna (Lm [m], n).

    phi1, phi2 em rad; beta em rad/m. Sem n, usa o que dá o menor Lm positivo.
    """
    soma = phi1 + phi2
    if n is None:
        n = int(np.floor(soma / (2 * np.pi)))
        if np.isclose((soma - 2 * np.pi * n) / (2 * beta), 0):
            n -= 1
    return (soma - 2 * np.pi * n) / (2 * beta), n


@dataclass
class ProjetoRessonador:
    f0: float
    refletor: cm.Perfil
    phi1: float        # rad, reflexão do refletor 1 vista da cavidade (S22)
    phi2: float        # rad, reflexão do refletor 2 vista da cavidade (S11)
    beta: float        # rad/m na cavidade em f0
    lm: float          # m
    n: int
    perfil: cm.Perfil  # estrutura completa


def projeto_ressonador(f0, k_amp, n_periodos, sub: Substrato, z_ref=50.0, perdas=False, dispersivo=False,
                       pts_por_periodo=64):
    """Dois refletores EBG senoidais (Bragg em f0) separados por uma cavidade de Lm (eq. 2.106).

    As fases dos refletores vêm do CMT em f0 com o mesmo modelo físico da análise.
    """
    ee = _eeff_ref(sub, z_ref)
    periodo = C0 / (2 * f0 * np.sqrt(ee))
    ref = cm.perfil("Senoidal", k_amp, periodo, n_periodos, z_ref, pts_por_periodo)
    f = np.array([f0])
    z0, eef, al = cm.parametros_linha(f, ref, sub, perdas, dispersivo)
    s11, _ = cm.resolver_cmt(f, ref.z, z0, eef, al)
    s22, _ = cm.resolver_cmt(f, ref.z, z0[:, ::-1], eef[:, ::-1], al[:, ::-1])
    w = sintese(z_ref, sub)
    _, ee_f, _ = parametros_linha(w, sub, f, False, dispersivo)
    beta = float(2 * np.pi * f0 * np.sqrt(ee_f[0, 0]) / C0)
    phi1, phi2 = float(np.angle(s22[0])), float(np.angle(s11[0]))
    lm, n = comprimento_cavidade(phi1, phi2, beta)
    # estrutura completa: refletor + cavidade uniforme + refletor
    n_cav = max(2, int(np.ceil(lm / periodo * pts_por_periodo)))
    z_cav = ref.z[-1] + np.linspace(0, lm, n_cav + 1)[1:]
    z = np.concatenate([ref.z, z_cav, z_cav[-1] + ref.z[1:]])
    z0 = np.concatenate([ref.z0, np.full(n_cav, z_ref), ref.z0[1:]])
    total = cm.Perfil("Ressonador de Bragg", z, z0, periodo, 2 * n_periodos)
    return ProjetoRessonador(f0, ref, phi1, phi2, beta, lm, n, total)


def fator_q(f, s21):
    """Pico de |S21|, QL = f0/BW-3dB (eq. 1.75) e Q0 = QL/(1 - 10^(-IL/20)) (eq. 1.74).

    Retorna (f_pico, IL [dB], QL, Q0); Q0 = inf se IL = 0 (sem perdas).
    """
    mag = np.abs(s21)
    i = int(np.argmax(mag))
    alvo = mag[i] / np.sqrt(2)
    a = i
    while a > 0 and mag[a] > alvo:
        a -= 1
    b = i
    while b < len(f) - 1 and mag[b] > alvo:
        b += 1
    if mag[a] > alvo or mag[b] > alvo:
        raise ValueError("Ressonância não resolvida na faixa: amplie a faixa ou aumente os pontos.")
    # interpolação linear das bordas de -3 dB
    fa = np.interp(alvo, [mag[a], mag[a + 1]], [f[a], f[a + 1]])
    fb = np.interp(alvo, [mag[b], mag[b - 1]], [f[b], f[b - 1]])
    il = -20 * np.log10(mag[i])
    ql = f[i] / (fb - fa)
    q0 = ql / (1 - 10 ** (-il / 20)) if il > 1e-9 else np.inf
    return f[i], il, ql, q0
