"""Linha microstrip com CSRR no terra e capacitância série (gap ou interdigital).

Referência: F. Martín (2015), seções 3.3.2.1 e 3.5.2.2 (eqs. 3.86-3.92).

Projeto em três níveis, do mais confiável ao menos confiável:
  1. Circuito (Fig. 3.38b): elementos para transição balanceada f0, ZB(f0) e zero
     de transmissão fz (eqs. 3.86, 3.91, 3.92). Exato para o modelo.
  2. Modelo melhorado (Fig. 3.40, eqs. 3.87-3.88): separa a capacitância série
     física Cs do gap da capacitância paralela Cpar = Cf + CL.
  3. Geometria inicial: raio do CSRR por estimativa quase-estática da ressonância
     (dualidade com o SRR, seção 3.3.2.1) e capacitor interdigital pela fórmula de
     Bahl. Erro típico de ±20%: o livro fecha o projeto com simulação EM e
     extração de parâmetros (Apêndice G) ou space mapping (Apêndice H).
"""
from dataclasses import dataclass

import numpy as np

from .extracao import abcd_csrr
from .microstrip import C0, Substrato, analise, sintese

MU0 = 4e-7 * np.pi
EPS0 = 8.8541878128e-12


def elipk(k):
    """Integral elíptica completa de 1ª espécie K(k) pela média aritmético-geométrica."""
    a, b = 1.0, np.sqrt(1 - k * k)
    for _ in range(40):
        a, b = (a + b) / 2, np.sqrt(a * b)
    return np.pi / (2 * a)


# ---------------------------------------------------------------------------
# Circuito
# ---------------------------------------------------------------------------
@dataclass
class CircuitoCSRR:
    L: float
    Cg: float
    C: float
    Lc: float
    Cc: float

    def params(self):
        return dict(L=self.L, Cg=self.Cg, C=self.C, Lc=self.Lc, Cc=self.Cc)

    def abcd(self, f):
        return abcd_csrr(f, self.params())

    @property
    def fz(self):
        """Zero de transmissão (eq. 3.86)."""
        return 1 / (2 * np.pi * np.sqrt(self.Lc * (self.C + self.Cc)))

    @property
    def fs(self):
        return 1 / (2 * np.pi * np.sqrt(self.L * self.Cg))

    @property
    def fp(self):
        """Ressonância do CSRR (ramo em derivação aberto)."""
        return 1 / (2 * np.pi * np.sqrt(self.Lc * self.Cc))

    @property
    def balanceada(self):
        return np.isclose(self.fs, self.fp, rtol=1e-3)

    def zb(self, f):
        """ZB pela eq. 3.90 (célula T da Fig. 3.38b)."""
        w = 2 * np.pi * np.asarray(f, dtype=float)
        ws, wp = 2 * np.pi * self.fs, 2 * np.pi * self.fp
        a = 1 - ws**2 / w**2
        b = 1 - wp**2 / w**2
        z2 = self.L / self.Cc * a / b - self.L**2 * w**2 / 4 * a**2 + self.L / self.C * a
        return np.sqrt(z2 + 0j)

    def modelo_melhorado(self):
        """(Cs, Cpar) da Fig. 3.40 pela inversão das eqs. 3.87-3.88."""
        cpar = self.Cg * self.C / (self.C + 2 * self.Cg)
        return self.Cg**2 / (self.C + 2 * self.Cg), cpar


def projeto_balanceado(f0, zb, fz, L):
    """Elementos da linha balanceada (eq. 3.91) com ZB(f0) = √(L/Cc) (eq. 3.92) e zero em fz (eq. 3.86).

    L é a indutância da linha na célula (escolhida pelo comprimento da célula).
    """
    if not 0 < fz < f0:
        raise ValueError("O zero de transmissão fz deve ficar abaixo da transição f0.")
    if L <= 0 or zb <= 0:
        raise ValueError("L e ZB devem ser positivos.")
    w0, wz = 2 * np.pi * f0, 2 * np.pi * fz
    cc = L / zb**2
    lc = 1 / (w0**2 * cc)
    c = 1 / (wz**2 * lc) - cc
    if c <= 0:
        raise ValueError("fz alto demais para esse L e ZB (C ≤ 0): reduza fz.")
    return CircuitoCSRR(L, 1 / (w0**2 * L), c, lc, cc)


def indutancia_linha(z_linha, sub: Substrato, comprimento_mm):
    """L ≈ Z0·√εeff·ℓ/c de um trecho de linha (aproximação de linha curta)."""
    w = sintese(z_linha, sub)
    _, ee = analise(w, sub)
    return z_linha * np.sqrt(ee) * comprimento_mm * 1e-3 / C0, float(w)


# ---------------------------------------------------------------------------
# Geometria
# ---------------------------------------------------------------------------
@dataclass
class GeometriaCSRR:
    r_ext: float    # mm, raio externo do anel externo
    c: float        # mm, largura das fendas
    d: float        # mm, distância entre as fendas (metal)
    abertura: float = 0.3   # mm, corte em cada anel

    @property
    def r0(self):
        """Raio médio (centro da fita metálica entre as fendas)."""
        return self.r_ext - self.c - self.d / 2

    def valida(self):
        if min(self.c, self.d, self.abertura) <= 0:
            raise ValueError("c, d e a abertura devem ser positivos.")
        if self.r_ext - 2 * self.c - self.d <= 0.1:
            raise ValueError("Raio externo pequeno demais para c e d (disco interno inexistente).")


def estimativa_csrr(geo: GeometriaCSRR, sub: Substrato):
    """(Lc, Cc, fo) aproximados do CSRR (seção 3.3.2.1).

    Lc = Lo/4, Lo = 2π·r0·L'CPW (fita d entre fendas c). Cc pela dualidade com o
    SRR: Cc = 4·ε·Ls/μ0, ε = ε0(εr + 1)/2, Ls = μ0·r0·(ln(8r0/c) - 1/2).
    Contra o livro: fo -16% (Fig. 3.39) e dentro da banda (Fig. 3.37). Os
    elementos isolados erram mais; use só fo, com calibração por EM/medição.
    """
    geo.valida()
    r0 = geo.r0 * 1e-3
    k = geo.d / (geo.d + 2 * geo.c)
    lpul = MU0 / 4 * elipk(np.sqrt(1 - k * k)) / elipk(k)
    lc = 2 * np.pi * r0 * lpul / 4
    ls = MU0 * r0 * (np.log(8 * geo.r0 / geo.c) - 0.5)
    cc = 4 * EPS0 * (sub.er + 1) / 2 * ls / MU0
    return lc, cc, 1 / (2 * np.pi * np.sqrt(lc * cc))


def raio_para_fo(fo, c, d, sub: Substrato, abertura=0.3, calibracao=1.0):
    """r_ext [mm] cuja ressonância estimada (× calibração) é fo [Hz] (bisseção)."""
    alvo = fo / calibracao
    lo, hi = 2 * c + d + 0.11, 60.0
    if estimativa_csrr(GeometriaCSRR(lo, c, d, abertura), sub)[2] < alvo:
        raise ValueError("fo alta demais para essas fendas: reduza c ou d.")
    if estimativa_csrr(GeometriaCSRR(hi, c, d, abertura), sub)[2] > alvo:
        raise ValueError("fo baixa demais (raio > 60 mm).")
    for _ in range(60):
        m = 0.5 * (lo + hi)
        if estimativa_csrr(GeometriaCSRR(m, c, d, abertura), sub)[2] > alvo:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def capacitancia_interdigital(n_dedos, comp_mm, largura_mm, sub: Substrato):
    """Capacitor interdigital, fórmula de Bahl (Lumped Elements for RF and Microwave Circuits).

    C [pF] = (εr + 1)·ℓ[µm]·[(N - 3)·A1 + A2], A1 = 4,409·tanh(0,55(h/W)^0,45)·1e-6,
    A2 = 9,92·tanh(0,52(h/W)^0,5)·1e-6 [pF/µm]; W = largura dos dedos, espaçamento = W.
    Válida para N ≥ 3. Retorna F.
    """
    if n_dedos < 3:
        raise ValueError("A fórmula de Bahl vale para N ≥ 3 dedos.")
    hw = sub.h_mm / largura_mm
    a1 = 4.409 * np.tanh(0.55 * hw**0.45) * 1e-6
    a2 = 9.92 * np.tanh(0.52 * hw**0.5) * 1e-6
    return (sub.er + 1) * comp_mm * 1e3 * ((n_dedos - 3) * a1 + a2) * 1e-12


def anel_c(cx, cy, r_int, r_ext, abertura, lado, n=72):
    """Fenda em C (anel com corte), corte centrado no ângulo `lado` [rad]."""
    meia = np.arcsin(min(abertura / 2 / r_ext, 0.99))
    meia_i = np.arcsin(min(abertura / 2 / r_int, 0.99))
    t_ext = np.linspace(lado + meia, lado + 2 * np.pi - meia, n)
    t_int = np.linspace(lado + 2 * np.pi - meia_i, lado + meia_i, n)
    ext = np.column_stack([cx + r_ext * np.cos(t_ext), cy + r_ext * np.sin(t_ext)])
    inn = np.column_stack([cx + r_int * np.cos(t_int), cy + r_int * np.sin(t_int)])
    return np.vstack([ext, inn])


def fendas_csrr(cx, cy, geo: GeometriaCSRR):
    """Dois anéis em C concêntricos, cortes em lados opostos (Fig. 3.16), ao longo do eixo x."""
    r1 = geo.r_ext
    r2 = geo.r_ext - geo.c - geo.d
    return [anel_c(cx, cy, r1 - geo.c, r1, geo.abertura, 0.0),
            anel_c(cx, cy, r2 - geo.c, r2, geo.abertura, np.pi)]
