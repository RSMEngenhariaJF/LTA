"""Ressonadores de linha em microstrip: Q descarregado, circuito RLC e acoplamento por gap.

Referência: D. M. Pozar, Microwave Engineering, 4ª ed., Wiley, 2012:
  - seção 6.2: λ/2 em aberto (≈ RLC paralelo, eq. 6.34) e λ/4 em curto; Q0 = β/(2α) (eq. 6.35, Ex. 6.2);
  - seção 6.6: ressonador acoplado por gap (eqs. 6.76-6.82, Ex. 6.6): acoplamento crítico com
    bc = √(π/(2Q0)), C = bc/(ωZ0), e ressonância deslocada pela raiz de tan βℓ + bc = 0
    (para o λ/4 em curto, cot βℓ = bc).
O gap físico vem do modelo de Garg-Bahl (descontinuidades.gap).
"""
from dataclasses import dataclass, field

import numpy as np

from .descontinuidades import (capacitancia_ponta, extensao_ponta_aberta, gap, gap_para_capacitancia,
                               pozar_alfa_condutor, q_radiacao)
from .microstrip import C0, NP_DB, Substrato, alfa_condutor, alfa_dieletrico, analise, sintese

TIPOS = ("λ/2 em aberto", "λ/4 em curto")


@dataclass
class Ressonador:
    tipo: str
    f0: float
    z0: float
    w: float          # mm
    eeff: float
    comp: float       # mm (físico, já descontada a ponta aberta se corrigido)
    alfa_c: float     # Np/m em f0
    alfa_d: float
    texto: list = field(default_factory=list)

    @property
    def beta(self):
        return 2 * np.pi * self.f0 * np.sqrt(self.eeff) / C0

    @property
    def alfa(self):
        return self.alfa_c + self.alfa_d

    @property
    def q0(self):
        return self.beta / (2 * self.alfa)

    @property
    def qc(self):
        return self.beta / (2 * self.alfa_c)

    @property
    def qd(self):
        return self.beta / (2 * self.alfa_d) if self.alfa_d > 0 else np.inf

    @property
    def comp_eletrico(self):
        """Comprimento elétrico nominal (π ou π/2) do ressonador."""
        return np.pi if self.tipo == TIPOS[0] else np.pi / 2

    def rlc(self):
        """(R, L, C) do RLC paralelo equivalente perto de f0 (eq. 6.34; λ/4 em curto: C = π/(4ω0Z0))."""
        w0 = 2 * np.pi * self.f0
        l_m = self.comp_eletrico / self.beta
        r = self.z0 / (self.alfa * l_m)
        c = np.pi / (2 * w0 * self.z0) if self.tipo == TIPOS[0] else np.pi / (4 * w0 * self.z0)
        return r, 1 / (w0**2 * c), c

    def zin(self, f):
        """Impedância de entrada do ressonador (sem acoplamento), modelo de linha com α ∝ √f."""
        f = np.asarray(f, dtype=float)
        a = self.alfa_c * np.sqrt(f / self.f0) + self.alfa_d * f / self.f0
        g = a + 1j * 2 * np.pi * f * np.sqrt(self.eeff) / C0
        l = self.comp_eletrico / self.beta
        t = np.tanh(g * l)
        return self.z0 / t if self.tipo == TIPOS[0] else self.z0 * t


def projeto(tipo, f0, z0, sub: Substrato, corrigir_ponta=True):
    """λ/2 em aberto ou λ/4 em curto em microstrip, com perdas do modelo da aba Substrato."""
    if tipo not in TIPOS:
        raise ValueError(f"Tipo desconhecido: {tipo}")
    w = float(sintese(z0, sub))
    _, ee = analise(w, sub)
    ee = float(ee)
    lg = C0 / (f0 * np.sqrt(ee)) * 1e3
    comp = lg / 2 if tipo == TIPOS[0] else lg / 4
    pontas = 2 if tipo == TIPOS[0] else 1
    dl = extensao_ponta_aberta(w, sub)
    if corrigir_ponta:
        comp -= pontas * dl
    ac = float(alfa_condutor(w, sub, f0)[0, 0])
    ad = float(alfa_dieletrico(w, sub, f0)[0, 0])
    r = Ressonador(tipo, f0, z0, w, ee, comp, ac, ad)
    r.texto = [f"{tipo} em {sub.nome}: W = {w:.4f} mm, εeff = {ee:.4f}, λg = {lg:.3f} mm",
               f"Comprimento físico = {comp:.3f} mm" + (f" (descontados {pontas}×Δl = {pontas * dl:.3f} mm "
                                                         "da ponta aberta)" if corrigir_ponta else
                                                         f" (sem corrigir Δl = {dl:.3f} mm por ponta)"),
               f"αc = {ac * NP_DB / 100:.4f} dB/cm   αd = {ad * NP_DB / 100:.4f} dB/cm",
               f"Q0 = β/(2α) = {r.q0:.1f}   (Qc = {r.qc:.1f}, Qd = {r.qd:.1f}; 1/Q0 = 1/Qc + 1/Qd)",
               f"Comparação: αc ≈ Rs/(Z0·W) do Pozar (eq. 3.199) daria Q0 = "
               f"{r.beta / (2 * (float(pozar_alfa_condutor(w, sub, f0)) + ad)):.1f}",
               f"Radiação (Garg-Bahl eq. 2.142, λ/2 aberto): Qr ≈ {float(q_radiacao(z0, sub, f0)):.0f} → "
               f"Q total ≈ {1 / (1 / r.q0 + 1 / float(q_radiacao(z0, sub, f0))):.1f} (não incluída no Q0 acima)"]
    return r


@dataclass
class AcoplamentoGap:
    bc: float          # susceptância normalizada do gap
    c: float           # F
    f_carregada: float  # ressonância com o gap [Hz]
    qe: float
    ql: float
    g: float           # coeficiente de acoplamento Q0/Qe
    s_gap: float = None   # mm (se realizável)
    cp_gap: float = 0.0
    texto: list = field(default_factory=list)


def acoplamento_gap(r: Ressonador, sub: Substrato, c=None):
    """Gap série entre linha de alimentação Z0 e o ressonador; c = None dá o acoplamento crítico."""
    w0 = 2 * np.pi * r.f0
    fator = np.pi / 2 if r.tipo == TIPOS[0] else np.pi / 4      # Qe = π/(2bc²) (eq. 6.81) e π/(4bc²)
    bc = np.sqrt(fator / r.q0) if c is None else w0 * c * r.z0
    c = bc / (w0 * r.z0) if c is None else c
    # tan βℓ + bc = 0 (λ/2 aberto) ou cot βℓ = bc (λ/4 em curto): menor raiz abaixo da nominal
    bl = np.pi - np.arctan(bc) if r.tipo == TIPOS[0] else np.arctan(1 / bc)
    f1 = r.f0 * bl / r.comp_eletrico
    qe = fator / bc**2
    ql = 1 / (1 / r.q0 + 1 / qe)
    a = AcoplamentoGap(bc, c, f1, qe, ql, r.q0 / qe)
    try:
        a.s_gap = gap_para_capacitancia(r.w, c, sub)
        a.cp_gap = gap(r.w, a.s_gap, sub)[1]
    except ValueError as e:
        a.texto.append(f"⚠ {e}")
    estado = "crítico" if np.isclose(a.g, 1, rtol=1e-3) else ("sobreacoplado" if a.g > 1 else "subacoplado")
    a.texto = [f"bc = {bc:.4f}  →  C = {c * 1e12:.4f} pF  (acoplamento {estado}, g = Q0/Qe = {a.g:.3f})",
               f"Qe = {qe:.1f}   QL = {ql:.1f}",
               f"Ressonância com o gap: {f1 / 1e9:.4f} GHz ({100 * (f1 / r.f0 - 1):+.2f}% de f0)"] + a.texto
    if a.s_gap is not None:
        a.texto.append(f"Gap (Garg-Bahl): s = {a.s_gap * 1e3:.1f} µm, Cp = {a.cp_gap * 1e15:.2f} fF por lado")
    return a


def zin_acoplado(f, r: Ressonador, a: AcoplamentoGap):
    """Zin vista da linha: gap série (C) + ressonador."""
    w = 2 * np.pi * np.asarray(f, dtype=float)
    return 1 / (1j * w * a.c) + r.zin(f)


def ponta_capacitancia(r: Ressonador, sub: Substrato):
    """Capacitância da ponta aberta [F] (informativo)."""
    return capacitancia_ponta(r.w, sub)
