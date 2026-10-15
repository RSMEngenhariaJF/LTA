"""Linhas de transmissão CRLH (mão esquerda/direita compostas) e D-CRLH.

Referência: F. Martín (2015), seções 3.4.2-3.4.3 (eqs. 3.51-3.66).

Célula CRLH (Fig. 3.20, T): braços série LR/2 + 2CL; derivação LL || CR.
Célula D-CRLH (Fig. 3.23, T): braços série (LR/2 || 2CL); derivação LL + CR em série.
Frequências: ωR = 1/√(LR·CR), ωL = 1/√(LL·CL), ωs = 1/√(LR·CL), ωp = 1/√(LL·CR).
"""
from dataclasses import dataclass

import numpy as np

from . import abcd

TIPOS = ("CRLH", "D-CRLH")
MU0 = 4e-7 * np.pi
EPS0 = 8.8541878128e-12


@dataclass
class CelulaCRLH:
    LR: float   # H
    CR: float   # F
    LL: float   # H
    CL: float   # F
    tipo: str = "CRLH"

    def __post_init__(self):
        if min(self.LR, self.CR, self.LL, self.CL) <= 0:
            raise ValueError("LR, CR, LL e CL devem ser positivos.")
        if self.tipo not in TIPOS:
            raise ValueError(f"Tipo desconhecido: {self.tipo}")

    @property
    def wR(self):
        return 1 / np.sqrt(self.LR * self.CR)          # eq. 3.53

    @property
    def wL(self):
        return 1 / np.sqrt(self.LL * self.CL)          # eq. 3.54

    @property
    def ws(self):
        return 1 / np.sqrt(self.LR * self.CL)          # eq. 3.55

    @property
    def wp(self):
        return 1 / np.sqrt(self.LL * self.CR)          # eq. 3.56

    @property
    def balanceada(self):
        return np.isclose(self.ws, self.wp, rtol=1e-3)

    @property
    def gap(self):
        """(ωG1, ωG2) da banda proibida entre LH e RH (eqs. 3.57-3.58)."""
        return min(self.ws, self.wp), max(self.ws, self.wp)

    @property
    def z_r(self):
        return np.sqrt(self.LR / self.CR)

    @property
    def z_l(self):
        return np.sqrt(self.LL / self.CL)

    def ramos(self, f):
        """(Zs total da célula, Yp) por frequência."""
        w = 2 * np.pi * np.asarray(f, dtype=float)
        if self.tipo == "CRLH":
            zs = 1j * w * self.LR + 1 / (1j * w * self.CL)
            yp = 1j * w * self.CR + 1 / (1j * w * self.LL)
        else:
            zs = 1j * w * self.LR / (1 - w**2 * self.LR * self.CL)
            yp = 1 / (1j * w * self.LL + 1 / (1j * w * self.CR))
        return zs, yp

    def abcd(self, f):
        zs, yp = self.ramos(f)
        return abcd.cascata(abcd.series(zs / 2), abcd.shunt(yp), abcd.series(zs / 2))

    def cos_bl(self, f):
        """cos(βl) analítico (eq. 3.51 para CRLH, 3.64 para D-CRLH)."""
        w = 2 * np.pi * np.asarray(f, dtype=float)
        if self.tipo == "CRLH":
            return 1 - w**2 / (2 * self.wR**2) * (1 - self.ws**2 / w**2) * (1 - self.wp**2 / w**2)
        return 1 - self.wL**2 / (2 * w**2) / ((1 - self.ws**2 / w**2) * (1 - self.wp**2 / w**2))

    def zb(self, f):
        """Impedância de Bloch da célula T (eqs. 3.52 / 3.65), ramo com Re ≥ 0 (eq. 2.26)."""
        return abcd.impedancia_bloch(self.abcd(f))[0]

    def sinal(self, f):
        """+1 onde a banda é de mão direita (RH), -1 onde é de mão esquerda (LH).

        Na banda passante de uma célula sem perdas, a reatância série X e a
        susceptância em derivação B têm o mesmo sinal (cos βl = 1 - XB/2); são
        ambas negativas na banda LH (eq. 3.60 generalizada).
        """
        zs, _ = self.ramos(f)
        return np.where(zs.imag < 0, -1.0, 1.0)

    def beta_l(self, f):
        """βl com sinal (LH < 0), βl ∈ [-π, π], e αl ≥ 0 [Np]."""
        al, bl = abcd.dispersao(self.abcd(f))
        return self.sinal(f) * bl, al

    def meio_efetivo(self, f, l):
        """(εr,eff, μr,eff) relativos, eqs. 3.61-3.62 (CRLH), com l = comprimento da célula [m]."""
        w = 2 * np.pi * np.asarray(f, dtype=float)
        eps = (self.CR / l - 1 / (w**2 * self.LL * l)) / EPS0
        mu = (self.LR / l - 1 / (w**2 * self.CL * l)) / MU0
        return eps, mu

    def bandas(self):
        """Bordas analíticas das bandas passantes do modelo CRLH [Hz] (cos βl = ±1).

        cos βl = -1 dá as bordas externas (ωCL, ωCR, eq. 3.66 no caso balanceado);
        cos βl = +1 dá o gap (ωs, ωp).
        """
        if self.tipo != "CRLH":
            raise ValueError("Bordas analíticas só implementadas para a CRLH.")
        # cos βl = -1: w^4/wR^2 - (ws^2 + wp^2 + 4 wR^2... ) -> polinômio em x = w^2
        a = 1 / self.wR**2
        b = -(self.ws**2 + self.wp**2) / self.wR**2 - 4
        c = self.ws**2 * self.wp**2 / self.wR**2
        x = np.roots([a, b, c])
        externas = np.sort(np.sqrt(x.real[x.real > 0]))
        return dict(f_cl=externas[0] / (2 * np.pi), f_cr=externas[-1] / (2 * np.pi),
                    f_g1=self.gap[0] / (2 * np.pi), f_g2=self.gap[1] / (2 * np.pi))


def balanceada_por_banda(f_cl, f_cr, zb):
    """CRLH balanceada a partir das bordas da banda (eq. 3.66) e de ZB (eq. 3.63).

    ωCR - ωCL = 2ωR e ωCL·ωCR = ω0² (transição), logo ωL = ω0²/ωR.
    """
    if not 0 < f_cl < f_cr:
        raise ValueError("É preciso 0 < f_CL < f_CR.")
    w_cl, w_cr = 2 * np.pi * f_cl, 2 * np.pi * f_cr
    w_r = (w_cr - w_cl) / 2
    w_l = w_cl * w_cr / w_r
    return CelulaCRLH(zb / w_r, 1 / (zb * w_r), zb / w_l, 1 / (zb * w_l))


def balanceada_por_fases(f1, fase1, f2, fase2, zb):
    """CRLH balanceada que dá βl = fase1 em f1 e fase2 em f2 (rad, βl < 0 na banda LH).

    Aproximação de grande comprimento de onda da eq. 3.59 balanceada:
    βl = ω/ωR - ωL/ω (sistema linear em 1/ωR e ωL). É o procedimento usado para
    componentes dual-band (cap. 4). Use fases pequenas (|βl| ≲ π/2) por célula.
    """
    w1, w2 = 2 * np.pi * f1, 2 * np.pi * f2
    m = np.array([[w1, -1 / w1], [w2, -1 / w2]])
    inv_wr, w_l = np.linalg.solve(m, [fase1, fase2])
    if inv_wr <= 0 or w_l <= 0:
        raise ValueError("Especificação inviável: exige f1 < f2 e fase em f1 menor que em f2 "
                         "(βl cresce com f; LH tem βl < 0).")
    w_r = 1 / inv_wr
    return CelulaCRLH(zb / w_r, 1 / (zb * w_r), zb / w_l, 1 / (zb * w_l))


def frequencia_transicao(c: CelulaCRLH):
    """ω0 = √(ωs·ωp) / 2π [Hz]; na balanceada, ω0 = ωs = ωp e ZB(ω0) = √(LR/CR) = √(LL/CL)."""
    return np.sqrt(c.ws * c.wp) / (2 * np.pi)
