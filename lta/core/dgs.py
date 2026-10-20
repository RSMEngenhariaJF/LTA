"""Estruturas com plano de terra defeituoso (DGS): modelo, linha periódica, filtro e extração.

Referência: R. Garg, I. Bahl e M. Bianchi, Microstrip Lines and Slotlines, 4ª ed., 2024, cap. 6:
  - célula DGS (haltere) como LC paralelo em série com a linha: C = ωc/(2Z0(ω0² - ωc²)),
    L = 1/(ω0²C) (eq. 6.8), a partir da rejeição f0 e do corte de 3 dB fc da simulação EM;
    perdas por R em paralelo, eq. 6.9;
  - extração do modelo π (RLC série em paralelo, Cp‖Rp nas derivações) a partir de S (eqs. 6.10-6.13);
  - linha DGS periódica: banda de rejeição e onda lenta (seções 6.2.1-6.2.2), analisada como
    célula ABCD (Martín, cap. 2);
  - passa-baixas com DGS e stubs (seção 6.4.1): a DGS faz o indutor série que exigiria uma
    linha muito estreita.
"""
from dataclasses import dataclass

import numpy as np

from . import abcd
from .microstrip import C0, Substrato, analise, parametros_linha, sintese


@dataclass
class CelulaDGS:
    L: float            # H
    C: float            # F
    R: float = np.inf   # Ω (perdas), em paralelo

    @property
    def f0(self):
        return 1 / (2 * np.pi * np.sqrt(self.L * self.C))

    def z(self, f):
        w = 2 * np.pi * np.asarray(f, dtype=float)
        y = 1 / (1j * w * self.L) + 1j * w * self.C + (0 if np.isinf(self.R) else 1 / self.R)
        return 1 / y

    def abcd(self, f):
        return abcd.series(self.z(f))


def lc_de_resposta(f0, fc, z0=50.0):
    """(L, C) da célula DGS a partir da rejeição f0 e do corte de 3 dB fc < f0 (eq. 6.8)."""
    if not 0 < fc < f0:
        raise ValueError("É preciso 0 < fc < f0.")
    w0, wc = 2 * np.pi * f0, 2 * np.pi * fc
    c = wc / (2 * z0 * (w0**2 - wc**2))
    return 1 / (w0**2 * c), c


def r_de_s11(s11_f0, f0, L, C, z0=50.0):
    """R da célula com perdas a partir de |S11| em f0 (eq. 6.9)."""
    w0 = 2 * np.pi * f0
    x = w0 * C - 1 / (w0 * L)
    return 2 * z0 / np.sqrt(1 / abs(s11_f0) ** 2 - (2 * z0 * x) ** 2 - 1)


def celula_periodica(f, cel: CelulaDGS, periodo_mm, w_mm, sub: Substrato, perdas=False):
    """Meia linha + DGS + meia linha (célula simétrica de período p)."""
    f = np.asarray(f, dtype=float)
    z, ee, a = (x[:, 0] for x in parametros_linha(w_mm, sub, f, perdas, False))
    th = (2 * np.pi * f * np.sqrt(ee) / C0 - 1j * a) * periodo_mm * 1e-3 / 2
    meia = abcd.linha(z, th)
    return abcd.cascata(meia, cel.abcd(f), meia)


def passa_baixas(fc, z0=50.0, n=3, ripple_db=0.01, aproximacao="Chebyshev"):
    """Elementos do passa-baixas L-C-L com DGS nos indutores e capacitores para stubs (seção 6.4.1).

    Retorna lista [('L'|'C', valor)] na ordem; os L serão células DGS e os C, stubs abertos."""
    from .filtros import prototipo
    g = prototipo(aproximacao, n, ripple_db)
    wc = 2 * np.pi * fc
    out = []
    for i, gi in enumerate(g[1:-1]):
        out.append(("L", gi * z0 / wc) if i % 2 == 0 else ("C", gi / (z0 * wc)))
    return out


def stub_para_capacitor(c, f, sub: Substrato, z_stub=50.0, pares=True):
    """Comprimento [mm] de stub aberto (ou par simétrico de stubs) que dá a capacitância c em f."""
    w = float(sintese(z_stub, sub))
    _, ee = analise(w, sub)
    b = 2 * np.pi * f * c / (2 if pares else 1)
    bl = np.arctan(b * z_stub)
    return float(bl * C0 / (2 * np.pi * f * np.sqrt(float(ee))) * 1e3), w


def abcd_passa_baixas(f, elementos, z0, sub: Substrato = None, stubs=None, f_proj=None, f0_dgs=None):
    """Rede ideal (L e C concentrados) ou com DGS (LC paralelo ressonante em f0_dgs) e stubs reais."""
    f = np.asarray(f, dtype=float)
    w = 2 * np.pi * f
    mats = []
    for k, (tipo, v) in enumerate(elementos):
        if tipo == "L":
            if f0_dgs:
                cd = 1 / ((2 * np.pi * f0_dgs) ** 2 * v)
                mats.append(CelulaDGS(v, cd).abcd(f))
            else:
                mats.append(abcd.series(1j * w * v))
        else:
            if stubs:
                comp, wst, pares = stubs[k]
                z, ee, _ = (x[:, 0] for x in parametros_linha(wst, sub, f, False, False))
                y = np.tan(2 * np.pi * f * np.sqrt(ee) / C0 * comp * 1e-3) / z * (2 if pares else 1)
                mats.append(abcd.shunt(1j * y))
            else:
                mats.append(abcd.shunt(1j * w * v))
    return abcd.cascata(*mats)


@dataclass
class ModeloPi:
    R: float
    L: float
    C: float
    Rp: float
    Cp: float
    f0: float


def extrair_pi(f, s, z0=50.0, f1=None):
    """Modelo π da DGS (eqs. 6.10-6.13): Ya = 1/B (RLC paralelo, ramo série) e Yb = (A - 1)/B
    (Cp‖Rp em cada derivação). f0 é onde Im(Ya) = 0; f1 (padrão: f0/2) fixa C e Cp."""
    m = abcd.s_para_abcd(s, z0)
    ya = 1 / m[:, 0, 1]
    yb = (m[:, 0, 0] - 1) / m[:, 0, 1]
    im = ya.imag
    trocas = np.flatnonzero(np.sign(im[:-1]) != np.sign(im[1:]))
    if not len(trocas):
        raise ValueError("Im(Ya) não cruza zero: a faixa medida não contém a ressonância da DGS.")
    i = trocas[0]
    f0 = f[i] - im[i] * (f[i + 1] - f[i]) / (im[i + 1] - im[i])
    f1 = f1 or f0 / 2
    j = int(np.argmin(np.abs(f - f1)))
    w0, w1 = 2 * np.pi * f0, 2 * np.pi * f[j]
    c = im[j] / (w0 * (w1 / w0 - w0 / w1))
    k0 = int(np.argmin(np.abs(f - f0)))
    return ModeloPi(R=float(1 / ya[k0].real), L=float(1 / (w0**2 * c)), C=float(c),
                    Rp=float(1 / yb[j].real) if yb[j].real > 0 else np.inf, Cp=float(yb[j].imag / w1),
                    f0=float(f0))


def abcd_pi(f, m: ModeloPi):
    w = 2 * np.pi * np.asarray(f, dtype=float)
    ya = 1 / m.R + 1j * w * m.C + 1 / (1j * w * m.L)
    yb = (0 if np.isinf(m.Rp) else 1 / m.Rp) + 1j * w * m.Cp
    return abcd.cascata(abcd.shunt(yb), abcd.series(1 / ya), abcd.shunt(yb))
