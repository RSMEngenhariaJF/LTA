"""Transição coaxial-microstrip (conector SMA): modelo, compensação e de-embedding.

Referência: R. Garg, I. Bahl e M. Bianchi, Microstrip Lines and Slotlines, 4ª ed., Artech House,
2024, seção 1.4.1 (Fig. 1.20: a interface é Ls em série e Cs em derivação; projetar minimizando
essas reatâncias ou compensando-as com √(Ls/Cs) = Z0) e Tabela 2.5 (conectores/transições com
VSWR típico de 1,065 em 3 GHz e 1,10 em 10 GHz num orçamento de erros).

Os valores de Ls e Cs dependem do conector e do pad: medir (ou simular em EM) e ajustar. Os
predefinidos abaixo são ordens de grandeza para começar, não dados de fabricante.
"""
from dataclasses import dataclass

import numpy as np

from . import abcd
from .microstrip import C0, Substrato, analise, parametros_linha, sintese


@dataclass
class Transicao:
    nome: str
    ls: float = 0.0     # H, indutância série da interface
    cs: float = 0.0     # F, capacitância em derivação da interface

    @property
    def c_excesso(self):
        """Capacitância líquida em excesso em Z0 = 50 Ω (positiva: capacitiva)."""
        return self.cs - self.ls / 50.0**2


def cs_de_vswr(vswr, f, z0=50.0):
    """Cs [F] de uma transição puramente capacitiva com o VSWR dado em f: |Γ| = x/√(1 + x²), x = ωCsZ0/2."""
    g = (vswr - 1) / (vswr + 1)
    x = g / np.sqrt(1 - g**2)
    return float(2 * x / (2 * np.pi * f * z0))


PREDEFINIDOS = {
    "Sem transição (ideal)": Transicao("ideal"),
    "Típica, VSWR 1,10 em 10 GHz (Tabela 2.5)": Transicao("tipica_10GHz", 0.0, cs_de_vswr(1.10, 10e9)),
    "Típica, VSWR 1,065 em 3 GHz (Tabela 2.5)": Transicao("tipica_3GHz", 0.0, cs_de_vswr(1.065, 3e9)),
}


def abcd_transicao(f, t: Transicao, lado_coaxial_primeiro=True):
    """ABCD da interface: coaxial → Ls série → Cs derivação → microstrip (ou o inverso)."""
    w = 2 * np.pi * np.asarray(f, dtype=float)
    serie = abcd.series(1j * w * t.ls)
    deriv = abcd.shunt(1j * w * t.cs)
    return abcd.cascata(serie, deriv) if lado_coaxial_primeiro else abcd.cascata(deriv, serie)


def compensacao(t: Transicao, sub: Substrato, z0=50.0, z_comp=None):
    """Trecho de linha (Zc, ℓ [mm]) junto ao pad que cancela a reatância líquida da transição.

    Um trecho curto de Zc equivale a L = Zc·θ/ω em série e C = θ/(ωZc) em derivação; para que
    trecho + transição pareçam uma linha de Z0 é preciso θ/ω·(Zc - Z0²/Zc) = Z0²·Cex, com
    Cex = Cs - Ls/Z0². Cex > 0 pede trecho de alta impedância (mais estreito); Cex < 0, de baixa.
    """
    cex = t.cs - t.ls / z0**2
    if abs(cex) < 1e-18:
        return None
    if z_comp is None:
        z_comp = 1.6 * z0 if cex > 0 else z0 / 1.6
    if (cex > 0) != (z_comp > z0):
        raise ValueError("Use Zc > Z0 para transição capacitiva e Zc < Z0 para indutiva.")
    w = float(sintese(z_comp, sub))
    _, ee = analise(w, sub)
    comp = z0**2 * cex * C0 / (np.sqrt(float(ee)) * (z_comp - z0**2 / z_comp)) * 1e3
    return z_comp, w, float(comp)


def abcd_linha_ms(f, w_mm, comp_mm, sub, perdas=False):
    f = np.asarray(f, dtype=float)
    z, ee, a = (x[:, 0] for x in parametros_linha(w_mm, sub, f, perdas, False))
    return abcd.linha(z, 2 * np.pi * f * np.sqrt(ee) / C0 * comp_mm * 1e-3 - 1j * a * comp_mm * 1e-3)


def dispositivo_de_teste(f, t: Transicao, sub: Substrato, comp_mm, z0=50.0, perdas=False, comp_trecho=None):
    """ABCD de um fixture coaxial-linha-coaxial: transição, [trecho de compensação], linha de Z0
    de comp_mm, [trecho], transição espelhada."""
    w0 = float(sintese(z0, sub))
    mats = [abcd_transicao(f, t, True)]
    if comp_trecho:
        mats.append(abcd_linha_ms(f, comp_trecho[1], comp_trecho[2], sub, perdas))
    mats.append(abcd_linha_ms(f, w0, comp_mm, sub, perdas))
    if comp_trecho:
        mats.append(abcd_linha_ms(f, comp_trecho[1], comp_trecho[2], sub, perdas))
    mats.append(abcd_transicao(f, t, False))
    return abcd.cascata(*mats)


def de_embutir(f, s_medido, t: Transicao, sub: Substrato, alim_mm=0.0, z0=50.0, perdas=False, comp_trecho=None):
    """Remove de S medido (2 portas) as duas transições e as linhas de alimentação de alim_mm:
    M_dut = A⁻¹·M_med·B⁻¹, com A = transição + [trecho] + alimentação e B o espelho."""
    w0 = float(sintese(z0, sub))
    a = [abcd_transicao(f, t, True)]
    if comp_trecho:
        a.append(abcd_linha_ms(f, comp_trecho[1], comp_trecho[2], sub, perdas))
    if alim_mm > 0:
        a.append(abcd_linha_ms(f, w0, alim_mm, sub, perdas))
    ma = abcd.cascata(*a)
    b = list(reversed([abcd_transicao(f, t, False)] + ([abcd_linha_ms(f, comp_trecho[1], comp_trecho[2], sub,
                                                                       perdas)] if comp_trecho else [])
                      + ([abcd_linha_ms(f, w0, alim_mm, sub, perdas)] if alim_mm > 0 else [])))
    mb = abcd.cascata(*b)
    m = abcd.s_para_abcd(s_medido, z0)
    return abcd.abcd_para_s(np.linalg.inv(ma) @ m @ np.linalg.inv(mb), z0)


def vswr(s11):
    g = np.abs(s11)
    return (1 + g) / (1 - np.minimum(g, 0.999999))
