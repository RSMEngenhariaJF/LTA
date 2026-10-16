"""Linhas diferenciais com supressão do modo comum por CSRR ou DS-CSRR.

Referência: F. Martín (2015), seção 6.3 (Fig. 6.12, eq. 6.11, Tabela 6.1) e
seção 3.5.2.3 (acoplamento entre ressonadores, eqs. 3.93-3.97).

Modo par (comum, parede magnética no plano de simetria): cada linha vê Le e o
acoplamento Ce ao meio CSRR, que tem indutância 2Lc e capacitância Cc/2; ressonadores
vizinhos se acoplam por CR (dois capacitores 2CR em série por célula, Fig. 3.44b).
Modo ímpar (diferencial, parede elétrica): linha comum com Lo e Co; o CSRR não é excitado.
"""
import numpy as np

from . import abcd
from .dispositivos import Rede
from .microstrip import C0

# Tabela 6.1: (f0 [GHz], Le [nH], Ce [pF], Lc [nH], Cc [pF], CR [pF], FBW [%])
TABELA_6_1 = {
    "CSRR 1,5 GHz": (1.5, 6.3, 1.1, 2.1, 3.2, 0.05, 30.7),
    "CSRR 3 GHz": (3.0, 4.2, 0.5, 1.0, 2.0, 0.03, 26.4),
    "CSRR 6 GHz": (6.0, 3.0, 0.2, 0.4, 1.4, 0.04, 20.7),
    "DS-CSRR 1,5 GHz": (1.5, 15.9, 1.6, 1.0, 8.6, 0.32, 32.6),
    "DS-CSRR 3 GHz": (3.0, 9.0, 0.8, 0.45, 4.8, 0.09, 29.4),
    "DS-CSRR 6 GHz": (6.0, 6.9, 0.3, 0.15, 3.6, 0.05, 23.6),
}


# Lado do ressonador (Fig. 6.13) + 0,2 mm entre ressonadores = período da célula [mm]
PERIODO_FIG_6_13 = {"CSRR 1,5 GHz": 7.5, "CSRR 3 GHz": 4.5, "CSRR 6 GHz": 2.8,
                    "DS-CSRR 1,5 GHz": 14.0, "DS-CSRR 3 GHz": 7.7, "DS-CSRR 6 GHz": 4.5}


def linha_impar(z_impar, eeff, periodo_mm):
    """(Lo, Co) por célula de um trecho de linha: L = Z·√εeff·ℓ/c, C = √εeff·ℓ/(Z·c)."""
    atraso = np.sqrt(eeff) * periodo_mm * 1e-3 / C0
    return z_impar * atraso, atraso / z_impar


def zero_comum(le, ce, lc, cc):
    """Zero de transmissão do modo comum: fz = 1/(2π√(2Lc(Ce + Cc/2))) (eq. 6.11)."""
    return 1 / (2 * np.pi * np.sqrt(2 * lc * (ce + cc / 2)))


def _ramos(f, le, ce, lc, cc, cr):
    w = 2 * np.pi * np.asarray(f, dtype=float)
    zs1 = 1j * w * le / 2                      # meia indutância da linha
    zs2 = 1 / (1j * w * 2 * cr)               # 2CR de cada lado do ressonador
    zp1 = 1 / (1j * w * ce)                    # acoplamento linha-ressonador
    zp2 = 1 / (1j * w * cc / 2 + 1 / (1j * w * 2 * lc))   # meio CSRR (2Lc ‖ Cc/2)
    return zs1, zs2, zp1, zp2


def dispersao_comum(f, le, ce, lc, cc, cr):
    """cosh(γl) dos dois modos de Bloch da célula de 4 portas (eqs. 3.96-3.97)."""
    zs1, zs2, zp1, zp2 = _ramos(f, le, ce, lc, cc, cr)
    a11 = 1 + zs1 / zp1
    a12 = -zs1 / zp1
    a21 = -zs2 / zp1
    a22 = 1 + zs2 / zp1 + zs2 / zp2
    raiz = np.sqrt((a11 - a22) ** 2 + 4 * a12 * a21 + 0j)
    return 0.5 * (a11 + a22 + raiz), 0.5 * (a11 + a22 - raiz)


def banda_rejeicao(le, ce, lc, cc, cr, n_pts=40001):
    """Faixa em que nenhum modo de Bloch propaga em torno de fz (estrutura infinita).

    Retorna (f_ini, f_fim, FBW = (f_fim - f_ini)/√(f_ini·f_fim)). É a banda máxima que um
    filtro com muitas células atinge (livro: seis células já ficam perto dela).
    """
    fz = zero_comum(le, ce, lc, cc)
    f = np.linspace(0.3 * fz, 3 * fz, n_pts)
    c1, c2 = dispersao_comum(f, le, ce, lc, cc, cr)
    propaga = lambda c: (np.abs(c.imag) < 1e-9) & (np.abs(c.real) <= 1)   # noqa: E731
    rej = ~(propaga(c1) | propaga(c2))
    i = int(np.argmin(np.abs(f - fz)))
    if not rej[i]:
        raise ValueError("Não há banda rejeitada em torno de fz com esses elementos.")
    a = b = i
    while a > 0 and rej[a - 1]:
        a -= 1
    while b < len(f) - 1 and rej[b + 1]:
        b += 1
    return f[a], f[b], (f[b] - f[a]) / np.sqrt(f[a] * f[b])


def s_comum(f, n, le, ce, lc, cc, cr, z0=50.0):
    """S do modo comum (meio circuito) de n células; capacitores de acoplamento das pontas abertos."""
    zs1, zs2, zp1, zp2 = _ramos(f, le, ce, lc, cc, cr)
    # nós: linha 0..n (fronteiras), meio de linha n+1..2n, ressonadores 2n+1..3n
    r = Rede(3 * n + 1, [0, n], f)
    for k in range(n):
        meio, res = n + 1 + k, 2 * n + 1 + k
        r.impedancia(k, meio, zs1)
        r.impedancia(meio, k + 1, zs1)
        r.impedancia(meio, res, zp1)
        r.terra(res, zp2)
        if k < n - 1:
            r.impedancia(res, res + 1, 2 * zs2)        # dois 2CR em série = CR entre vizinhos
    return r.s(z0)


def s_diferencial(f, n, lo, co, z0=50.0):
    """S do modo diferencial (meio circuito): n células T de linha comum (Lo, Co)."""
    m = abcd.potencia(abcd.celula(f, "T", ls=lo, cp=co), n)
    return abcd.abcd_para_s(m, z0)
