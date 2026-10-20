"""Filtro combline com capacitores concentrados (Steer, Microwave and RF Design, seção 10.14).

Procedimento (generalizado do exemplo de 3ª ordem do livro, Figs. 10-53 a 10-69):
  1. protótipo passa-baixas g; cada indutor série vira inversor-capacitor-inversor (J = 1), de modo
     que todos os elementos são capacitores em derivação Ck = gk separados por inversores;
  2. passa-faixa (α = 1/Δ) e escala para Z0: ressonadores C# = α·gk/(ω0·Z0), L# = Z0/(α·gk·ω0),
     inversores K = Z0;
  3. cada ressonador LC vira capacitor C0 + stub em curto Zs ressonante em fr = 2f0 (45° em f0),
     igualando admitância e derivada em f0: Zs = (1 + π/2)/(2ω0C#), C0 = 1/(ω0·Zs);
  4. escala dos nós para que todos os ressonadores fiquem iguais ao primeiro: ak = C1/Ck,
     K'k,k+1 = K/√(ak·ak+1), Zs,k → Zs,k/ak;
  5. cada inversor interno vira π de stubs em curto (série +jK', derivações -jK'); as derivações
     negativas se somam aos stubs dos ressonadores: 1/Zk = 1/Zs,k - nk/K';
  6. escala de impedância interna s (stub central = Zalvo): stubs × s, capacitores ÷ s, nível R' = s·Z0;
  7. inversores de entrada/saída K = √(Z0·R') realizados por Cb série e Ca (negativo) absorvido
     no ressonador da ponta: ωCb = 1/√(K² - Z0²), Ca = -Cb/(1 + (ωCbZ0)²).
A realização física são N linhas acopladas em curto numa ponta (vias) e com Ck na outra; as
linhas são estimadas por pares (Z0e ≈ Zk, 1/Z0o = 1/Zk + 2/Zk,k+1) e devem ser ajustadas em EM.
"""
from dataclasses import dataclass, field

import numpy as np

from . import filtros
from .dispositivos import Rede


@dataclass
class ProjetoCombline:
    f0: float
    z0: float
    c: list                    # capacitores dos ressonadores [F]
    z_stub: list               # stubs em curto em derivação [Ω] (λ/4 em fr)
    z_acop: list               # stubs série de acoplamento entre ressonadores [Ω]
    cb: float                  # capacitor série de entrada/saída [F]
    fr: float                  # frequência comensurada (stubs λ/4)
    nivel: float               # impedância interna R' [Ω]
    texto: list = field(default_factory=list)


def projeto(f0, fbw, n=3, aproximacao="Chebyshev", ripple_db=0.0432, z0=50.0, z_alvo=80.0, ripple_eps=None):
    """Combline de ordem n; ripple_eps (ε) tem prioridade sobre ripple_db (ε = 0,1 ↔ 0,0432 dB)."""
    if ripple_eps is not None:
        ripple_db = 10 * np.log10(1 + ripple_eps**2)
    g = filtros.prototipo(aproximacao, n, ripple_db)
    if not np.isclose(g[-1], 1.0):
        raise ValueError("Use ordem ímpar (Chebyshev) ou Butterworth: o combline exige terminações iguais.")
    w0 = 2 * np.pi * f0
    alfa = 1 / fbw
    gk = np.asarray(g[1:-1], dtype=float)
    c_lumped = alfa * gk / (w0 * z0)
    zs = (1 + np.pi / 2) / (2 * w0 * c_lumped)                    # passo 3
    c0 = 1 / (w0 * zs)
    a = c0[0] / c0                                                # passo 4
    zs = zs / a
    c0 = c0 * a
    k = z0 / np.sqrt(a[:-1] * a[1:])
    zk = []
    for i in range(n):                                            # passo 5
        y = 1 / zs[i] - (1 / k[i - 1] if i > 0 else 0) - (1 / k[i] if i < n - 1 else 0)
        zk.append(1 / y)
    zk = np.array(zk)
    s = z_alvo / zk[n // 2]                                       # passo 6
    zk, kk, c0 = zk * s, k * s, c0 / s
    nivel = z0 * s
    kio = np.sqrt(z0 * nivel)                                     # passo 7
    cb = 1 / (w0 * np.sqrt(kio**2 - z0**2))
    ca = -cb / (1 + (w0 * cb * z0) ** 2)
    c = c0.copy()
    c[0] += ca
    c[-1] += ca
    p = ProjetoCombline(f0, z0, list(c), list(zk), list(kk), float(cb), 2 * f0, float(nivel))
    p.texto = [f"Protótipo {aproximacao} N = {n}" + (f", ripple {ripple_db:.4f} dB" if aproximacao == "Chebyshev"
                                                    else "") + ": g = " + ", ".join(f"{x:.5f}" for x in gk),
               f"Stubs em curto λ/4 em fr = {2 * f0 / 1e9:g} GHz (45° em f0); nível interno R' = {nivel:.2f} Ω",
               "Ressonadores: " + "; ".join(f"C{i + 1} = {ci * 1e12:.4f} pF, Z{i + 1} = {zi:.3f} Ω"
                                            for i, (ci, zi) in enumerate(zip(c, zk))),
               "Acoplamentos (stubs série): " + ", ".join(f"{x:.3f} Ω" for x in kk),
               f"Capacitores de entrada/saída: Cb = {cb * 1e12:.4f} pF (Ca = {ca * 1e12:.4f} pF absorvido)"]
    return p


def s_combline(f, p: ProjetoCombline, c_extra=None):
    """S do circuito de stubs (ideais, TEM) com capacitores: nós 0 e n+1 são as portas."""
    f = np.asarray(f, dtype=float)
    n = len(p.c)
    w = 2 * np.pi * f
    t = np.tan(np.pi / 2 * f / p.fr)
    r = Rede(n + 2, [0, n + 1], f)
    r.impedancia(0, 1, 1 / (1j * w * p.cb))
    r.impedancia(n, n + 1, 1 / (1j * w * p.cb))
    for i in range(n):
        y = 1j * w * p.c[i] + 1 / (1j * p.z_stub[i] * t)
        r.terra(i + 1, 1 / y)
        if i < n - 1:
            r.impedancia(i + 1, i + 2, 1j * p.z_acop[i] * t)
    return r.s(p.z0)


def pares_acoplados(p: ProjetoCombline):
    """(Z0e, Z0o) estimados para cada par de linhas vizinhas (meio homogêneo, ver docstring)."""
    out = []
    for i in range(len(p.c) - 1):
        zs = 0.5 * (p.z_stub[i] + p.z_stub[i + 1])
        out.append((zs, 1 / (1 / zs + 2 / p.z_acop[i])))
    return out
