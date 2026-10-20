"""Casamento de impedância: rede L, stub simples e duplo, transformadores multisseção e tapers.

Referência: D. M. Pozar, Microwave Engineering, 4ª ed., Wiley, 2012, cap. 5:
  - rede L (seção 5.1, eqs. 5.3 e 5.6, Ex. 5.1);
  - stub simples em derivação e em série (5.2, Ex. 5.2-5.3) e stub duplo (5.3, eq. 5.22, Ex. 5.4);
  - transformadores binomial e Chebyshev (5.6-5.7, eqs. 5.53-5.64, Ex. 5.6-5.7);
  - tapers exponencial, triangular e Klopfenstein (5.8, eqs. 5.69-5.78, Ex. 5.8);
  - limite de Bode-Fano (5.9, Fig. 5.23);
  - redes Pi e T de três elementos com Q especificado (M. Steer, Microwave and RF Design, 2010,
    seção 8.6, Ex. 8.7): duas redes L em sequência através de uma resistência virtual RV.

Comprimentos elétricos em comprimentos de onda (λ) da linha; impedâncias em Ω.
A carga é modelada como R em série com L ou C, com ZL(f0) dada (como nos exemplos do livro).
"""
from dataclasses import dataclass, field
from math import comb

import numpy as np

from . import abcd


# ---------------------------------------------------------------------------
# Carga e utilidades
# ---------------------------------------------------------------------------
def carga_serie(zl0, f0):
    """ZL(f) de R em série com o reativo que dá ZL(f0) = zl0 (L se X > 0, C se X < 0)."""
    r, x = zl0.real, zl0.imag
    w0 = 2 * np.pi * f0

    def zl(f):
        w = 2 * np.pi * np.asarray(f, dtype=float)
        if x > 0:
            return r + 1j * w * (x / w0)
        if x < 0:
            return r + 1 / (1j * w * (-1 / (w0 * x)))
        return np.full(w.shape, r, dtype=complex)
    return zl


def descrever_carga(zl0, f0):
    x = zl0.imag
    w0 = 2 * np.pi * f0
    if x > 0:
        return f"R = {zl0.real:g} Ω em série com L = {x / w0 * 1e9:.4g} nH"
    if x < 0:
        return f"R = {zl0.real:g} Ω em série com C = {-1 / (w0 * x) * 1e12:.4g} pF"
    return f"R = {zl0.real:g} Ω"


def gama(zin, z0):
    return (zin - z0) / (zin + z0)


def zin_abcd(m, zl):
    """Impedância de entrada de uma rede ABCD terminada em zl."""
    return (m[..., 0, 0] * zl + m[..., 0, 1]) / (m[..., 1, 0] * zl + m[..., 1, 1])


def _transformar(z, z0, bl):
    """Impedância vista a uma distância elétrica bl [rad] de z numa linha z0."""
    t = np.tan(bl)
    return z0 * (z + 1j * z0 * t) / (z0 + 1j * z * t)


def _comp_reativo(valor, f0, serie):
    """Elemento concentrado que realiza X (série) ou B (derivação) em f0: ('L'|'C', valor SI)."""
    w0 = 2 * np.pi * f0
    if serie:
        return ("L", valor / w0) if valor > 0 else ("C", -1 / (w0 * valor))
    return ("C", valor / w0) if valor > 0 else ("L", -1 / (w0 * valor))


def fmt_comp(comp):
    tipo, v = comp
    return f"L = {v * 1e9:.4g} nH" if tipo == "L" else f"C = {v * 1e12:.4g} pF"


def _impedancia_comp(comp, f):
    tipo, v = comp
    w = 2 * np.pi * np.asarray(f, dtype=float)
    return 1j * w * v if tipo == "L" else 1 / (1j * w * v)


# ---------------------------------------------------------------------------
# Rede L (elementos concentrados)
# ---------------------------------------------------------------------------
@dataclass
class SolucaoL:
    topologia: str          # "derivação na carga" (Fig. 5.2a) ou "série na carga" (Fig. 5.2b)
    x: float                # reatância série [Ω]
    b: float                # susceptância em derivação [S]
    serie: tuple = None     # ('L'|'C', valor)
    derivacao: tuple = None

    def abcd(self, f):
        zs = _impedancia_comp(self.serie, f)
        yp = 1 / _impedancia_comp(self.derivacao, f)
        if self.topologia.startswith("derivação"):
            return abcd.cascata(abcd.series(zs), abcd.shunt(yp))     # fonte - X - B - carga
        return abcd.cascata(abcd.shunt(yp), abcd.series(zs))         # fonte - B - X - carga


def rede_l(zl, z0, f0):
    """Todas as soluções reais da rede L (eqs. 5.3a-b e 5.6a-b)."""
    rl, xl = zl.real, zl.imag
    if rl <= 0:
        raise ValueError("A carga precisa ter parte real positiva.")
    sol = []
    disc = rl**2 + xl**2 - z0 * rl
    if disc >= 0:                                    # B junto à carga (zL fora do círculo 1 + jx... ou dentro)
        for s in (1, -1):
            b = (xl + s * np.sqrt(rl / z0) * np.sqrt(disc)) / (rl**2 + xl**2)
            if abs(b) < 1e-15:
                continue
            x = 1 / b + xl * z0 / rl - z0 / (b * rl)
            if abs(x) > 1e-12:
                sol.append(SolucaoL("derivação na carga", x, b, _comp_reativo(x, f0, True),
                                    _comp_reativo(b, f0, False)))
    if rl < z0:                                      # X junto à carga
        for s in (1, -1):
            x = s * np.sqrt(rl * (z0 - rl)) - xl
            b = s * np.sqrt((z0 - rl) / rl) / z0
            if abs(x) > 1e-12 and abs(b) > 1e-15:
                sol.append(SolucaoL("série na carga", x, b, _comp_reativo(x, f0, True),
                                    _comp_reativo(b, f0, False)))
    if not sol:
        raise ValueError("Sem solução de rede L com dois elementos para esta carga.")
    return sol


@dataclass
class SolucaoTres:
    """Rede de três elementos: Pi (derivação, série, derivação) ou T (série, derivação, série),
    na ordem fonte → carga; valores em Ω (X série) ou S (B derivação)."""
    topologia: str                # "Pi" ou "T"
    valores: tuple                # (B1, X2, B3) no Pi; (X1, B2, X3) no T
    comps: tuple                  # elementos ('L'|'C', valor) ou None se o valor se anula
    rv: float
    q: float

    def abcd(self, f):
        ms = []
        for i, comp in enumerate(self.comps):
            if comp is None:
                continue
            z = _impedancia_comp(comp, f)
            serie = (i == 1) if self.topologia == "Pi" else (i != 1)
            ms.append(abcd.series(z) if serie else abcd.shunt(1 / z))
        return abcd.cascata(*ms) if len(ms) > 1 else ms[0]

    def descricao(self):
        nomes = ("derivação", "série", "derivação") if self.topologia == "Pi" else ("série", "derivação", "série")
        return ", ".join(f"{n}: {fmt_comp(c)}" if c else f"{n}: —" for n, c in zip(nomes, self.comps))


def q_minimo(zl, z0):
    """Q mínimo de uma rede de três elementos (o da rede L): √(Rmaior/Rmenor - 1)."""
    rl = zl.real
    return float(np.sqrt(max(rl, z0) / min(rl, z0) - 1))


def _tres(zl, z0, f0, q, pi):
    rl = zl.real
    if rl <= 0:
        raise ValueError("A carga precisa ter parte real positiva.")
    qmin = q_minimo(zl, z0)
    if q <= qmin + 1e-9:
        raise ValueError(f"Q precisa ser maior que o da rede L ({qmin:.3f}).")
    rv = max(rl, z0) / (q**2 + 1) if pi else min(rl, z0) * (q**2 + 1)
    # lado da carga: zl → RV; lado da fonte: RV → z0 (o elemento do meio é a soma dos dois)
    lado_carga = [s for s in rede_l(zl, rv, f0) if s.topologia.startswith("derivação" if pi else "série")]
    lado_fonte = [s for s in rede_l(complex(rv), z0, f0) if s.topologia.startswith("série" if pi else "derivação")]
    sol = []
    for a in lado_fonte:
        for b in lado_carga:
            if pi:
                v = (a.b, a.x + b.x, b.b)
                comps = (_comp_reativo(v[0], f0, False), _comp_reativo(v[1], f0, True) if abs(v[1]) > 1e-9 else None,
                         _comp_reativo(v[2], f0, False))
            else:
                v = (a.x, a.b + b.b, b.x)
                comps = (_comp_reativo(v[0], f0, True), _comp_reativo(v[1], f0, False) if abs(v[1]) > 1e-12 else None,
                         _comp_reativo(v[2], f0, True))
            sol.append(SolucaoTres("Pi" if pi else "T", v, comps, rv, q))
    if not sol:
        raise ValueError("Sem solução de três elementos para esta carga e este Q.")
    return sol


def rede_pi(zl, z0, f0, q):
    """Redes Pi (até 4 soluções) com Q = √(Rmaior/RV - 1), RV = Rmaior/(Q² + 1) (Steer, eq. 8.72)."""
    return _tres(zl, z0, f0, q, True)


def rede_t(zl, z0, f0, q):
    """Redes T (até 4 soluções) com RV = Rmenor·(Q² + 1) > Rmaior."""
    return _tres(zl, z0, f0, q, False)


# ---------------------------------------------------------------------------
# Stubs
# ---------------------------------------------------------------------------
def comprimento_stub(reativo, z0, curto, serie):
    """Comprimento [λ] de um stub que realiza X (série) ou B (derivação), em (0, 0,5)."""
    u = reativo / z0 if serie else reativo * z0      # x = X/Z0 ou b = B·Z0
    # série: aberto x = -cot βl, curto x = tan βl; derivação: aberto b = tan βl, curto b = -cot βl
    tangente = curto if serie else not curto
    bl = np.arctan(u) if tangente else np.pi / 2 + np.arctan(u)
    return float(np.mod(bl, np.pi) / (2 * np.pi))


@dataclass
class SolucaoStub:
    d: float          # distância da carga ao stub [λ]
    l: float          # comprimento do stub [λ]
    reativo: float    # susceptância (derivação) ou reatância (série) do stub


def _raizes_periodicas(g, n=4000):
    """Raízes de g(d) em d ∈ [0, 0,5) por varredura e bisseção."""
    d = np.linspace(0, 0.5, n + 1)
    y = g(d)
    out = []
    for i in range(n):
        if np.isfinite(y[i]) and np.isfinite(y[i + 1]) and y[i] * y[i + 1] < 0:
            a, b = d[i], d[i + 1]
            for _ in range(60):
                m = 0.5 * (a + b)
                if g(np.array([a]))[0] * g(np.array([m]))[0] <= 0:
                    b = m
                else:
                    a = m
            out.append(0.5 * (a + b))
        elif y[i] == 0:
            out.append(d[i])
    return out


def stub_simples(zl, z0, serie=False, curto=False):
    """Soluções (d, l) do stub simples (Pozar, seção 5.2). d e l em λ, as duas mais curtas."""
    if zl.real <= 0:
        raise ValueError("A carga precisa ter parte real positiva.")
    if serie:
        g = lambda d: (_transformar(zl, z0, 2 * np.pi * d).real - z0) / z0
    else:
        g = lambda d: (1 / _transformar(zl, z0, 2 * np.pi * d)).real * z0 - 1
    sol = []
    for d in _raizes_periodicas(g):
        z = _transformar(zl, z0, 2 * np.pi * d)
        reat = -z.imag if serie else -(1 / z).imag
        sol.append(SolucaoStub(float(d), comprimento_stub(reat, z0, curto, serie), float(reat)))
    if not sol:
        raise ValueError("Não foi encontrada solução de stub simples.")
    return sol[:2]


def abcd_stub_ideal(f, f0, z0, l_lambda, curto, serie):
    """Stub TEM ideal de comprimento l (em λ em f0), em derivação ou em série."""
    bl = 2 * np.pi * l_lambda * np.asarray(f, dtype=float) / f0
    if serie:
        z = 1j * z0 * np.tan(bl) if curto else -1j * z0 / np.tan(bl)
        return abcd.series(z)
    y = -1j / (z0 * np.tan(bl)) if curto else 1j * np.tan(bl) / z0
    return abcd.shunt(y)


def abcd_stub_simples(f, f0, z0, s: SolucaoStub, curto, serie):
    """Fonte - stub - linha d - carga (linhas TEM ideais)."""
    linha_d = abcd.linha(z0, 2 * np.pi * s.d * np.asarray(f, dtype=float) / f0)
    return abcd.cascata(abcd_stub_ideal(f, f0, z0, s.l, curto, serie), linha_d)


@dataclass
class SolucaoStubDuplo:
    b1: float         # susceptância normalizada do stub 1 (junto à carga)
    b2: float
    l1: float         # [λ]
    l2: float


def stub_duplo(zl, z0, d_lambda=0.125, curto=False, d0_lambda=0.0):
    """Stub duplo em derivação (Pozar, eq. 5.22): stub 1 a d0 da carga, stub 2 a d do stub 1."""
    y_l = z0 / _transformar(zl, z0, 2 * np.pi * d0_lambda)        # normalizada no plano do stub 1
    gl, bl = y_l.real, y_l.imag
    t = np.tan(2 * np.pi * d_lambda)
    if abs(t) < 1e-9:
        raise ValueError("Espaçamento entre stubs de 0 ou λ/2 não permite casamento.")
    lim = (1 + t**2) / t**2
    if gl > lim:
        raise ValueError(f"Carga na região proibida (g = {gl:.3f} > {lim:.3f}): altere d0 (distância "
                         "da carga ao 1º stub) ou o espaçamento.")
    raiz = np.sqrt(max((1 + t**2) * gl - gl**2 * t**2, 0.0))
    sol = []
    for s in (1, -1):
        b1 = -bl + (1 + s * raiz) / t
        y1 = gl + 1j * (bl + b1)
        y2 = 1 / _transformar(1 / y1, 1.0, 2 * np.pi * d_lambda)
        b2 = -y2.imag
        sol.append(SolucaoStubDuplo(float(b1), float(b2), comprimento_stub(b1 / z0, z0, curto, False),
                                    comprimento_stub(b2 / z0, z0, curto, False)))
    return sol


def abcd_stub_duplo(f, f0, z0, s: SolucaoStubDuplo, d_lambda, curto, d0_lambda=0.0):
    f = np.asarray(f, dtype=float)
    return abcd.cascata(abcd_stub_ideal(f, f0, z0, s.l2, curto, False),
                        abcd.linha(z0, 2 * np.pi * d_lambda * f / f0),
                        abcd_stub_ideal(f, f0, z0, s.l1, curto, False),
                        abcd.linha(z0, 2 * np.pi * d0_lambda * f / f0))


# ---------------------------------------------------------------------------
# Transformadores multisseção (λ/4 em f0)
# ---------------------------------------------------------------------------
def binomial(z0, zl, n):
    """Z1..ZN do transformador binomial (eq. 5.53, ln Z(n+1) = ln Zn + 2^-N·C(N,n)·ln(ZL/Z0))."""
    if n < 1:
        raise ValueError("Use ao menos 1 seção.")
    z, out = np.log(z0), []
    for k in range(n):
        z += comb(n, k) * np.log(zl / z0) / 2**n
        out.append(float(np.exp(z)))
    return out


def banda_binomial(z0, zl, n, gama_m):
    """Δf/f0 do binomial para |Γ| ≤ Γm (eq. 5.55, com A ≈ 2^-(N+1)·ln(ZL/Z0))."""
    a = abs(np.log(zl / z0)) / 2 ** (n + 1)
    x = 0.5 * (gama_m / a) ** (1 / n)
    if x >= 1:
        return 2.0
    return float(2 - 4 / np.pi * np.arccos(x))


def chebyshev(z0, zl, n, gama_m):
    """Z1..ZN do transformador de Chebyshev (eqs. 5.61-5.63) e θm [rad].

    Os Γn saem da expansão de A·T_N(sec θm·cos θ) em cos((N-2n)θ), com A = Γm·sinal(ln ZL/Z0).
    """
    if n < 1:
        raise ValueError("Use ao menos 1 seção.")
    r = abs(np.log(zl / z0)) / (2 * gama_m)
    if r <= 1:
        raise ValueError("Γm maior que o descasamento: não é preciso transformador.")
    sec_m = np.cosh(np.arccosh(r) / n)
    a = gama_m * np.sign(np.log(zl / z0))
    th = np.linspace(0, np.pi, 4097)
    x = sec_m * np.cos(th)
    tn = np.where(np.abs(x) <= 1, np.cos(n * np.arccos(np.clip(x, -1, 1))),
                  np.sign(x) ** n * np.cosh(n * np.arccosh(np.maximum(np.abs(x), 1))))
    func = a * tn
    gam = []
    for k in range(n + 1):
        m = n - 2 * k
        if m < 0:
            break
        coef = np.trapezoid(func * np.cos(m * th), th) / np.pi * (2 if m > 0 else 1)
        gam.append(coef / 2 if m > 0 else coef)
    gam = gam + gam[: n + 1 - len(gam)][::-1]                    # simetria Γn = Γ(N-n)
    z, out = np.log(z0), []
    for k in range(n):
        z += 2 * gam[k]
        out.append(float(np.exp(z)))
    return out, float(np.arccos(1 / sec_m))


def banda_chebyshev(theta_m):
    return float(2 - 4 * theta_m / np.pi)                        # eq. 5.64


def abcd_multissecao(f, f0, zs):
    """Cascata de seções λ/4 ideais em f0 (da fonte para a carga)."""
    th = np.pi / 2 * np.asarray(f, dtype=float) / f0
    return abcd.cascata(*[abcd.linha(z, th) for z in zs])


# ---------------------------------------------------------------------------
# Tapers
# ---------------------------------------------------------------------------
TAPERS = ("Exponencial", "Triangular", "Klopfenstein")


def _i1_sobre_x(u):
    """I1(u)/u pela série de potências (estável para todo u ≥ 0 usado aqui, u ≤ ~15)."""
    u = np.asarray(u, dtype=float)
    q = (u / 2) ** 2
    termo = np.full_like(u, 0.5)
    soma = termo.copy()
    for k in range(1, 80):
        termo = termo * q / (k * (k + 1))
        soma = soma + termo
    return soma


def phi_klopfenstein(x, a, n=2001):
    """φ(x, A) = ∫0^x I1(A√(1-y²))/(A√(1-y²)) dy (eq. 5.75), por trapézios; x em [-1, 1]."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    out = np.empty_like(x)
    for i, xi in enumerate(x):
        y = np.linspace(0, abs(xi), n)
        out[i] = np.sign(xi) * np.trapezoid(_i1_sobre_x(a * np.sqrt(np.clip(1 - y**2, 0, None))), y)
    return out


def perfil_taper(tipo, z0, zl, n_pts=401, gama_m=0.02):
    """(z/L, Z(z)) do taper; Klopfenstein com ripple Γm na banda passante (βL ≥ A)."""
    s = np.linspace(0, 1, n_pts)
    ln = np.log(zl / z0)
    if tipo == "Exponencial":                                         # eq. 5.68-5.69
        return s, z0 * np.exp(s * ln)
    if tipo == "Triangular":                                          # eq. 5.71
        e = np.where(s <= 0.5, 2 * s**2, 4 * s - 2 * s**2 - 1)
        return s, z0 * np.exp(e * ln)
    if tipo == "Klopfenstein":                                        # eqs. 5.74-5.77
        g0 = 0.5 * ln
        if gama_m <= 0 or gama_m >= abs(g0):
            raise ValueError("Γm deve estar entre 0 e |Γ0| = |ln(ZL/Z0)|/2.")
        a = np.arccosh(abs(g0) / gama_m)
        lnz = 0.5 * np.log(z0 * zl) + g0 / np.cosh(a) * a**2 * phi_klopfenstein(2 * s - 1, a)
        return s, np.exp(lnz)
    raise ValueError(f"Taper desconhecido: {tipo}")


def parametro_a(z0, zl, gama_m):
    return float(np.arccosh(abs(0.5 * np.log(zl / z0)) / gama_m))


def gama_taper_analitico(tipo, bl, z0, zl, gama_m=0.02):
    """|Γ(βL)| pela teoria das pequenas reflexões (eqs. 5.70, 5.73, 5.76)."""
    bl = np.asarray(bl, dtype=float)
    ln = abs(np.log(zl / z0))
    if tipo == "Exponencial":
        return 0.5 * ln * np.abs(np.sinc(bl / np.pi))
    if tipo == "Triangular":
        return 0.5 * ln * np.abs(np.sinc(bl / (2 * np.pi))) ** 2
    a = parametro_a(z0, zl, gama_m)
    g0 = 0.5 * ln
    arg = bl**2 - a**2
    v = np.where(arg >= 0, np.cos(np.sqrt(np.abs(arg))), np.cosh(np.sqrt(np.abs(arg))))
    return g0 * np.abs(v) / np.cosh(a)


def abcd_perfil(f, comprimentos_m, zs, eeffs, alfas=None):
    """Cascata de trechos curtos (Z, εeff, α por trecho): taper ou degraus em linha real."""
    f = np.asarray(f, dtype=float)
    k0 = 2 * np.pi * f / 299_792_458.0
    m = None
    for i, (dl, z, ee) in enumerate(zip(comprimentos_m, zs, eeffs)):
        th = k0 * np.sqrt(ee) * dl
        if alfas is not None:
            th = th - 1j * alfas[i] * dl
        mi = abcd.linha(z, th)
        m = mi if m is None else m @ mi
    return m


# ---------------------------------------------------------------------------
# Bode-Fano
# ---------------------------------------------------------------------------
@dataclass
class BodeFano:
    tipo: str
    gama_min: float
    rl_max_db: float
    texto: list = field(default_factory=list)


def bode_fano(zl0, f0, f1, f2):
    """Menor |Γ| possível, constante em [f1, f2] e 1 fora dela, para R série com L ou C (Fig. 5.23)."""
    r, x = zl0.real, zl0.imag
    w0, w1, w2 = (2 * np.pi * v for v in (f0, f1, f2))
    if not (0 < w1 < w2):
        raise ValueError("Banda inválida para Bode-Fano.")
    if x < 0:                                                     # RC série: ∫ω⁻²·ln(1/|Γ|)dω ≤ πRC
        c = -1 / (w0 * x)
        k = np.pi * r * c
        integral = 1 / w1 - 1 / w2
        tipo = "RC série: ∫ ω⁻² ln(1/|Γ|) dω ≤ π·R·C"
    elif x > 0:                                                   # RL série: ∫ln(1/|Γ|)dω ≤ πR/L
        l = x / w0
        k = np.pi * r / l
        integral = w2 - w1
        tipo = "RL série: ∫ ln(1/|Γ|) dω ≤ π·R/L"
    else:
        return BodeFano("resistiva", 0.0, np.inf, ["Carga resistiva: casamento perfeito em banda larga é "
                                                   "possível (ex.: transformadores)."])
    g = float(np.exp(-k / integral))
    rl = -20 * np.log10(g) if g > 0 else np.inf
    return BodeFano(tipo, g, rl, [f"Bode-Fano ({tipo}):",
                                  f"  em {f1 / 1e9:g}-{f2 / 1e9:g} GHz, |Γ| ≥ {g:.4f} "
                                  f"(perda de retorno ≤ {rl:.2f} dB) para qualquer rede sem perdas."])
