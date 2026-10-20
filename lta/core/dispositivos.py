"""Dispositivos de N portas montados com linhas (convencionais ou artificiais).

Referências: F. Martín (2015), seções 4.2.1 e 4.2.2 (híbridas, divisores, banda) e 2.5.2
(linha com carga capacitiva, slow-wave); D. M. Pozar, Microwave Engineering, 4ª ed., 2012,
seções 7.2 (junção T), 7.3 (Wilkinson, inclusive desigual, eq. 7.37), 7.5 (branch-line) e
7.8 (anel híbrido).

Rede: nós ligados por bipolos descritos por ABCD (Nf, 2, 2), todos com o terra
comum. A matriz Y nodal é montada, os nós internos são eliminados (redução de
Kron) e S sai das ondas de potência com referência real Rk em cada porta:
S = R^-1/2 (I - R·Y)(I + R·Y)^-1 R^1/2 (igual a (Y0 - Y)(Y0 + Y)^-1 se todas forem Z0).
"""
from dataclasses import dataclass

import numpy as np

from . import abcd
from . import carga_capacitiva as cc
from .crlh import CelulaCRLH
from .microstrip import C0, Substrato, analise, parametros_linha, sintese
from .multibanda import crlh_dual


class Rede:
    def __init__(self, n_nos, portas, f):
        self.f = np.asarray(f, dtype=float)
        self.n = n_nos
        self.portas = list(portas)
        self.y = np.zeros((self.f.size, n_nos, n_nos), dtype=complex)

    def bipolo(self, a, b, m):
        """Liga o quadripolo de matriz ABCD m entre os nós a (porta 1) e b (porta 2)."""
        A, B, C, D = m[:, 0, 0], m[:, 0, 1], m[:, 1, 0], m[:, 1, 1]
        det = A * D - B * C
        B = np.where(np.abs(B) < 1e-12, 1e-12, B)          # bipolo em curto (ex.: ressonância série)
        self.y[:, a, a] += D / B
        self.y[:, a, b] += -det / B
        self.y[:, b, a] += -1 / B
        self.y[:, b, b] += A / B

    def multiporta(self, nos, s, z0=50.0):
        """Liga um bloco de N portas (S (Nf,N,N) com referência z0) aos nós dados: Y = (I - S)(I + S)⁻¹/z0."""
        n = len(nos)
        eye = np.eye(n)
        y = (eye - s) @ np.linalg.inv(eye + s) / z0
        for i, a in enumerate(nos):
            for j, b in enumerate(nos):
                self.y[:, a, b] += y[:, i, j]

    def impedancia(self, a, b, z):
        """Impedância z (escalar ou por frequência) entre os nós a e b."""
        y = np.broadcast_to(1 / np.asarray(z, dtype=complex), self.f.shape)
        self.y[:, a, a] += y
        self.y[:, b, b] += y
        self.y[:, a, b] -= y
        self.y[:, b, a] -= y

    def terra(self, a, z):
        """Impedância z (escalar ou por frequência) do nó a para o terra."""
        self.y[:, a, a] += np.broadcast_to(1 / np.asarray(z, dtype=complex), self.f.shape)

    def s(self, z0=50.0):
        """S com referência z0 (escalar) ou uma impedância real por porta (sequência).

        As portas são terminadas nas suas referências (G = 1/R) e o sistema nodal completo é
        invertido: S = 2·G^1/2·Zpp·G^1/2 - I. Evita a eliminação dos nós internos, que fica
        singular quando eles só veem linhas λ/4 ideais (exatamente em f0).
        """
        p = self.portas
        r = np.broadcast_to(np.asarray(z0, dtype=float), (len(p),))
        yt = self.y.copy()
        yt[:, p, p] += 1 / r
        zpp = np.linalg.inv(yt)[:, p][:, :, p]
        g = 1 / np.sqrt(r)
        return 2 * g[:, None] * zpp * g[None, :] - np.eye(len(p))


# ---------------------------------------------------------------------------
# Linhas de λ/4 em diferentes implementações
# ---------------------------------------------------------------------------
def linha_convencional(z, fase0, f0):
    """Linha TEM ideal de impedância z e comprimento elétrico fase0 [rad] em f0."""
    return lambda f: abcd.linha(z, fase0 * np.asarray(f, dtype=float) / f0)


def linha_concentrada(z, fase0, f0):
    """Equivalente π concentrado de uma linha z, fase0 [rad] em f0 (Steer, seção 7.13): ramo série
    X = Z·sen θ e derivações B = (1 - cos θ)/(Z·sen θ). Com θ = 90° é passa-baixas (L = Z/ω0,
    C = 1/(ω0·Z)); 270° é realizado como -90° passa-altas (C série, L em derivação).
    Retorna (f -> ABCD, (('L'|'C', série), ('L'|'C', derivação)))."""
    th = np.angle(np.exp(1j * fase0))                   # (-π, π]: 270° → -90°
    if abs(np.sin(th)) < 1e-6:
        raise ValueError("Fase sem equivalente π concentrado (0° ou 180°).")
    w0 = 2 * np.pi * f0
    x = z * np.sin(th)
    b = (1 - np.cos(th)) / (z * np.sin(th))
    serie = ("L", x / w0) if x > 0 else ("C", -1 / (w0 * x))
    deriv = ("C", b / w0) if b > 0 else ("L", -1 / (w0 * b))

    def imp(comp, w):
        return 1j * w * comp[1] if comp[0] == "L" else 1 / (1j * w * comp[1])

    def m(f):
        w = 2 * np.pi * np.asarray(f, dtype=float)
        y = 1 / imp(deriv, w)
        return abcd.cascata(abcd.shunt(y), abcd.series(imp(serie, w)), abcd.shunt(y))
    return m, (serie, deriv)


def linha_crlh_dual(z, f1, f2, fase1=-np.pi / 2, fase2=np.pi / 2):
    """Célula CRLH com fase1 em f1 e fase2 em f2 e impedância z (eqs. 4.28-4.32)."""
    c = crlh_dual(f1, fase1, z, f2, fase2, z)
    return c.abcd, c


def linha_microstrip(z, fase0, f0, sub: Substrato, perdas=False, dispersivo=False):
    """Linha microstrip real (Hammerstad-Jensen, com perdas/dispersão opcionais) de fase0 em f0.

    Retorna (função f -> ABCD, (W [mm], comprimento [mm], εeff)).
    """
    w = float(sintese(z, sub))
    _, ee0 = analise(w, sub)
    comp = fase0 * C0 / (2 * np.pi * f0 * np.sqrt(float(ee0)))

    def m(f):
        zf, ee, a = (x[:, 0] for x in parametros_linha(w, sub, f, perdas, dispersivo))
        th = 2 * np.pi * np.asarray(f, dtype=float) * np.sqrt(ee) / C0 * comp - 1j * a * comp
        return abcd.linha(zf, th)
    return m, (w, comp * 1e3, float(ee0))


CARGAS_SW = ("Stubs abertos", "Patches")


@dataclass
class GeometriaSW:
    """Célula slow-wave física: linha Zo (w_linha) com patch centrado ou stub aberto no centro."""
    tipo: str            # CARGAS_SW
    w_linha: float       # mm
    eeff_linha: float
    periodo: float       # mm
    w_carga: float       # mm: largura do patch ou do stub
    l_carga: float       # mm: comprimento do patch (ao longo da linha) ou do stub (transversal)
    eeff_carga: float
    z_carga: float       # Ω


@dataclass
class ProjetoSlowWave:
    """Linha com carga capacitiva periódica que substitui um trecho de fase0 com impedância z."""
    z: float
    fase0: float
    n: int
    f0: float
    carga: cc.ProjetoCarga
    geo: GeometriaSW
    comp_convencional: float     # mm, linha comum de mesma fase e impedância

    @property
    def comp(self):
        return self.n * self.geo.periodo

    @property
    def reducao(self):
        return self.comp / self.comp_convencional

    @property
    def zo(self):
        return self.carga.zo


def projeto_slow_wave(z, fase0, f0, n, zo, sub: Substrato, z_carga=20.0, tipo=CARGAS_SW[0]):
    """Célula distribuída meia linha (Zo) + Cls + meia linha com fase0/n e Bloch z exatos em f0.

    Com B = ω0·Cls, a célula simétrica tem cos φ = cos θ - (B·Zo/2)·sen θ e
    ZB² = Zo·(sen θ - Zo·B·(1 - cos θ)/2)/(sen θ/Zo + B·(1 + cos θ)/2); θ (fase da linha
    descarregada por célula) sai por bisseção em (0, φ). O fator slow-wave é ≈ Zo/z.
    Cls é realizada por um patch de impedância z_carga (capacitância excedente, como em
    carga_capacitiva.geometria) ou por um stub aberto de z_carga com tan(βℓ)/Z = ω0·Cls.
    """
    if n < 1:
        raise ValueError("Use ao menos 1 célula.")
    if zo <= z:
        raise ValueError(f"A linha não carregada (Zo = {zo:g} Ω) precisa ter impedância maior que a "
                         f"desejada ({z:.2f} Ω).")
    phi = abs(fase0) / n
    if phi >= np.pi:
        raise ValueError("Fase por célula ≥ 180°: aumente o número de células.")

    def b_de(th):
        return 2 * (np.cos(th) - np.cos(phi)) / (zo * np.sin(th))

    def zb(th):
        b = b_de(th)
        num = np.sin(th) - zo * b * (1 - np.cos(th)) / 2
        den = np.sin(th) / zo + b * (1 + np.cos(th)) / 2
        return np.sqrt(max(zo * num / den, 0.0))

    lo, hi = 1e-6 * phi, phi * (1 - 1e-9)
    for _ in range(100):
        m = 0.5 * (lo + hi)
        if zb(m) < z:
            lo = m
        else:
            hi = m
    th = 0.5 * (lo + hi)
    w0 = 2 * np.pi * f0
    tau = th / w0
    carga = cc.ProjetoCarga(zo * tau, tau / zo, b_de(th) / w0)
    w_l = float(sintese(zo, sub))
    _, ee_l = analise(w_l, sub)
    periodo = C0 * tau / np.sqrt(float(ee_l)) * 1e3
    w_p = float(sintese(z_carga, sub))
    _, ee_p = analise(w_p, sub)
    ee_p = float(ee_p)
    if tipo == CARGAS_SW[1]:
        c_m = np.sqrt(ee_p) / (z_carga * C0) - np.sqrt(float(ee_l)) / (zo * C0)
        if c_m <= 0:
            raise ValueError("A impedância do patch deve ser menor que Zo.")
        l_c = carga.Cls / c_m * 1e3
        if l_c >= 0.9 * periodo:
            raise ValueError(f"O patch ({l_c:.2f} mm) não cabe no período ({periodo:.2f} mm): aumente o "
                             "número de células, reduza Zo, baixe a impedância do patch ou use stubs.")
    else:
        bl = np.arctan(w0 * carga.Cls * z_carga)
        l_c = bl * C0 / (w0 * np.sqrt(ee_p)) * 1e3
        if w_p >= 0.8 * periodo:
            raise ValueError(f"O stub ({w_p:.2f} mm de largura) não cabe no período ({periodo:.2f} mm): "
                             "aumente a impedância do stub ou use menos células.")
    geo = GeometriaSW(tipo, w_l, float(ee_l), float(periodo), w_p, float(l_c), ee_p, float(z_carga))
    w_c = float(sintese(z, sub))
    _, ee_c = analise(w_c, sub)
    comp_conv = abs(fase0) * C0 / (w0 * np.sqrt(float(ee_c))) * 1e3
    return ProjetoSlowWave(z, fase0, int(n), f0, carga, geo, float(comp_conv))


def linha_slow_wave(p: ProjetoSlowWave, sub: Substrato, perdas=False, dispersivo=False):
    """f -> ABCD das n células com a linha microstrip real (dispersão e perdas opcionais).

    Patch: Cls concentrada (como carga_capacitiva.celula_microstrip). Stub: Y = j·tan(βℓ)/Z do
    stub aberto real, que coincide com ω0·Cls em f0.
    """
    g = p.geo

    def m(f):
        f = np.asarray(f, dtype=float)
        z_l, ee, a = (x[:, 0] for x in parametros_linha(g.w_linha, sub, f, perdas, dispersivo))
        th = 2 * np.pi * f * np.sqrt(ee) / C0 * g.periodo * 1e-3 - 1j * a * g.periodo * 1e-3
        meia = abcd.linha(z_l, th / 2)
        if g.tipo == CARGAS_SW[1]:
            y = 1j * 2 * np.pi * f * p.carga.Cls
        else:
            z_s, ee_s, a_s = (x[:, 0] for x in parametros_linha(g.w_carga, sub, f, perdas, dispersivo))
            bl = 2 * np.pi * f * np.sqrt(ee_s) / C0 * g.l_carga * 1e-3 - 1j * a_s * g.l_carga * 1e-3
            y = np.tanh(1j * bl) / z_s
        return abcd.potencia(abcd.cascata(meia, abcd.shunt(y), meia), p.n)
    return m


def linha_plh(z, fase0, f0, n=1):
    """N células PLH (T: braços 2CL, derivação LL) com fase total fase0 < 0 em f0."""
    if fase0 >= 0:
        raise ValueError("A linha PLH tem fase negativa.")
    w = 2 * np.pi * f0
    fc = fase0 / n
    cl = -1 / (2 * w * z * np.tan(fc / 2))      # braço -1/(2ωCL) = Z·tan(φ/2)
    ll = -z / (w * np.sin(fc))                  # ωLL = -Z/sen φ
    return lambda f: abcd.potencia(abcd.celula(f, "T", cs=cl, lp=ll), n), (cl, ll)


# ---------------------------------------------------------------------------
# Dispositivos (portas numeradas a partir de 0 internamente; 1..N no texto)
# ---------------------------------------------------------------------------
def branch_line(f, l_serie, l_derivacao, z0=50.0):
    """Híbrida de quadratura: 1-2 e 4-3 com l_serie (Z0/√2), 1-4 e 2-3 com l_derivacao (Z0).

    Entrada 1, saída direta 2, acoplada 3, isolada 4 (livro, Fig. 4.16).
    """
    r = Rede(4, [0, 1, 2, 3], f)
    r.bipolo(0, 1, l_serie(f))
    r.bipolo(3, 2, l_serie(f))
    r.bipolo(0, 3, l_derivacao(f))
    r.bipolo(1, 2, l_derivacao(f))
    return r.s(z0)


def rat_race(f, l90, l270, z0=50.0):
    """Anel 1-2 (90°), 2-4 (90°), 4-3 (90°), 3-1 (270° ou -90° LH), todos Z0·√2.

    Porta 4 = Σ (saídas 2 e 3 em fase), porta 1 = Δ (saídas em oposição) (livro, Fig. 4.7).
    """
    r = Rede(4, [0, 1, 2, 3], f)
    r.bipolo(0, 1, l90(f))
    r.bipolo(1, 3, l90(f))
    r.bipolo(3, 2, l90(f))
    r.bipolo(2, 0, l270(f))
    return r.s(z0)


def wilkinson(f, l_braco, z0=50.0):
    """Divisor de Wilkinson: dois braços (Z0·√2, 90°) e resistor 2·Z0 entre as saídas."""
    r = Rede(3, [0, 1, 2], f)
    r.bipolo(0, 1, l_braco(f))
    r.bipolo(0, 2, l_braco(f))
    r.impedancia(1, 2, 2 * z0)
    return r.s(z0)


def wilkinson_desigual(f, l2, l3, k2, z0=50.0, t2=None, t3=None):
    """Wilkinson com P3/P2 = k2 (eq. 7.37): braços Z02 e Z03, R = Z0(K + 1/K).

    Sem transformadores (t2, t3 = None), as portas 2 e 3 têm referência R2 = Z0·K e R3 = Z0/K;
    com eles (λ/4 de √(Z0·R2) e √(Z0·R3)), todas as portas ficam em Z0.
    """
    k = np.sqrt(k2)
    r = Rede(5 if t2 else 3, [0, 3, 4] if t2 else [0, 1, 2], f)
    r.bipolo(0, 1, l2(f))
    r.bipolo(0, 2, l3(f))
    r.impedancia(1, 2, z0 * (k + 1 / k))
    if t2:
        r.bipolo(1, 3, t2(f))
        r.bipolo(2, 4, t3(f))
        return r.s(z0)
    return r.s([z0, z0 * k, z0 / k])


def impedancias_wilkinson(z0, k2):
    """(Z02, Z03, R, R2, R3) do Wilkinson desigual (eq. 7.37)."""
    k = np.sqrt(k2)
    z03 = z0 * np.sqrt((1 + k2) / k**3)
    return k2 * z03, z03, z0 * (k + 1 / k), z0 * k, z0 / k


def impedancias_juncao_t(z0, p3_p2):
    """(Z2, Z3) das saídas da junção T sem perdas casada na entrada (Pozar, Ex. 7.1)."""
    return z0 * (1 + p3_p2), z0 * (1 + p3_p2) / p3_p2


def juncao_t(f, z0, z2, z3, t2=None, t3=None):
    """Junção T sem perdas: entrada Z0 e saídas Z2 e Z3 no mesmo nó (portas com essas referências)
    ou, com transformadores λ/4 para Z0 (t2, t3), todas as portas em Z0."""
    r = Rede(3, [0, 1, 2], f)
    if t2:
        r.bipolo(0, 1, t2(f))
        r.bipolo(0, 2, t3(f))
        return r.s(z0)
    r.impedancia(0, 1, 1e-9)                    # nó único (ligações de 1 nΩ)
    r.impedancia(0, 2, 1e-9)
    return r.s([z0, z2, z3])


def divisor_inversor(f, l_inversor, z0=50.0):
    """Divisor em Y: inversor (Z0/√2, 90°) da porta 1 a um nó com as portas 2 e 3 (Figs. 4.11 e 4.15)."""
    r = Rede(3, [0, 1, 2], f)
    r.bipolo(0, 1, l_inversor(f))
    r.impedancia(1, 2, 1e-9)                    # portas 2 e 3 no mesmo nó (ligação de 1 nΩ)
    return r.s(z0)


# ---------------------------------------------------------------------------
# Limites de banda (seção 4.2.1)
# ---------------------------------------------------------------------------
def inclinacao_dual(fase0, w0, n=1):
    """dφ/dω da linha dual (PLH) de N células com fase total fase0 < 0 (eq. 4.9).

    Por célula: cos φ = 1 - 1/(2ω²LLCL) ⇒ dφ/dω = (2/ω)·tan(|φ|/2); N → ∞ dá |φ0|/ω0 (eq. 4.11).
    """
    return 2 * n / w0 * np.tan(abs(fase0) / (2 * n))


def fase_limite_uma_celula():
    """|φ| acima do qual 1 célula dual tem menor dφ/dω que a linha convencional de 2π - |φ|.

    Raiz de 2·tan(x/2) = 2π - x; o livro dá x ≈ 0,7π (equivalente a 1,3π convencional).
    """
    a, b = 0.1, np.pi - 1e-9
    for _ in range(200):
        m = 0.5 * (a + b)
        if 2 * np.tan(m / 2) - (2 * np.pi - m) > 0:
            b = m
        else:
            a = m
    return 0.5 * (a + b)


def atraso_crlh_zero(c: CelulaCRLH, n=1):
    """Atraso de grupo de N células CRLH balanceadas na transição: τ = 2N/ωR (eq. 4.15)."""
    return 2 * n / c.wR


# ---------------------------------------------------------------------------
# Métricas
# ---------------------------------------------------------------------------
def faixa_em_torno(f, ok, f0):
    """Faixa contígua de f em que ok é verdadeiro e que contém f0: (f_ini, f_fim) ou None."""
    i = int(np.argmin(np.abs(f - f0)))
    if not ok[i]:
        return None
    a = i
    while a > 0 and ok[a - 1]:
        a -= 1
    b = i
    while b < len(f) - 1 and ok[b + 1]:
        b += 1
    return f[a], f[b]


def diferenca_fase(s_a, s_b):
    """∠s_a - ∠s_b em graus, em (-180°, 180°]."""
    return np.degrees(np.angle(s_a * np.conj(s_b)))
