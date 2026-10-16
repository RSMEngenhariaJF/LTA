"""Dispositivos de N portas montados com linhas (convencionais ou artificiais).

Referência: F. Martín (2015), seções 4.2.1 e 4.2.2 (híbridas, divisores, banda).

Rede: nós ligados por bipolos descritos por ABCD (Nf, 2, 2), todos com o terra
comum. A matriz Y nodal é montada, os nós internos são eliminados (redução de
Kron) e S sai de S = (Y0 - Y)(Y0 + Y)^-1 com a mesma Z0 em todas as portas.
"""
import numpy as np

from . import abcd
from .crlh import CelulaCRLH
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
        self.y[:, a, a] += D / B
        self.y[:, a, b] += -det / B
        self.y[:, b, a] += -1 / B
        self.y[:, b, b] += A / B

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
        p = self.portas
        q = [i for i in range(self.n) if i not in p]
        y = self.y[:, p][:, :, p]
        if q:
            ypq = self.y[:, p][:, :, q]
            yqp = self.y[:, q][:, :, p]
            yqq = self.y[:, q][:, :, q]
            y = y - ypq @ np.linalg.solve(yqq, yqp)
        eye = np.eye(len(p)) / z0
        return (eye - y) @ np.linalg.inv(eye + y)


# ---------------------------------------------------------------------------
# Linhas de λ/4 em diferentes implementações
# ---------------------------------------------------------------------------
def linha_convencional(z, fase0, f0):
    """Linha TEM ideal de impedância z e comprimento elétrico fase0 [rad] em f0."""
    return lambda f: abcd.linha(z, fase0 * np.asarray(f, dtype=float) / f0)


def linha_crlh_dual(z, f1, f2, fase1=-np.pi / 2, fase2=np.pi / 2):
    """Célula CRLH com fase1 em f1 e fase2 em f2 e impedância z (eqs. 4.28-4.32)."""
    c = crlh_dual(f1, fase1, z, f2, fase2, z)
    return c.abcd, c


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
