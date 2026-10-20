"""Linhas acopladas simétricas: modos par/ímpar, acopladores, Lange e versão artificial (LTA).

Referências:
  - Pozar, Microwave Engineering, 4ª ed., 2012: seções 7.6 (acoplador de linhas acopladas,
    eqs. 7.81-7.91, Ex. 7.7-7.8) e 7.7 (Lange, eqs. 7.97-7.100);
  - microstrip acoplada: Garg e Bahl, "Characteristics of coupled microstriplines", IEEE
    Trans. MTT-27(7), 1979, pp. 700-705 (capacitâncias par/ímpar, precisão ~3% para
    0,2 ≤ W/h ≤ 2, 0,05 ≤ S/h ≤ 2);
  - microstrip acoplada (padrão): Kirschning e Jansen, IEEE Trans. MTT-32(1), 1984, com
    espessura, dispersão, perdas de Wheeler e ponta aberta por modo, conforme R. Garg, I. Bahl e
    M. Bianchi, Microstrip Lines and Slotlines, 4ª ed., 2024, eqs. 4.79-4.113 e 4.155;
  - stripline acoplada: S. B. Cohn, "Shielded coupled-strip transmission line", IRE Trans.
    MTT-3, 1955 (exata para t = 0; Pozar, Fig. 7.29);
  - LTA: as microstrips acopladas têm εeff do modo par maior que o do ímpar (Pozar, p. 375),
    o que degrada a diretividade e cria a banda espúria em 2f0 dos filtros de linhas
    acopladas (Martín 2015, seção 4.2.3.2). Carregar periodicamente o modo ímpar com
    capacitores entre as linhas (o modo par não os vê) iguala as velocidades de fase:
    é uma linha artificial com carga capacitiva (Martín, seção 2.5.2) aplicada a um só modo.

Convenção das 4 portas (Pozar, Fig. 7.31): 1 = entrada (linha a, esquerda), 2 = direta
(linha a, direita), 3 = acoplada (linha b, esquerda), 4 = isolada (linha b, direita).
"""
from dataclasses import dataclass, field

import numpy as np

from . import abcd
from .csrr import elipk
from .microstrip import C0, ETA0, Substrato, alfa_condutor, analise, resistencia_superficial
from .linhas_planares import alfa_dieletrico

EPS0 = 8.8541878128e-12
TIPOS = ("Microstrip", "Stripline")


def _razao(k):
    """K(k)/K(k')."""
    return elipk(k) / elipk(np.sqrt(1 - k * k))


# ---------------------------------------------------------------------------
# Análise
# ---------------------------------------------------------------------------
MODELOS = ("Kirschning-Jansen", "Garg-Bahl")


@dataclass
class Modos:
    ze: float
    zo: float
    eee: float     # εeff do modo par
    eeo: float     # εeff do modo ímpar
    disp: object = None   # f -> (εee(f), εeo(f)), se a dispersão estiver ativa

    def eeff(self, f):
        """(εee, εeo) na(s) frequência(s) f (estáticos sem dispersão)."""
        if self.disp is None:
            return self.eee, self.eeo
        return self.disp(f)

    @property
    def z0(self):
        return float(np.sqrt(self.ze * self.zo))

    @property
    def c(self):
        """Acoplamento de tensão em f0, (Ze - Zo)/(Ze + Zo) (eq. 7.81)."""
        return (self.ze - self.zo) / (self.ze + self.zo)


def _garg_bahl_c(w, s, h, er):
    """(Ce, Co) [F/m] das microstrips acopladas (Garg-Bahl). w, s, h em mm."""
    u, g = w / h, s / h
    sub = Substrato(er, h)
    z1, e1 = analise(w, sub)
    cp = EPS0 * er * u
    cf = 0.5 * (np.sqrt(e1) / (C0 * z1) - cp)
    a = np.exp(-0.1 * np.exp(2.33 - 2.53 * u))
    cf_l = cf / (1 + a / g * np.tanh(8 * g)) * (er / e1) ** 0.25
    k = g / (g + 2 * u)
    cga = EPS0 / _razao(k)                                    # ε0·K(k')/K(k)
    cgd = EPS0 * er / np.pi * np.log(1 / np.tanh(np.pi * g / 4)) + \
        0.65 * cf * (0.02 * np.sqrt(er) / g + 1 - er ** -2)
    return cp + cf + cf_l, cp + cf + cga + cgd


def microstrip_acoplada(w, s, sub: Substrato):
    """Modos de duas microstrips de largura w separadas por s [mm] (Garg-Bahl, quase-estático)."""
    if w <= 0 or s <= 0:
        raise ValueError("W e S devem ser positivos.")
    ce, co = _garg_bahl_c(w, s, sub.h_mm, sub.er)
    ce_a, co_a = _garg_bahl_c(w, s, sub.h_mm, 1.0)
    return Modos(float(1 / (C0 * np.sqrt(ce * ce_a))), float(1 / (C0 * np.sqrt(co * co_a))),
                 float(ce / ce_a), float(co / co_a))


def _kj_estatico(u, g, er):
    """(εee(0), εeo(0), Ze/Z0 e Zo/Z0 sem a espessura) de Kirschning-Jansen (eqs. 4.96-4.101)."""
    sub = Substrato(er, 1.0)
    z0, ee0 = (float(x) for x in analise(u, sub))
    v = u * (20 + g**2) / (10 + g**2) + g * np.exp(-g)
    ae = 1 + np.log((v**4 + (v / 52) ** 2) / (v**4 + 0.432)) / 49 + np.log(1 + (v / 18.1) ** 3) / 18.7
    be = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053
    eee = (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / v) ** (-ae * be)
    ao = 0.7287 * (ee0 - 0.5 * (er + 1)) * (1 - np.exp(-0.179 * u))
    bo = 0.747 * er / (0.15 + er)
    co = bo - (bo - 0.207) * np.exp(-0.414 * u)
    do = 0.593 + 0.694 * np.exp(-0.562 * u)
    eeo = ee0 + (0.5 * (er + 1) - ee0 + ao) * np.exp(-co * g**do)
    q1 = 0.8695 * u**0.194
    q2 = 1 + 0.7519 * g + 0.189 * g**2.31
    q3 = 0.1975 + (16.6 + (8.4 / g) ** 6) ** -0.387 + np.log(g**10 / (1 + (g / 3.4) ** 10)) / 241
    q4 = 2 * q1 / q2 / (u**q3 * np.exp(-g) + (2 - np.exp(-g)) * u**-q3)
    q5 = 1.794 + 1.14 * np.log(1 + 0.638 / (g + 0.517 * g**2.43))
    q6 = 0.2305 + np.log(g**10 / (1 + (g / 5.8) ** 10)) / 281.3 + np.log(1 + 0.598 * g**1.154) / 5.1
    q7 = (10 + 190 * g**2) / (1 + 82.3 * g**3)
    q8 = np.exp(-(6.5 + 0.95 * np.log(g) + (g / 0.15) ** 5))
    q9 = np.log(q7) * (q8 + 1 / 16.5)
    q10 = q4 - q5 / q2 * np.exp(q6 * np.log(u) / u**q9)
    ze = z0 * np.sqrt(ee0 / eee) / (1 - q4 * np.sqrt(ee0) * z0 / 377)
    zo = z0 * np.sqrt(ee0 / eeo) / (1 - q10 * np.sqrt(ee0) * z0 / 377)
    return eee, eeo, ze, zo


def _kj_dispersao(u, g, er, h_mm, eee0, eeo0):
    """f -> (εee(f), εeo(f)) de Kirschning-Jansen (eqs. 4.79, 4.81-4.82; fn = f·h ≤ 25 GHz·mm)."""
    def disp(f):
        fn = np.asarray(f, dtype=float) * 1e-9 * h_mm
        p1 = 0.27488 + (0.6315 + 0.525 / (1 + 0.0157 * fn) ** 20) * u - 0.065683 * np.exp(-8.7513 * u)
        p2 = 0.33622 * (1 - np.exp(-0.03442 * er))
        p3 = 0.0363 * np.exp(-4.6 * u) * (1 - np.exp(-((fn / 38.7) ** 4.97)))
        p4 = 1 + 2.751 * (1 - np.exp(-((er / 15.916) ** 8)))
        p5 = 0.746 + 0.334 * np.exp(-3.3 * (er / 15) ** 3)
        p6 = p5 * np.exp(-((fn / 18) ** 0.368))
        p7 = 1 + 4.069 * p6 * g**0.479 * np.exp(-1.347 * g**0.595 - 0.17 * g**2.5)
        fe = p1 * p2 * ((p3 * p4 + 0.1844 * p7) * fn) ** 1.5763
        p8 = 0.7168 * (1 + 1.076 / (1 + 0.0576 * (er - 1)))
        p9 = p8 - 0.7913 * (1 - np.exp(-((fn / 20) ** 1.424))) * np.arctan(2.481 * (er / 8) ** 0.946)
        p10 = 0.242 * (er - 1) ** 0.55
        p11 = 0.6366 * (np.exp(-0.3401 * fn) - 1) * np.arctan(1.263 * (u / 3) ** 1.629)
        p12 = p9 + (1 - p9) / (1 + 1.183 * u**1.376)
        p13 = 1.695 * p10 / (0.414 + 1.605 * p10)
        p14 = 0.8928 + 0.1072 * (1 - np.exp(-0.42 * (fn / 20) ** 3.215))
        p15 = np.abs(1 - 0.8928 * (1 + p11) * p12 * np.exp(-p13 * g**1.092) / p14)
        fo = p1 * p2 * ((p3 * p4 + 0.1844) * fn * p15) ** 1.5763
        return er - (er - eee0) / (1 + fe), er - (er - eeo0) / (1 + fo)
    return disp


def _kj_modos(w, s, h, er, t_mm=0.0):
    """(Ze, Zo, εee, εeo) com espessura pelas capacitâncias (eqs. 4.102, 4.108-4.109):
    ΔCf da linha isolada com t e Cgt = 2ε0·t/S somada ao modo ímpar (em ar e no dielétrico)."""
    u, g = w / h, s / h
    eee, eeo, ze, zo = _kj_estatico(u, g, er)
    if t_mm <= 0:
        return ze, zo, eee, eeo
    cap = lambda z, e: (np.sqrt(e) / (C0 * z), 1 / (C0 * z * np.sqrt(e)))      # noqa: E731
    ce, cae = cap(ze, eee)
    co, cao = cap(zo, eeo)
    sub0, subt = Substrato(er, h), Substrato(er, h, t_mm * 1e3, incluir_t=True)
    z_0, e_0 = (float(x) for x in analise(w, sub0))
    z_t, e_t = (float(x) for x in analise(w, subt))
    dcf = 0.5 * (cap(z_t, e_t)[0] - cap(z_0, e_0)[0])
    dcfa = 0.5 * (cap(z_t, e_t)[1] - cap(z_0, e_0)[1])
    cgt = 2 * EPS0 * t_mm / s
    ce, cae, co, cao = ce + dcf, cae + dcfa, co + dcf + cgt, cao + dcfa + cgt
    return (1 / (C0 * np.sqrt(ce * cae)), 1 / (C0 * np.sqrt(co * cao)), ce / cae, co / cao)


def microstrip_acoplada_kj(w, s, sub: Substrato, dispersivo=False):
    """Modos das microstrips acopladas por Kirschning-Jansen (erro ≤ 0,7% em εeff e 0,6% em Z para
    0,1 ≤ W/h ≤ 10, 0,1 ≤ S/h ≤ 10, 1 ≤ εr ≤ 18), com espessura se sub.incluir_t."""
    if w <= 0 or s <= 0:
        raise ValueError("W e S devem ser positivos.")
    t = sub.t_um * 1e-3 if sub.incluir_t else 0.0
    ze, zo, eee, eeo = _kj_modos(w, s, sub.h_mm, sub.er, t)
    md = Modos(float(ze), float(zo), float(eee), float(eeo))
    if dispersivo:
        md.disp = _kj_dispersao(w / sub.h_mm, s / sub.h_mm, sub.er, sub.h_mm, md.eee, md.eeo)
    return md


def alfa_condutor_modos(w, s, sub: Substrato, f, md: Modos = None):
    """(αce, αco) [Np/m] pela regra da indutância incremental de Wheeler (eqs. 4.110-4.113):
    α = Rs/(2·η0·Z)·∂Za/∂n, com todas as paredes metálicas recuando dn (W - 2dn, S + 2dn,
    h + 2dn, t - 2dn); Za = impedância do modo em ar (Kirschning-Jansen com espessura)."""
    t = max(sub.t_um, 1e-3) * 1e-3
    md = md or microstrip_acoplada_kj(w, s, sub)
    za = _kj_modos(w, s, sub.h_mm, 1.0, t)[:2]
    d = 1e-4 * min(w, s, sub.h_mm, t * 50)
    za_d = _kj_modos(w - 2 * d, s + 2 * d, sub.h_mm + 2 * d, 1.0, max(t - 2 * d, 1e-6))[:2]
    rs = resistencia_superficial(np.asarray(f, dtype=float), sub)
    dn = d * 1e-3
    return tuple(rs / (2 * ETA0 * z) * (b - a) / dn for z, a, b in zip((md.ze, md.zo), za, za_d))


def ponta_aberta_acoplada(w, s, sub: Substrato):
    """(Δle, Δlo) [mm] da extremidade aberta de linhas acopladas (eq. 4.155, ~5%)."""
    from .descontinuidades import extensao_ponta_aberta
    u, g, er, h = w / sub.h_mm, s / sub.h_mm, sub.er, sub.h_mm
    l1, l2 = extensao_ponta_aberta(w, sub), extensao_ponta_aberta(2 * w, sub)
    r1 = 1.187 * (1 - np.exp(-0.069 * u**2.1))
    r2 = 0.343 * u**0.6187 + 0.45 * er / (1 + er) * u ** (1.357 + 1.65 / (1 + 0.7 * er))
    r3 = 0.2974 * (1 - np.exp(-r2))
    r4 = ((0.271 + 0.0281 * er) * g ** (1.167 * er / (0.66 + er))
          + 1.025 * er / (0.687 + er) * g ** (0.958 * er / (0.706 + er)))
    dle = l1 + (l2 - l1 + 0.0198 * h * g**r1) * np.exp(-0.328 * g**2.244)
    dlo = l1 * (1 - np.exp(-r4)) + h * r3 * np.exp(-r4)
    return float(dle), float(dlo)


def stripline_acoplada(w, s, b, er):
    """Modos da stripline acoplada lado a lado (Cohn, t = 0); w, s, b em mm."""
    if w <= 0 or s <= 0 or b <= 0:
        raise ValueError("W, S e b devem ser positivos.")
    x1, x2 = np.pi * w / (2 * b), np.pi * (w + s) / (2 * b)
    ke = np.tanh(x1) * np.tanh(x2)
    ko = np.tanh(x1) / np.tanh(x2)
    fator = 30 * np.pi / np.sqrt(er)
    return Modos(float(fator / _razao(ke)), float(fator / _razao(ko)), er, er)


def analisar(tipo, w, s, sub: Substrato, b=0.0, modelo=MODELOS[0], dispersivo=False):
    """Modos da microstrip (Kirschning-Jansen ou Garg-Bahl) ou da stripline acoplada (Cohn)."""
    if tipo == TIPOS[0]:
        if modelo == MODELOS[1]:
            return microstrip_acoplada(w, s, sub)
        return microstrip_acoplada_kj(w, s, sub, dispersivo)
    if b <= 0:
        raise ValueError("Informe a distância entre planos b da stripline.")
    return stripline_acoplada(w, s, b, sub.er)


def sintetizar(tipo, ze, zo, sub: Substrato, b=0.0, modelo=MODELOS[0]):
    """(W, S) [mm] para os modos (ze, zo): bisseção dupla (√(ZeZo) em W e Ze/Zo em S)."""
    if not ze > zo > 0:
        raise ValueError("É preciso Z0e > Z0o > 0.")
    escala = sub.h_mm if tipo == TIPOS[0] else b
    if escala <= 0:
        raise ValueError("Informe a distância entre planos b da stripline.")
    alvo_z, alvo_r = np.sqrt(ze * zo), ze / zo
    kj = tipo == TIPOS[0] and modelo == MODELOS[0]
    s_min, s_max = (0.02, 50.0) if kj else (1e-3, 20.0)      # KJ converge para a linha isolada
    w_min, w_max = (0.05, 20.0) if kj else (1e-3, 30.0)

    def largura(s):
        lo, hi = np.log(w_min * escala), np.log(w_max * escala)
        for _ in range(60):
            m = 0.5 * (lo + hi)
            md = analisar(tipo, np.exp(m), s, sub, b, modelo)
            if md.z0 > alvo_z:
                lo = m
            else:
                hi = m
        return float(np.exp(0.5 * (lo + hi)))

    def razao(s):
        md = analisar(tipo, largura(s), s, sub, b, modelo)
        return md.ze / md.zo

    lo, hi = np.log(s_min * escala), np.log(s_max * escala)
    r_lo, r_hi = razao(np.exp(lo)), razao(np.exp(hi))
    if not r_hi < alvo_r < r_lo:
        raise ValueError(f"Z0e/Z0o = {alvo_r:.3f} fora do alcance ({r_hi:.3f} a {r_lo:.3f}) "
                         f"com separações de {s_min:g} a {s_max:g} vezes h (ou b).")
    for _ in range(50):
        m = 0.5 * (lo + hi)
        if razao(np.exp(m)) > alvo_r:
            lo = m
        else:
            hi = m
    s = float(np.exp(0.5 * (lo + hi)))
    return largura(s), s


def impedancias_acoplador(z0, c):
    """(Z0e, Z0o) para acoplamento de tensão c (eq. 7.87)."""
    return z0 * np.sqrt((1 + c) / (1 - c)), z0 * np.sqrt((1 - c) / (1 + c))


def impedancias_lange(z0, c):
    """(Z0e, Z0o) do par de linhas do Lange de 4 dedos (eq. 7.100)."""
    r = np.sqrt(9 - 8 * c**2)
    ze = (4 * c - 3 + r) / (2 * c * np.sqrt((1 - c) / (1 + c))) * z0
    zo = (4 * c + 3 - r) / (2 * c * np.sqrt((1 + c) / (1 - c))) * z0
    return float(ze), float(zo)


def modos_lange(ze, zo):
    """(Ze4, Zo4) do Lange de 4 linhas a partir do par (eq. 7.97)."""
    return (zo + ze) / (3 * zo + ze) * ze, (zo + ze) / (3 * ze + zo) * zo


# ---------------------------------------------------------------------------
# Multisseção binomial (eq. 7.90)
# ---------------------------------------------------------------------------
def coeficientes_binomial(c0, n):
    """C1..CN (simétricos) do acoplador de N seções (N ímpar) maximamente plano em θ = π/2.

    C(θ) = 2 senθ [C1 cos(N-1)θ + C2 cos(N-3)θ + ... + CM/2], M = (N+1)/2; as derivadas
    pares de ordem 2..2(M-1) se anulam em π/2 e C(π/2) = c0.
    """
    if n < 1 or n % 2 == 0:
        raise ValueError("Use número ímpar de seções.")
    m = (n + 1) // 2
    a = np.zeros((m, m))
    for k in range(m):
        if k < m - 1:
            freq = [(n - 2 * k, 1.0), (n - 2 * k - 2, -1.0)]       # 2 sen θ cos(qθ) = sen(q+1)θ - sen(q-1)θ
        else:
            freq = [(1, 1.0)]                                       # (CM/2)·2 senθ
        for lin in range(m):
            a[lin, k] = sum(sgn * (-1) ** lin * fq ** (2 * lin) * np.sin(fq * np.pi / 2) for fq, sgn in freq)
    rhs = np.zeros(m)
    rhs[0] = c0
    cs = np.linalg.solve(a, rhs)
    return list(cs) + list(cs[:-1][::-1])


# ---------------------------------------------------------------------------
# Seções e redes de 4 portas
# ---------------------------------------------------------------------------
@dataclass
class Secao:
    """Trecho de linhas acopladas de comprimento comp [mm], com modos e geometria."""
    modos: Modos
    comp: float
    w: float = 0.0
    s: float = 0.0
    # LTA: capacitores entre as linhas (só o modo ímpar os vê, como 2·Cc para o terra)
    n_cel: int = 0
    cc: float = 0.0
    tipo: str = "Microstrip"
    # extensão elétrica de cada modo pela ponta aberta (filtros de linhas acopladas) [mm]
    dl_e: float = 0.0
    dl_o: float = 0.0
    # capacitor entre as linhas em cada ponta (Steer, seção 9.10): o modo ímpar vê 2·c ao terra
    c_pontas: float = 0.0

    def gamas(self, f, sub: Substrato = None, perdas=False):
        """(γe, γo) [1/m] dos modos, com dispersão (se ativa em modos.disp) e perdas: dielétrico
        por modo e condutor por Wheeler (microstrip) ou pela linha isolada (outros casos)."""
        f = np.asarray(f, dtype=float)
        k0 = 2 * np.pi * f / C0
        eee, eeo = self.modos.eeff(f)
        ge = 1j * k0 * np.sqrt(eee)
        go = 1j * k0 * np.sqrt(eeo)
        if perdas and sub is not None:
            if self.tipo == TIPOS[0] and self.w > 0 and self.s > 0:
                ace, aco = alfa_condutor_modos(self.w, self.s, sub, f, self.modos)
            elif self.w > 0:
                ace = alfa_condutor(self.w, sub, f, z0=self.modos.ze)[:, 0]
                aco = alfa_condutor(self.w, sub, f, z0=self.modos.zo)[:, 0]
            else:
                ace = aco = 0.0
            if sub.er > 1:
                ace = ace + alfa_dieletrico(f, sub.er, eee, sub.tand)
                aco = aco + alfa_dieletrico(f, sub.er, eeo, sub.tand)
            ge, go = ge + ace, go + aco
        return ge, go

    def abcd_modos(self, f, sub=None, perdas=False):
        ge, go = self.gamas(f, sub, perdas)
        le, lo = (self.comp + self.dl_e) * 1e-3, (self.comp + self.dl_o) * 1e-3
        me = abcd.linha(self.modos.ze, -1j * ge * le)
        w = 2 * np.pi * np.asarray(f, dtype=float)
        if self.n_cel <= 0 or self.cc <= 0:
            mo = abcd.linha(self.modos.zo, -1j * go * lo)
        else:
            meia = abcd.linha(self.modos.zo, -1j * go * lo / (2 * self.n_cel))
            cel = abcd.cascata(meia, abcd.shunt(1j * w * 2 * self.cc), meia)
            mo = abcd.potencia(cel, self.n_cel)
        if self.c_pontas > 0:
            ponta = abcd.shunt(1j * w * 2 * self.c_pontas)
            mo = abcd.cascata(ponta, mo, ponta)
        return me, mo


def s4_de_modos(me, mo, z0=50.0):
    """S de 4 portas [1: a esq., 2: a dir., 3: b esq., 4: b dir.] a partir dos ABCD par/ímpar."""
    se, so = abcd.abcd_para_s(me, z0), abcd.abcd_para_s(mo, z0)
    s = np.empty(se.shape[:-2] + (4, 4), dtype=complex)
    for p in range(2):
        for q in range(2):
            soma, dif = (se[..., p, q] + so[..., p, q]) / 2, (se[..., p, q] - so[..., p, q]) / 2
            s[..., p, q] = s[..., p + 2, q + 2] = soma
            s[..., p + 2, q] = s[..., p, q + 2] = dif
    return s


def acoplador(f, secoes, z0=50.0, sub=None, perdas=False):
    """S (Nf,4,4) de seções acopladas em cascata (descontinuidades entre seções desprezadas)."""
    me = mo = None
    for sec in secoes:
        a, b = sec.abcd_modos(f, sub, perdas)
        me = a if me is None else me @ a
        mo = b if mo is None else mo @ b
    return s4_de_modos(me, mo, z0)


def terminar(s, abertas):
    """Reduz S terminando as portas `abertas` em circuito aberto (Γ = 1)."""
    n = s.shape[-1]
    q = list(abertas)
    p = [i for i in range(n) if i not in q]
    spp = s[..., p, :][..., :, p]
    spq = s[..., p, :][..., :, q]
    sqp = s[..., q, :][..., :, p]
    sqq = s[..., q, :][..., :, q]
    eye = np.eye(len(q))
    return spp + spq @ np.linalg.solve(eye - sqq, sqp)


def cascata_s(a, b):
    """Cascata de dois quadripolos pelas matrizes S (robusta onde S21 = 0)."""
    d = 1 - a[..., 1, 1] * b[..., 0, 0]
    s = np.empty(a.shape, dtype=complex)
    s[..., 0, 0] = a[..., 0, 0] + a[..., 0, 1] * b[..., 0, 0] * a[..., 1, 0] / d
    s[..., 0, 1] = a[..., 0, 1] * b[..., 0, 1] / d
    s[..., 1, 0] = b[..., 1, 0] * a[..., 1, 0] / d
    s[..., 1, 1] = b[..., 1, 1] + b[..., 1, 0] * a[..., 1, 1] * b[..., 0, 1] / d
    return s


# ---------------------------------------------------------------------------
# Projeto
# ---------------------------------------------------------------------------
def comprimento_quarto(modos: Modos, f0):
    """Comprimento [mm] de 90° na média das fases dos dois modos."""
    raiz = 0.5 * (np.sqrt(modos.eee) + np.sqrt(modos.eeo))
    return float(C0 / (4 * f0 * raiz) * 1e3)


def compensacao(modos: Modos, comp_mm, n_cel, eeo_alvo=None):
    """Cc [F] por célula que leva o εeff do modo ímpar ao alvo (padrão: εee, mesma velocidade).

    Capacitância do modo ímpar por metro: Co = √εeo/(c·Zo); carga necessária
    ΔC = Co·(alvo/εeo - 1), distribuída como 2·Cc por célula.
    """
    alvo = modos.eee if eeo_alvo is None else eeo_alvo
    if alvo <= modos.eeo:
        return 0.0
    co = np.sqrt(modos.eeo) / (C0 * modos.zo)
    dc = co * (alvo / modos.eeo - 1) * comp_mm * 1e-3
    return float(dc / (2 * n_cel))


def _em_f0(md: Modos, f0):
    """Modos com εeff avaliados em f0 (para projetar comprimentos e compensação com dispersão)."""
    eee, eeo = md.eeff(f0)
    return Modos(md.ze, md.zo, float(np.ravel(eee)[0]), float(np.ravel(eeo)[0]))


def secao_projeto(tipo, ze, zo, f0, sub: Substrato, b=0.0, lta=False, n_cel=8, modelo=MODELOS[0],
                  dispersivo=False):
    """Seção λ/4 com modos (ze, zo); com lta (microstrip), sintetiza para Zo' = Zo·√(εee/εeo) e
    calcula os Cc que baixam Zo' ao alvo e igualam as velocidades (ponto fixo).

    Retorna (Secao, aviso). Se Zo' ≥ Ze (acoplamento muito fraco), a seção fica sem compensação.
    """
    lta = lta and tipo == TIPOS[0]
    alvo_o, aviso = zo, ""
    for _ in range(8 if lta else 1):
        w, s = sintetizar(tipo, ze, alvo_o, sub, b, modelo)
        md = analisar(tipo, w, s, sub, b, modelo, dispersivo)
        m0 = _em_f0(md, f0)
        if not lta:
            break
        novo = zo * np.sqrt(m0.eee / m0.eeo)
        if novo >= 0.999 * ze:
            lta, aviso = False, "acoplamento fraco demais para compensar: seção convencional"
            w, s = sintetizar(tipo, ze, zo, sub, b, modelo)
            md = analisar(tipo, w, s, sub, b, modelo, dispersivo)
            m0 = _em_f0(md, f0)
            break
        alvo_o = novo
    if lta:
        comp = float(C0 / (4 * f0 * np.sqrt(m0.eee)) * 1e3)      # os dois modos com a fase do par
        return Secao(md, comp, w, s, n_cel, compensacao(m0, comp, n_cel), tipo), aviso
    return Secao(md, comprimento_quarto(m0, f0), w, s, tipo=tipo), aviso


def secao_capacitores_pontas(tipo, ze, zo, f0, sub: Substrato, b=0.0, modelo=MODELOS[0], dispersivo=False):
    """Seção λ/4 do modo par com capacitores entre as linhas nas duas pontas (Steer, seção 9.10).

    O modo ímpar, mais rápido, tem só θx = 90°·√(εeo/εee) no comprimento do λ/4 par; uma linha
    Zx = Zo/sen θx com derivações B = cos θx/(Zx·sen θx) nas pontas equivale em f0 à linha Zo de
    90°, e B = 2ω0·c dá o capacitor c entre as linhas. A geometria é ressintetizada para (Ze, Zx)
    até convergir. Retorna (Secao, aviso)."""
    alvo_o, aviso = zo, ""
    for _ in range(10):
        w, s = sintetizar(tipo, ze, alvo_o, sub, b, modelo)
        md = analisar(tipo, w, s, sub, b, modelo, dispersivo)
        m0 = _em_f0(md, f0)
        thx = np.pi / 2 * np.sqrt(m0.eeo / m0.eee)
        novo = zo / np.sin(thx)
        if novo >= 0.999 * ze:
            aviso = "acoplamento fraco demais para compensar: seção convencional"
            return secao_projeto(tipo, ze, zo, f0, sub, b, modelo=modelo, dispersivo=dispersivo)[0], aviso
        if abs(novo - alvo_o) < 1e-6 * zo:
            break
        alvo_o = novo
    comp = float(C0 / (4 * f0 * np.sqrt(m0.eee)) * 1e3)
    zx = m0.zo
    c = np.cos(thx) / (zx * np.sin(thx)) / (2 * 2 * np.pi * f0)
    return Secao(md, comp, w, s, tipo=tipo, c_pontas=float(max(c, 0.0))), aviso


@dataclass
class ProjetoAcoplador:
    tipo: str
    z0: float
    f0: float
    cs: list                 # acoplamentos de tensão por seção
    secoes: list
    lta: bool = False
    texto: list = field(default_factory=list)


def projeto_acoplador(tipo, z0, c_db, f0, sub: Substrato, n=1, b=0.0, lta=False, n_cel=8, lange=False,
                      modelo=MODELOS[0], dispersivo=False, pontas=False):
    """Acoplador de linhas acopladas (1 ou N seções binomiais), Lange, versão LTA ou com
    capacitores nas pontas (pontas=True, Steer 9.10; não combinado com LTA).

    LTA: a geometria é sintetizada para o modo ímpar não carregado com Zo' = Zo·√(εee/εeo)
    e os capacitores Cc baixam Zo' ao alvo e igualam as velocidades (iteração de ponto fixo).
    """
    c0 = 10 ** (-abs(c_db) / 20)
    if lange:
        cs = [c0]
    else:
        cs = coeficientes_binomial(c0, n) if n > 1 else [c0]
    secoes, txt = [], []
    for i, c in enumerate(cs):
        ze, zo = impedancias_lange(z0, c) if lange else impedancias_acoplador(z0, c)
        if pontas and not lange:
            sec, aviso = secao_capacitores_pontas(tipo, ze, zo, f0, sub, b, modelo, dispersivo)
        else:
            sec, aviso = secao_projeto(tipo, ze, zo, f0, sub, b, lta and not lange, n_cel, modelo, dispersivo)
        md = sec.modos
        if lange:
            sec.modos = Modos(*modos_lange(md.ze, md.zo), md.eee, md.eeo, md.disp)
        secoes.append(sec)
        txt.append(f"  seção {i + 1}: C = {c:.4f} ({-20 * np.log10(c):.2f} dB)  Z0e = {ze:.2f} Ω  "
                   f"Z0o = {zo:.2f} Ω" + (f"  ({aviso})" if aviso else ""))
        txt.append(f"     W = {sec.w:.4f} mm  S = {sec.s:.4f} mm  ℓ = {sec.comp:.3f} mm  εe,par = {md.eee:.3f}  "
                   f"εe,ímpar = {md.eeo:.3f}" + (f"  Cc = {sec.cc * 1e15:.1f} fF × {sec.n_cel}" if sec.n_cel else "")
                   + (f"  C pontas = {sec.c_pontas * 1e15:.1f} fF" if sec.c_pontas else ""))
    nome = "Lange" if lange else ("capacitores nas pontas" if pontas else ("LTA compensado" if lta else "convencional"))
    return ProjetoAcoplador(nome, z0, f0, cs, secoes, lta, txt)


def metricas(s):
    """(acoplamento, diretividade, direta, perda de retorno) em dB para entrada na porta 1."""
    c = -abcd.db(s[..., 2, 0])
    d = abcd.db(s[..., 2, 0]) - abcd.db(s[..., 3, 0])
    t = -abcd.db(s[..., 1, 0])
    rl = -abcd.db(s[..., 0, 0])
    return c, d, t, rl
