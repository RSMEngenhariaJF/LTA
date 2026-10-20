"""Balun de Marchand com duas seções λ/4 de linhas acopladas.

Referência: R. Garg, I. Bahl e M. Bianchi, Microstrip Lines and Slotlines, 4ª ed., 2024, seção 4.9
(o Marchand com seções λ/4 acopladas é o balun planar mais usado, de banda de múltiplas oitavas);
análise por modos par/ímpar (linhas_acopladas) montados numa rede nodal (dispositivos.Rede).

Topologia: a linha a vai da entrada (porta 1) pela seção 1 e pela seção 2 até uma ponta aberta;
as linhas b têm as pontas externas aterradas e as internas são as saídas balanceadas (portas 2 e 3).
A carga balanceada ZL é representada por duas portas de ZL/2 ao terra.
"""
from dataclasses import dataclass, field

import numpy as np

from . import linhas_acopladas as la
from .dispositivos import Rede
from .microstrip import C0, Substrato


def s_marchand(f, sec1: la.Secao, sec2: la.Secao, z_s, z_l, sub=None, perdas=False, z_ref=50.0, caps=None):
    """S de 3 portas [entrada (z_s), saída + (z_l/2), saída - (z_l/2)].

    caps: {nó: C [F]} capacitores ao terra nos nós 1 (meio da linha a), 2 (ponta aberta) e
    3/4 (saídas), usados no balun miniaturizado.
    """
    f = np.asarray(f, dtype=float)
    r = Rede(7, [0, 3, 4], f)
    # nós: 0 entrada, 1 meio da linha a, 2 ponta aberta, 3 saída da seção 1, 4 saída da seção 2,
    # 5 e 6 pontas aterradas das linhas b; portas da seção: [a esq., a dir., b esq., b dir.]
    for sec, nos in ((sec1, [0, 1, 5, 3]), (sec2, [1, 2, 4, 6])):
        me, mo = sec.abcd_modos(f, sub, perdas)
        r.multiporta(nos, la.s4_de_modos(me, mo, z_ref), z_ref)
    r.terra(5, 1e-6)
    r.terra(6, 1e-6)
    for no, c in (caps or {}).items():
        if c > 0:
            r.terra(no, 1 / (1j * 2 * np.pi * f * c))
    return r.s([z_s, z_l / 2, z_l / 2])


def _secoes_ideais(c, z_lin, f0):
    ze, zo = la.impedancias_acoplador(z_lin, c)
    comp = C0 / (4 * f0) * 1e3
    return [la.Secao(la.Modos(ze, zo, 1.0, 1.0), comp, tipo="ideal") for _ in range(2)]


def acoplamento_otimo(f0, z_s, z_l, z_lin=None):
    """Acoplamento de tensão C das duas seções (TEM ideais, Z0 = √(Z0e·Z0o) = z_lin) que casa a
    entrada em f0 (mínimo de |S11| por seção áurea)."""
    z_lin = z_lin or z_s
    f = np.array([f0])

    def s11(c):
        return abs(s_marchand(f, *_secoes_ideais(c, z_lin, f0), z_s, z_l)[0, 0, 0])
    a, b = 0.05, 0.95
    r = (np.sqrt(5) - 1) / 2
    for _ in range(80):
        c1, c2 = b - r * (b - a), a + r * (b - a)
        if s11(c1) < s11(c2):
            b = c2
        else:
            a = c1
    c = 0.5 * (a + b)
    return float(c), float(s11(c))


@dataclass
class ProjetoBalun:
    c: float
    ze: float
    zo: float
    secoes: list
    texto: list = field(default_factory=list)


def projeto(f0, z_s, z_l, sub: Substrato = None, z_lin=None, modelo=la.MODELOS[0], dispersivo=False,
            lta=False, n_cel=6):
    """Marchand: acoplamento ótimo em TEM ideal (C = 1/√(1 + ZL/Zs) com Z0 = Zs) e síntese das
    seções em microstrip (se sub), opcionalmente com compensação LTA do modo ímpar."""
    z_lin = z_lin or z_s
    c, resto = acoplamento_otimo(f0, z_s, z_l, z_lin)
    ze, zo = la.impedancias_acoplador(z_lin, c)
    p = ProjetoBalun(c, ze, zo, _secoes_ideais(c, z_lin, f0))
    p.texto = [f"Acoplamento ótimo C = {c:.4f} ({-20 * np.log10(c):.2f} dB), |S11(f0)| ideal = {resto:.2e}",
               f"Z0e = {ze:.2f} Ω, Z0o = {zo:.2f} Ω (Z0 das linhas = {z_lin:g} Ω)"]
    if sub is not None:
        try:
            sec, aviso = la.secao_projeto(la.TIPOS[0], ze, zo, f0, sub, lta=lta, n_cel=n_cel, modelo=modelo,
                                          dispersivo=dispersivo)
        except ValueError as e:
            raise ValueError(f"{e} Acoplamento de {-20 * np.log10(c):.1f} dB é forte demais para linhas lado a "
                             "lado: aumente ZL, use Z0 das linhas diferente, ou linhas acopladas pela face "
                             "(broadside) / interdigitadas (Lange).") from None
        p.secoes = [sec, la.Secao(sec.modos, sec.comp, sec.w, sec.s, sec.n_cel, sec.cc, sec.tipo)]
        p.texto.append(f"Microstrip: W = {sec.w:.4f} mm, S = {sec.s:.4f} mm, ℓ = {sec.comp:.3f} mm por seção"
                       + (f" ({aviso})" if aviso else ""))
        if sec.n_cel:
            p.texto.append(f"LTA: {sec.n_cel} capacitores Cc = {sec.cc * 1e15:.2f} fF por seção")
    return p


def metricas(s):
    """(|S11| dB, |S21| dB, |S31| dB, desequilíbrio de amplitude dB, fase S21 - S31 em graus)."""
    db = lambda x: 20 * np.log10(np.maximum(np.abs(x), 1e-15))      # noqa: E731
    fase = np.degrees(np.angle(s[:, 1, 0] * np.conj(s[:, 2, 0])))
    return db(s[:, 0, 0]), db(s[:, 1, 0]), db(s[:, 2, 0]), db(s[:, 1, 0]) - db(s[:, 2, 0]), fase


@dataclass
class BalunMiniaturizado:
    c: float                 # acoplamento de tensão das seções
    theta0: float            # comprimento elétrico das seções em f0 [graus]
    caps: dict               # {nó: C [F]}
    secoes: list
    erro: float
    texto: list = field(default_factory=list)
    caps_ideais: dict = field(default_factory=dict)      # projeto com linhas TEM ideais (referência)
    secoes_ideais: list = field(default_factory=list)


def _secoes_theta(c, z_lin, f0, theta0, modos=None):
    ze, zo = la.impedancias_acoplador(z_lin, c)
    md = modos or la.Modos(ze, zo, 1.0, 1.0)
    raiz = 0.5 * (np.sqrt(md.eee) + np.sqrt(md.eeo))
    comp = theta0 / 360 * C0 / (f0 * raiz) * 1e3
    return [la.Secao(md, comp, tipo="ideal") for _ in range(2)]


def projeto_miniaturizado(f0, z_s, z_l, theta0=45.0, z_lin=None, sub=None, modelo=la.MODELOS[0], c_fixo=None):
    """Marchand com seções de θ0 < 90° em f0 carregadas por capacitores (Steer, seção 13.2: linhas
    λ/4 na frequência comensurada fr > f0 e stubs abertos aproximados por capacitores).

    Os capacitores do meio da linha a (nó 1), da ponta aberta (nó 2) e das saídas (nós 3 e 4) são
    ajustados por Nelder-Mead para |S11| = 0, |S21| = |S31| e fase de 180° em f0: primeiro com o
    acoplamento do Marchand convencional (com um termo pequeno que favorece capacitâncias
    menores) e depois refinando também o acoplamento, que fica um pouco mais forte (a menos que
    c_fixo seja dado); com linhas TEM ideais e depois, se sub, com as microstrips acopladas
    sintetizadas (mesmo θ0 na média dos modos, só os capacitores reajustados).
    """
    from .otimizacao import nelder_mead
    if not 10 <= theta0 <= 90:
        raise ValueError("Use θ0 entre 10° e 90°.")
    z_lin = z_lin or z_s
    w0 = 2 * np.pi * f0
    c_ref = 1 / (w0 * z_lin)
    f = np.array([f0])
    c = c_fixo if c_fixo else acoplamento_otimo(f0, z_s, z_l, z_lin)[0]
    if not 0 < c < 1:
        raise ValueError("O acoplamento deve estar entre 0 e 1.")
    ze, zo = la.impedancias_acoplador(z_lin, c)
    secs = _secoes_theta(c, z_lin, f0, theta0)

    def caps_de(y):
        e = np.exp(np.clip(y, -30, 5))
        return {1: c_ref * e[0], 2: c_ref * e[1], 3: c_ref * e[2], 4: c_ref * e[3]}

    def erro(y, secoes, peso=1e-4):
        s = s_marchand(f, *secoes, z_s, z_l, caps=caps_de(y))[0]
        bal = abs(s[1, 0]) - abs(s[2, 0])
        fase = np.angle(-s[1, 0] * np.conj(s[2, 0]))                 # 0 quando a diferença é 180°
        reg = peso * float(np.sum(np.exp(np.clip(y, -30, 5))))
        return abs(s[0, 0]) ** 2 + bal**2 + 0.1 * fase**2 + reg

    melhor = None
    for y0 in ([-2.0, -2.0, -2.0, -2.0], [-4.0, 0.0, -1.0, -1.0], [0.0, -4.0, 0.0, 0.0]):
        y, v = nelder_mead(lambda q: erro(q, secs), y0, passo=0.7, max_iter=600)
        if melhor is None or v < melhor[1]:
            melhor = (y, v)
    y, v = melhor
    if c_fixo:
        y, v = nelder_mead(lambda q: erro(q, secs, 1e-8), y, passo=0.1, max_iter=600)   # refino quase sem penalidade
    else:
        def erro_c(q):
            return erro(q[:4], _secoes_theta(1 / (1 + np.exp(-q[4])), z_lin, f0, theta0), 1e-8)
        q, v = nelder_mead(erro_c, list(y) + [np.log(c / (1 - c))], passo=0.2, tol=1e-14, max_iter=1200)
        y, c = q[:4], float(1 / (1 + np.exp(-q[4])))
        ze, zo = la.impedancias_acoplador(z_lin, c)
        secs = _secoes_theta(c, z_lin, f0, theta0)
    caps = caps_de(y)
    s0 = s_marchand(f, *secs, z_s, z_l, caps=caps)[0]
    texto = [f"Seções de {theta0:g}° em f0 (λ/4 em fr = {90 / theta0 * f0 / 1e9:.3f} GHz), acoplamento "
             f"C = {c:.4f} ({-20 * np.log10(c):.2f} dB): Z0e = {ze:.2f} Ω, Z0o = {zo:.2f} Ω",
             "Capacitores (TEM ideal): meio da linha a {:.4f} pF, ponta aberta {:.4f} pF, saídas {:.4f} / {:.4f} pF"
             .format(caps[1] * 1e12, caps[2] * 1e12, caps[3] * 1e12, caps[4] * 1e12),
             f"|S11(f0)| = {20 * np.log10(max(abs(s0[0, 0]), 1e-12)):.1f} dB"]
    caps_id, secs_id = dict(caps), list(secs)
    if sub is not None:
        try:
            sec0, aviso = la.secao_projeto(la.TIPOS[0], ze, zo, f0, sub, modelo=modelo)
        except ValueError as e:
            raise ValueError(f"{e} Acoplamento de {-20 * np.log10(c):.1f} dB forte demais para microstrip lado a "
                             "lado: aumente ZL, mude Z0 das linhas ou fixe um acoplamento menor.") from None
        md = sec0.modos
        comp = theta0 / 90 * la.comprimento_quarto(md, f0)
        secs = [la.Secao(md, comp, sec0.w, sec0.s, tipo=sec0.tipo) for _ in range(2)]
        y, v = nelder_mead(lambda q: erro(q, secs, 1e-8), y, passo=0.3, max_iter=800)
        caps = caps_de(y)
        s0 = s_marchand(f, *secs, z_s, z_l, caps=caps)[0]
        texto += [f"Microstrip: W = {sec0.w:.4f} mm, S = {sec0.s:.4f} mm, ℓ = {comp:.3f} mm por seção "
                  f"({100 * theta0 / 90:.0f}% do λ/4)" + (f" ({aviso})" if aviso else ""),
                  "Capacitores reajustados: meio {:.4f} pF, ponta {:.4f} pF, saídas {:.4f} / {:.4f} pF; |S11(f0)| = {:.1f} dB"
                  .format(caps[1] * 1e12, caps[2] * 1e12, caps[3] * 1e12, caps[4] * 1e12,
                          20 * np.log10(max(abs(s0[0, 0]), 1e-12)))]
        if sec0.s < 0.1:
            texto.append(f"Atenção: S = {sec0.s * 1e3:.0f} µm está abaixo do mínimo típico de PCB (~100 µm); "
                         "considere linhas acopladas pela face (broadside) ou Lange.")
    return BalunMiniaturizado(float(c), theta0, caps, secs, float(v), texto, caps_id, secs_id)
