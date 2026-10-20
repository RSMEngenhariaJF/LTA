"""Filtros com linhas de transmissão convencionais (Pozar, cap. 8) e variante com LTA.

Referência: D. M. Pozar, Microwave Engineering, 4ª ed., Wiley, 2012:
  - passa-baixas por degraus de impedância (seção 8.6, eq. 8.86, Ex. 8.6);
  - passa-baixas com stubs: transformação de Richards e identidades de Kuroda (8.4-8.5, Ex. 8.5);
  - passa-faixa com linhas acopladas em paralelo (8.7, eqs. 8.108 e 8.121, Ex. 8.7);
  - rejeita-faixa com stubs λ/4 abertos e passa-faixa com stubs λ/4 em curto (8.8, eqs. 8.130-8.131);
  - passa-faixa com ressonadores série acoplados por gaps (8.8, eqs. 8.132-8.135, Ex. 8.9) e com
    ressonadores em derivação acoplados por capacitores (eqs. 8.136-8.141, Ex. 8.10).
Protótipos g: filtros.prototipo (Butterworth e Chebyshev, Pozar Tabelas 8.3-8.4).

LTA: no passa-faixa de linhas acopladas em microstrip, a diferença de velocidade dos modos
cria a banda espúria em 2f0 (Martín 2015, seção 4.2.3.2). A opção compensada carrega o modo
ímpar com capacitores periódicos entre as linhas (linhas_acopladas.compensacao).
"""
from dataclasses import dataclass, field

import numpy as np

from . import abcd
from . import linhas_acopladas as la
from .microstrip import C0, Substrato


# ---------------------------------------------------------------------------
# Elementos ideais comensuráveis (θ = θ0·f/f0)
# ---------------------------------------------------------------------------
def _theta(f, f0, theta0):
    return theta0 * np.asarray(f, dtype=float) / f0


def abcd_ue(f, f0, z, theta0):
    return abcd.linha(z, _theta(f, f0, theta0))


def abcd_stub_derivacao(f, f0, z, theta0, curto=False):
    t = np.tan(_theta(f, f0, theta0))
    return abcd.shunt(-1j / (z * t) if curto else 1j * t / z)


def abcd_stub_serie(f, f0, z, theta0, curto=True):
    t = np.tan(_theta(f, f0, theta0))
    return abcd.series(1j * z * t if curto else -1j * z / t)


# ---------------------------------------------------------------------------
# Passa-baixas por degraus de impedância
# ---------------------------------------------------------------------------
@dataclass
class Trecho:
    z: float             # Ω
    bl: float            # comprimento elétrico em fc [rad]
    elemento: str        # "L" ou "C"


def degraus(g, r0, zl, zh, comeca_com_c=True):
    """Trechos de linha (eq. 8.86): βℓ = g·R0/Zh (indutor) e βℓ = g·Zℓ/R0 (capacitor)."""
    if not zh > r0 > zl > 0:
        raise ValueError("É preciso Zh > R0 > Zℓ > 0.")
    out = []
    for i, gi in enumerate(g[1:-1]):
        capacitor = (i % 2 == 0) == comeca_com_c
        out.append(Trecho(zl, gi * zl / r0, "C") if capacitor else Trecho(zh, gi * r0 / zh, "L"))
    return out


def correcao_degraus(segs, sub, w_entrada=None, w_saida=None):
    """Comprimentos físicos [(comp, W)] corrigidos pelos degraus de largura (Garg-Bahl-Bianchi,
    seção 3.4.3): em cada degrau a linha larga fica Δl mais longa e a estreita Δl mais curta, então
    a larga é encurtada e a estreita alongada de Δl. Com w_entrada/w_saida (alimentações, que não
    mudam), os degraus das pontas corrigem só o trecho adjacente. Retorna (trechos, lista de Δl)."""
    from .descontinuidades import degrau
    comps = [c for c, _ in segs]
    dls = []
    for i in range(len(segs) - 1):
        (_, w1), (_, w2) = segs[i], segs[i + 1]
        dl = degrau(w1, w2, sub)
        dls.append(dl)
        larga, estreita = (i, i + 1) if w1 > w2 else (i + 1, i)
        comps[larga] -= dl
        comps[estreita] += dl
    for w_ext, i in ((w_entrada, 0), (w_saida, len(segs) - 1)):
        if w_ext and abs(w_ext - segs[i][1]) > 1e-9:
            dl = degrau(w_ext, segs[i][1], sub)
            dls.append(dl)
            comps[i] += -dl if segs[i][1] > w_ext else dl
    if min(comps) <= 0:
        raise ValueError("Trecho curto demais para a correção de degrau.")
    return [(c, w) for c, (_, w) in zip(comps, segs)], dls


def abcd_degraus_microstrip(f, segs, sub, perdas=False, com_degraus=True, w_entrada=None, w_saida=None):
    """Cascata dos trechos em microstrip real; com_degraus inclui o efeito de cada degrau como
    deslocamento de plano de referência (+Δl na larga, -Δl na estreita); w_entrada/w_saida são as
    larguras das alimentações (o deslocamento do lado delas entra como trecho de linha negativo
    ou positivo da própria alimentação)."""
    from .descontinuidades import degrau
    from .microstrip import parametros_linha
    f = np.asarray(f, dtype=float)
    segs = ([(0.0, w_entrada)] if w_entrada else []) + list(segs) + ([(0.0, w_saida)] if w_saida else [])
    comps = [c for c, _ in segs]
    if com_degraus:
        for i in range(len(segs) - 1):
            (_, w1), (_, w2) = segs[i], segs[i + 1]
            if abs(w1 - w2) < 1e-12:
                continue
            dl = degrau(w1, w2, sub)
            larga, estreita = (i, i + 1) if w1 > w2 else (i + 1, i)
            comps[larga] += dl
            comps[estreita] -= dl
    m = None
    for comp, (_, w) in zip(comps, segs):
        z, ee, a = (q[:, 0] for q in parametros_linha(w, sub, f, perdas, False))
        mi = abcd.linha(z, 2 * np.pi * f * np.sqrt(ee) / C0 * comp * 1e-3 - 1j * a * comp * 1e-3)
        m = mi if m is None else m @ mi
    return m


def abcd_degraus_ideal(f, fc, trechos):
    return abcd.cascata(*[abcd.linha(t.z, t.bl * np.asarray(f, dtype=float) / fc) for t in trechos])


def abcd_concentrado(f, fc, g, r0, comeca_com_c=True):
    """Protótipo concentrado escalado (comparação, Fig. 8.41)."""
    w = 2 * np.pi * np.asarray(f, dtype=float)
    wc = 2 * np.pi * fc
    mats = []
    for i, gi in enumerate(g[1:-1]):
        if (i % 2 == 0) == comeca_com_c:
            mats.append(abcd.shunt(1j * w * gi / (r0 * wc)))
        else:
            mats.append(abcd.series(1j * w * gi * r0 / wc))
    return abcd.cascata(*mats)


# ---------------------------------------------------------------------------
# Passa-baixas com stubs (Richards + Kuroda)
# ---------------------------------------------------------------------------
@dataclass
class ElementoK:
    tipo: str     # "UE" (linha unitária), "P" (stub aberto em derivação), "S" (stub em curto em série)
    z: float      # impedância normalizada


def richards_kuroda(g, comeca_com_c=False):
    """Elementos λ/8 em fc após Richards (L → stub série Z = L, C → stub em derivação Z = 1/C) e
    Kuroda (UE(Z1) + stub série(Z2) ≡ stub em derivação(n²Z1) + UE(n²Z2), n² = 1 + Z1/Z2).

    Só stubs série adjacentes às pontas podem ser convertidos com as linhas redundantes: vale
    para N ≤ 3 começando com L (L, L-C, L-C-L) e N ≤ 2 começando com C (C, C-L).
    """
    els = []
    for i, gi in enumerate(g[1:-1]):
        capacitor = (i % 2 == 0) == comeca_com_c
        els.append(ElementoK("P", 1 / gi) if capacitor else ElementoK("S", gi))
    if els[0].tipo == "S":
        z1, z2 = 1.0, els[0].z
        n2 = 1 + z1 / z2
        els[0:1] = [ElementoK("P", n2 * z1), ElementoK("UE", n2 * z2)]
    if els[-1].tipo == "S":
        z1, z2 = g[-1], els[-1].z                     # linha redundante casada à carga (g(N+1))
        n2 = 1 + z1 / z2
        els[-1:] = [ElementoK("UE", n2 * z2), ElementoK("P", n2 * z1)]
    if any(e.tipo == "S" for e in els):
        raise ValueError("Kuroda com linhas redundantes só converte stubs série nas pontas: use N ≤ 3 "
                         "começando com L (ou N ≤ 2 com C), ou o filtro de degraus de impedância.")
    return els


def abcd_elementos(f, fc, els, z0=1.0, theta0=np.pi / 4):
    mats = []
    for e in els:
        z = e.z * z0
        if e.tipo == "UE":
            mats.append(abcd_ue(f, fc, z, theta0))
        elif e.tipo == "P":
            mats.append(abcd_stub_derivacao(f, fc, z, theta0))
        else:
            mats.append(abcd_stub_serie(f, fc, z, theta0))
    return abcd.cascata(*mats)


def abcd_richards(f, fc, g, z0, comeca_com_c=False):
    """Rede após Richards, antes de Kuroda (stubs série e em derivação, λ/8 em fc)."""
    els = []
    for i, gi in enumerate(g[1:-1]):
        capacitor = (i % 2 == 0) == comeca_com_c
        els.append(ElementoK("P", 1 / gi) if capacitor else ElementoK("S", gi))
    return abcd_elementos(f, fc, els, z0)


# ---------------------------------------------------------------------------
# Passa-faixa com linhas acopladas
# ---------------------------------------------------------------------------
def inversores_acoplados(g, fbw):
    """Z0·Jn, n = 1..N+1 (eq. 8.121)."""
    n = len(g) - 2
    j = [np.sqrt(np.pi * fbw / (2 * g[1]))]
    for k in range(2, n + 1):
        j.append(np.pi * fbw / (2 * np.sqrt(g[k - 1] * g[k])))
    j.append(np.sqrt(np.pi * fbw / (2 * g[n] * g[n + 1])))
    return np.array(j)


def modos_inversor(z0, jz):
    """(Z0e, Z0o) da seção de linhas acopladas com inversor J (eq. 8.108)."""
    return z0 * (1 + jz + jz**2), z0 * (1 - jz + jz**2)


@dataclass
class ProjetoAcopladas:
    jz: np.ndarray
    ze: list
    zo: list
    secoes: list = field(default_factory=list)      # la.Secao (geometria física), se sintetizada
    texto: list = field(default_factory=list)


def projeto_linhas_acopladas(g, f0, fbw, z0, sub: Substrato = None, tipo=la.TIPOS[0], b=0.0,
                             lta=False, n_cel=6, corrigir_ponta=False, modelo=la.MODELOS[0], dispersivo=False):
    """Passa-faixa de N+1 seções λ/4 acopladas; com sub, sintetiza W e S de cada seção.

    Microstrip: cada modo fica eletricamente mais longo pela ponta aberta acoplada (Δle, Δlo,
    Garg-Bahl-Bianchi eq. 4.155), sempre incluída no modelo; corrigir_ponta desconta a média
    dos dois do comprimento físico.
    """
    jz = inversores_acoplados(g, fbw)
    zs = [modos_inversor(z0, j) for j in jz]
    p = ProjetoAcopladas(jz, [z[0] for z in zs], [z[1] for z in zs])
    if sub is None:
        return p
    for i, (ze, zo) in enumerate(zs):
        sec, aviso = la.secao_projeto(tipo, ze, zo, f0, sub, b, lta, n_cel, modelo, dispersivo)
        if tipo == la.TIPOS[0]:
            sec.dl_e, sec.dl_o = la.ponta_aberta_acoplada(sec.w, sec.s, sub)
            if corrigir_ponta:
                sec.comp -= 0.5 * (sec.dl_e + sec.dl_o)
            if sec.n_cel:
                # LTA: iguala o comprimento elétrico total dos modos, ponta aberta incluída:
                # (ℓ + Δlo)·√εeo' = (ℓ + Δle)·√εee
                m0 = la._em_f0(sec.modos, f0)
                alvo = m0.eee * ((sec.comp + sec.dl_e) / (sec.comp + sec.dl_o)) ** 2
                sec.cc = la.compensacao(m0, sec.comp, sec.n_cel, alvo)
        p.secoes.append(sec)
        p.texto.append(f"  {i + 1}: Z0·J = {jz[i]:.4f}  Z0e = {ze:.2f}  Z0o = {zo:.2f} Ω  W = {sec.w:.4f}  "
                       f"S = {sec.s:.4f}  ℓ = {sec.comp:.3f} mm"
                       + (f"  Cc = {sec.cc * 1e15:.1f} fF × {sec.n_cel}" if sec.n_cel else "")
                       + (f"  ({aviso})" if aviso else ""))
    return p


def s_linhas_acopladas(f, secoes, z0=50.0, sub=None, perdas=False):
    """S (2 portas) do filtro: cada seção usa a entrada na linha a (esq.) e a saída na b (dir.),
    com as outras duas pontas abertas (Pozar, Fig. 8.44)."""
    s = None
    for sec in secoes:
        me, mo = sec.abcd_modos(f, sub, perdas)
        s2 = la.terminar(la.s4_de_modos(me, mo, z0), [1, 2])
        s = s2 if s is None else la.cascata_s(s, s2)
    return s


def secoes_ideais(f0, ze, zo):
    """Seções TEM ideais (mesma velocidade nos modos, εeff = 1) de 90° em f0."""
    comp = C0 / (4 * f0) * 1e3
    return [la.Secao(la.Modos(a, b, 1.0, 1.0), comp) for a, b in zip(ze, zo)]


# ---------------------------------------------------------------------------
# Stubs λ/4 (rejeita-faixa aberto, passa-faixa em curto)
# ---------------------------------------------------------------------------
def stubs_quarto(g, fbw, z0, rejeita=True):
    """Z0n dos stubs (eq. 8.130: 4Z0/(π·gn·Δ); eq. 8.131: π·Z0·Δ/(4·gn))."""
    gn = np.asarray(g[1:-1], dtype=float)
    if not np.isclose(g[-1], 1.0):
        raise ValueError("Estas fórmulas exigem g(N+1) = 1 (Chebyshev de ordem ímpar ou Butterworth).")
    return list(4 * z0 / (np.pi * gn * fbw)) if rejeita else list(np.pi * z0 * fbw / (4 * gn))


def abcd_stubs_quarto(f, f0, zs, z0, rejeita=True):
    mats = []
    for i, z in enumerate(zs):
        if i:
            mats.append(abcd_ue(f, f0, z0, np.pi / 2))
        mats.append(abcd_stub_derivacao(f, f0, z, np.pi / 2, curto=not rejeita))
    return abcd.cascata(*mats)


# ---------------------------------------------------------------------------
# Ressonadores acoplados por capacitores
# ---------------------------------------------------------------------------
@dataclass
class ProjetoGaps:
    jz: np.ndarray
    b: np.ndarray        # susceptâncias [S]
    c: np.ndarray        # capacitores série [F]
    theta: np.ndarray    # comprimentos elétricos dos ressonadores [rad]


def ressonadores_serie_gap(g, f0, fbw, z0):
    """Ressonadores λ/2 série entre gaps (eqs. 8.121, 8.134, 8.135)."""
    jz = inversores_acoplados(g, fbw)
    j = jz / z0
    b = j / (1 - jz**2)
    c = b / (2 * np.pi * f0)
    fi = np.arctan(2 * z0 * b)
    theta = np.pi - 0.5 * (fi[:-1] + fi[1:])
    return ProjetoGaps(jz, b, c, theta)


def abcd_gaps_ideal(f, f0, p: ProjetoGaps, z0, cps=None):
    """Capacitores série + linhas (θ = θi·f/f0); cps: capacitâncias em derivação dos gaps (Π)."""
    f = np.asarray(f, dtype=float)
    w = 2 * np.pi * f
    mats = []
    for i, c in enumerate(p.c):
        if cps is not None:
            mats.append(abcd.shunt(1j * w * cps[i]))
        mats.append(abcd.series(1 / (1j * w * c)))
        if cps is not None:
            mats.append(abcd.shunt(1j * w * cps[i]))
        if i < len(p.theta):
            mats.append(abcd.linha(z0, p.theta[i] * f / f0))
    return abcd.cascata(*mats)


@dataclass
class ProjetoDerivacao:
    jz: np.ndarray
    c: np.ndarray        # capacitores de acoplamento C01..CN,N+1 [F]
    dc: np.ndarray       # ΔCn [F] (eq. 8.138)
    l: np.ndarray        # comprimentos dos stubs em curto [λ] (eq. 8.141)


def ressonadores_derivacao(g, f0, fbw, z0):
    """Stubs em curto (≈λ/4, impedância Z0) acoplados por capacitores série (eqs. 8.136-8.141)."""
    n = len(g) - 2
    w0 = 2 * np.pi * f0
    jz = [np.sqrt(np.pi * fbw / (4 * g[1]))]
    for k in range(1, n):
        jz.append(np.pi * fbw / (4 * np.sqrt(g[k] * g[k + 1])))
    jz.append(np.sqrt(np.pi * fbw / (4 * g[n] * g[n + 1])))
    jz = np.array(jz)
    j = jz / z0
    c = j / w0
    c[0] = j[0] / (w0 * np.sqrt(1 - jz[0] ** 2))
    c[-1] = j[-1] / (w0 * np.sqrt(1 - jz[-1] ** 2))
    dc = -(c[:-1] + c[1:])
    l = 0.25 + z0 * w0 * dc / (2 * np.pi)
    return ProjetoDerivacao(jz, c, dc, l)


def abcd_derivacao_ideal(f, f0, p: ProjetoDerivacao, z0):
    f = np.asarray(f, dtype=float)
    w = 2 * np.pi * f
    mats = [abcd.series(1 / (1j * w * p.c[0]))]
    for i, l in enumerate(p.l):
        mats.append(abcd_stub_derivacao(f, f0, z0, 2 * np.pi * l, curto=True))
        mats.append(abcd.series(1 / (1j * w * p.c[i + 1])))
    return abcd.cascata(*mats)
