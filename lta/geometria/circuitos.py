"""Layouts de circuitos com linhas (Pozar, caps. 5-8): casamento, filtros, acopladores e divisores.

Mesmas convenções de layout.py (mm, porta 1 em x = 0, portas nas faces x = xmin/xmax).
Os comprimentos são medidos nas linhas de centro; as curvas podem ser chanfradas (Garg, Bahl e
Bianchi, eq. 3.115); junções T e degraus não são compensados no desenho e devem ser ajustados em
simulação EM.
"""
import numpy as np

from ..core.microstrip import Substrato
from .layout import Layout, Porta, poligono_de_perfil, simplificar_perfil


def _ret(x0, y0, x1, y1):
    xa, xb = sorted((x0, x1))
    ya, yb = sorted((y0, y1))
    return np.array([[xa, ya], [xb, ya], [xb, yb], [xa, yb]], dtype=float)


def _cota_h(x1, x2, y, texto):
    return (x1, y, x2, y, texto)


def _cota_v(x, y1, y2, texto):
    return (x, y1, x, y2, texto)


def caminho(pontos, w, chanfro=None):
    """Trilha de largura w ao longo de segmentos horizontais/verticais.

    chanfro: None (cantos retos) ou M [%] do chanfro das curvas de 90° (Garg-Bahl-Bianchi,
    eq. 3.115: M = X/d·100, d = diagonal do quadrado do canto, X = trecho cortado a partir do
    canto externo). O corte é um triângulo de catetos 2·W·M/100 no canto externo.
    """
    pts = [np.asarray(p, dtype=float) for p in pontos]
    if chanfro and len(pts) > 2:
        return _caminho_chanfrado(pts, w, chanfro)
    polis = []
    for i, (a, b) in enumerate(zip(pts[:-1], pts[1:])):
        e0 = w / 2 if i > 0 else 0.0
        e1 = w / 2 if i < len(pts) - 2 else 0.0
        if abs(a[1] - b[1]) < 1e-12:                          # horizontal
            s = np.sign(b[0] - a[0]) or 1.0
            polis.append(_ret(a[0] - s * e0, a[1] - w / 2, b[0] + s * e1, a[1] + w / 2))
        elif abs(a[0] - b[0]) < 1e-12:                        # vertical
            s = np.sign(b[1] - a[1]) or 1.0
            polis.append(_ret(a[0] - w / 2, a[1] - s * e0, a[0] + w / 2, b[1] + s * e1))
        else:
            raise ValueError("Caminho só com segmentos horizontais e verticais.")
    return polis


def _caminho_chanfrado(pts, w, m):
    leg = 2 * w * m / 100
    extra = max(leg - w, 0.0)
    n = len(pts) - 1
    polis = []
    for i, (a, b) in enumerate(zip(pts[:-1], pts[1:])):
        d = (b - a) / np.linalg.norm(b - a)
        ini = a + d * (w / 2 + extra) if i > 0 else a
        fim = b - d * (w / 2 + extra) if i < n - 1 else b
        if np.dot(fim - ini, d) <= 1e-9:
            continue
        polis.append(_ret(ini[0] - abs(d[1]) * w / 2, ini[1] - abs(d[0]) * w / 2,
                          fim[0] + abs(d[1]) * w / 2, fim[1] + abs(d[0]) * w / 2))
    for i in range(1, n):
        p = pts[i]
        ain = (p - pts[i - 1]) / np.linalg.norm(p - pts[i - 1])
        out = (pts[i + 1] - p) / np.linalg.norm(pts[i + 1] - p)
        loc = [(-w / 2 - extra, -w / 2), (w / 2 - leg, -w / 2), (w / 2, -w / 2 + leg), (w / 2, w / 2 + extra),
               (-w / 2, w / 2 + extra), (-w / 2, w / 2), (-w / 2 - extra, w / 2)]
        q = np.array([p + x * ain + y * out for x, y in loc])
        manter = np.r_[True, np.any(np.abs(np.diff(q, axis=0)) > 1e-12, axis=1)]
        polis.append(q[manter])
    return polis


def chanfro_para(w, sub, f0=None):
    """(M [%], Δb [mm]) da curva chanfrada ótima: M pela eq. 3.115 e Δb pela eq. 3.116 (se f0)."""
    from ..core.descontinuidades import chanfro_otimo, encurtamento_chanfro
    return chanfro_otimo(w, sub), (encurtamento_chanfro(w, sub, f0) if f0 else 0.0)


def trechos_eixo(origem, eixo, trechos):
    """Retângulos de [(comprimento, largura), ...] ao longo de +x ('x') ou -y ('-y') a partir de origem."""
    x, y = origem
    polis = []
    for comp, larg in trechos:
        if eixo == "x":
            polis.append(_ret(x, y - larg / 2, x + comp, y + larg / 2))
            x += comp
        else:
            polis.append(_ret(x - larg / 2, y, x + larg / 2, y - comp))
            y -= comp
    return polis


def braco_slow_wave(origem, eixo, geo, n, lado=1):
    """Polígonos de n células slow-wave (dispositivos.GeometriaSW) ao longo de +x ou -y.

    Patches: trechos largos centrados em cada célula. Stubs: linha uniforme e um stub aberto
    no centro de cada célula, para o lado `lado` (+1: +y no eixo x, +x no eixo -y).
    """
    if geo.tipo == "Patches":
        resto = (geo.periodo - geo.l_carga) / 2
        tr = []
        for _ in range(int(n)):
            tr += [(resto, geo.w_linha), (geo.l_carga, geo.w_carga), (resto, geo.w_linha)]
        return trechos_eixo(origem, eixo, tr)
    comp = n * geo.periodo
    polis = trechos_eixo(origem, eixo, [(comp, geo.w_linha)])
    x, y = origem
    for k in range(int(n)):
        c = (k + 0.5) * geo.periodo
        if eixo == "x":
            y0 = y + lado * geo.w_linha / 2
            polis.append(_ret(x + c - geo.w_carga / 2, y0, x + c + geo.w_carga / 2, y0 + lado * geo.l_carga))
        else:
            x0 = x + lado * geo.w_linha / 2
            polis.append(_ret(x0, y - c - geo.w_carga / 2, x0 + lado * geo.l_carga, y - c + geo.w_carga / 2))
    return polis


def _margem(sub):
    return max(5 * sub.h_mm, 2.0)


def _placa(topo, sub, x_fim, extra_y=0.0):
    ys = np.concatenate([p[:, 1] for p in topo])
    m = _margem(sub) + extra_y
    return (0.0, float(ys.min() - m), float(x_fim), float(ys.max() + m))


# ---------------------------------------------------------------------------
# Linha principal com stubs (casamento por stubs, filtros com stubs)
# ---------------------------------------------------------------------------
def linha_com_stubs(nome, sub: Substrato, segmentos, stubs, w_in, w_out, l_alim=5.0, l_saida=None,
                    raio_via=0.2, info=None):
    """Linha ao longo de y = 0 com stubs para cima.

    segmentos: [(comprimento, largura)] da linha principal, a partir de x = l_alim.
    stubs: [dict(x=posição ao longo da linha principal, w, comp, curto)], comprimento medido a
    partir do eixo da linha principal. l_saida: alimentação após a linha (None = l_alim; 0 =
    porta 2 no fim da linha principal, como no plano da carga).
    """
    l_saida = l_alim if l_saida is None else l_saida
    trechos = ([(l_alim, w_in)] if l_alim > 0 else []) + list(segmentos) + ([(l_saida, w_out)] if l_saida > 0
                                                                              else [])
    topo = trechos_eixo((0.0, 0.0), "x", trechos)
    x_fim = sum(c for c, _ in trechos)
    vias, cotas = [], []
    y_max = 0.0
    for k, st in enumerate(stubs):
        xc = l_alim + st["x"]
        topo.append(_ret(xc - st["w"] / 2, 0.0, xc + st["w"] / 2, st["comp"]))
        topo_y = st["comp"]
        if st.get("curto"):
            lado = max(2 * raio_via + 0.4, st["w"])
            topo.append(_ret(xc - lado / 2, topo_y, xc + lado / 2, topo_y + lado))
            vias.append((xc, topo_y + lado / 2, raio_via))
            topo_y += lado
        y_max = max(y_max, topo_y)
        if k < 4:
            cotas.append(_cota_v(xc + st["w"] / 2 + 0.3, 0.0, st["comp"], f"stub {k + 1} = {st['comp']:.3f}"))
    for k, (comp, larg) in enumerate(segmentos[:3]):
        x0 = l_alim + sum(c for c, _ in segmentos[:k])
        cotas.append(_cota_h(x0, x0 + comp, -max(w for _, w in trechos) * 0.8 - 0.3 * k,
                             f"ℓ{k + 1} = {comp:.3f} (W {larg:.3f})"))
    w_max = max(w for _, w in trechos)
    m = _margem(sub)
    placa = (0.0, -w_max / 2 - m, float(x_fim), max(y_max, w_max / 2) + m)
    portas = [Porta(0.0, 0.0, w_in, "xmin"), Porta(float(x_fim), 0.0, w_out if l_saida > 0 else
                                                   segmentos[-1][1], "xmax")]
    xs, ws = [], []
    x = 0.0
    for comp, larg in trechos:
        xs += [x, x + comp]
        ws += [larg, larg]
        x += comp
    lay = Layout(nome, sub, topo, placa, portas=portas, info=info or {}, vias=vias,
                 perfil=(np.array(xs), np.array(ws)))
    lay.celula = (max(0.0, l_alim - 1), min(x_fim, l_alim + sum(c for c, _ in segmentos) + 1))
    lay.cotas = cotas
    return lay


def linha_perfil(nome, sub: Substrato, x, w, w_in, w_out, l_alim=5.0, info=None, celula=None):
    """Trilha com perfil W(x) (taper ou degraus) e alimentações de larguras w_in e w_out."""
    x = np.asarray(x, dtype=float) + l_alim
    xs = np.concatenate([[0.0, l_alim], x, [x[-1], x[-1] + l_alim]])
    ws = np.concatenate([[w_in, w_in], np.asarray(w, dtype=float), [w_out, w_out]])
    xs, ws = simplificar_perfil(xs, ws)
    m = _margem(sub) + float(np.max(ws)) / 2
    comp = float(xs[-1])
    lay = Layout(nome, sub, [poligono_de_perfil(xs, ws)], (0.0, -m, comp, m),
                 portas=[Porta(0.0, 0.0, w_in, "xmin"), Porta(comp, 0.0, w_out, "xmax")], info=info or {},
                 perfil=(xs, ws))
    lay.celula = celula or (l_alim, float(x[-1]))
    i_max, i_min = int(np.argmax(ws)), int(np.argmin(ws))
    lay.cotas = [_cota_h(l_alim, float(x[-1]), float(np.max(ws)) * 0.7, f"L = {x[-1] - l_alim:.3f}"),
                 _cota_v(float(xs[i_max]), -ws[i_max] / 2, ws[i_max] / 2, f"Wmax = {ws[i_max]:.3f}"),
                 _cota_v(float(xs[i_min]), -ws[i_min] / 2, ws[i_min] / 2, f"Wmin = {ws[i_min]:.3f}")]
    return lay


# ---------------------------------------------------------------------------
# Linhas acopladas
# ---------------------------------------------------------------------------
def acoplador_linhas(nome, sub: Substrato, secoes, w0, l_alim=5.0, info=None, capacitores=None, chanfro=False,
                     pontas=False):
    """Acoplador de 4 portas: secoes = [(comprimento, W, S)] em cascata ao longo de x.

    Linha a em cima (portas 1 à esquerda e 2 à direita), linha b embaixo (3 à esquerda,
    4 à direita); as pontas dobram 90° para separar as portas. capacitores: n por seção
    (posições marcadas como cotas; LTA com capacitores SMD/interdigitais entre as linhas);
    pontas: marca os capacitores entre as linhas nas extremidades das seções (Steer 9.10).
    """
    d = max(5 * sub.h_mm, 3 * w0, 2.0)
    x1 = l_alim + w0
    topo, x = [], x1
    for comp, w, s in secoes:
        topo += [_ret(x, s / 2, x + comp, s / 2 + w), _ret(x, -s / 2 - w, x + comp, -s / 2)]
        x += comp
    x2 = x
    c0, c1 = secoes[0], secoes[-1]
    ya0, ya1 = c0[2] / 2 + c0[1] / 2, c1[2] / 2 + c1[1] / 2
    total = x2 + w0 + l_alim
    mf = chanfro_para(w0, sub)[0] if chanfro else None
    topo += caminho([(0.0, ya0 + d), (x1 - w0 / 2, ya0 + d), (x1 - w0 / 2, ya0 - c0[1] / 2 + 1e-9)], w0, mf)
    topo += caminho([(0.0, -ya0 - d), (x1 - w0 / 2, -ya0 - d), (x1 - w0 / 2, -ya0 + c0[1] / 2 - 1e-9)], w0, mf)
    topo += caminho([(total, ya1 + d), (x2 + w0 / 2, ya1 + d), (x2 + w0 / 2, ya1 - c1[1] / 2 + 1e-9)], w0, mf)
    topo += caminho([(total, -ya1 - d), (x2 + w0 / 2, -ya1 - d), (x2 + w0 / 2, -ya1 + c1[1] / 2 - 1e-9)], w0,
                    mf)
    placa = _placa(topo, sub, total)
    portas = [Porta(0.0, ya0 + d, w0, "xmin"), Porta(total, ya1 + d, w0, "xmax"),
              Porta(0.0, -ya0 - d, w0, "xmin"), Porta(total, -ya1 - d, w0, "xmax")]
    lay = Layout(nome, sub, topo, placa, portas=portas, info=info or {})
    lay.celula = (x1 - w0, min(x2, x1 + 3 * secoes[0][0]) + w0)
    comp, w, s = secoes[0]
    lay.cotas = [_cota_h(x1, x1 + comp, s / 2 + w + 0.4 * d, f"ℓ = {comp:.3f}"),
                 _cota_v(x1 + 0.3 * comp, -s / 2, s / 2, f"S = {s:.4f}"),
                 _cota_v(x1 + 0.6 * comp, s / 2, s / 2 + w, f"W = {w:.4f}")]
    if capacitores:
        xc, posicoes = x1, []
        for (comp, w, s), n in zip(secoes, capacitores):
            posicoes += [xc + (k + 0.5) * comp / n for k in range(n)]
            xc += comp
        s0 = secoes[0][2]
        lay.cotas += [_cota_v(x, -s0 / 2, s0 / 2, f"Cc {k + 1} (x = {x:.2f})")
                      for k, x in enumerate(posicoes[:7])]
        lay.info["capacitores_entre_linhas_x_mm"] = ", ".join(f"{x:.3f}" for x in posicoes)
    if pontas:
        xs, xc = [x1], x1
        for comp, _, _ in secoes:
            xc += comp
            xs.append(xc)
        s0 = secoes[0][2]
        lay.cotas += [_cota_v(x, -s0 / 2, s0 / 2, f"C ponta {k + 1}") for k, x in enumerate(xs)]
        lay.info["capacitores_pontas_x_mm"] = ", ".join(f"{x:.3f}" for x in xs)
    return lay


def filtro_linhas_acopladas(nome, sub: Substrato, secoes, w0, l_alim=5.0, info=None):
    """Passa-faixa de linhas acopladas em paralelo (escada): secoes = [(comprimento, W, S)].

    A linha b da seção k continua como linha a da seção k+1 (ressonador λ/2 de x_k a x_(k+2)).
    """
    topo, cotas = [], []
    x, y_a = l_alim, 0.0
    w_first = secoes[0][1]
    topo.append(_ret(0.0, -w0 / 2, l_alim, w0 / 2))
    for k, (comp, w, s) in enumerate(secoes):
        y_b = y_a - (w / 2 + s + w / 2)
        topo.append(_ret(x, y_a - w / 2, x + comp, y_a + w / 2))
        topo.append(_ret(x, y_b - w / 2, x + comp, y_b + w / 2))
        if k < 3:
            cotas.append(_cota_h(x, x + comp, y_a + w / 2 + 0.5, f"ℓ{k + 1} = {comp:.3f}"))
            cotas.append(_cota_v(x + comp * 0.5, y_b + w / 2, y_a - w / 2, f"S{k + 1} = {s:.4f}"))
        x += comp
        y_a = y_b
    total = x + l_alim
    topo.append(_ret(x, y_a - w0 / 2, total, y_a + w0 / 2))
    placa = _placa(topo, sub, total)
    lay = Layout(nome, sub, topo, placa, portas=[Porta(0.0, 0.0, w0, "xmin"), Porta(total, y_a, w0, "xmax")],
                 info=info or {})
    lay.celula = (l_alim - 1, l_alim + sum(c for c, _, _ in secoes[:2]) + 1)
    cotas.append(_cota_v(l_alim + 0.1, -w_first / 2, w_first / 2, f"W1 = {w_first:.4f}"))
    lay.cotas = cotas
    return lay


def filtro_gaps(nome, sub: Substrato, w0, gaps, comprimentos, l_alim=5.0, info=None):
    """Ressonadores série em linha (largura w0) separados por gaps: gaps[N+1], comprimentos[N]."""
    if len(gaps) != len(comprimentos) + 1:
        raise ValueError("São necessários N+1 gaps para N ressonadores.")
    topo, cotas, x = [_ret(0.0, -w0 / 2, l_alim, w0 / 2)], [], l_alim
    for k, s in enumerate(gaps):
        x += s
        if k < len(comprimentos):
            topo.append(_ret(x, -w0 / 2, x + comprimentos[k], w0 / 2))
            if k < 3:
                cotas.append(_cota_h(x, x + comprimentos[k], w0, f"ℓ{k + 1} = {comprimentos[k]:.3f}"))
            x += comprimentos[k]
    total = x + l_alim
    topo.append(_ret(x, -w0 / 2, total, w0 / 2))
    placa = _placa(topo, sub, total)
    lay = Layout(nome, sub, topo, placa, portas=[Porta(0.0, 0.0, w0, "xmin"), Porta(total, 0.0, w0, "xmax")],
                 info=info or {}, perfil=(np.array([0.0, total]), np.array([w0, w0])))
    lay.celula = (l_alim - 1, l_alim + gaps[0] + comprimentos[0] + gaps[1] + 1)
    cotas.append(_cota_h(l_alim, l_alim + gaps[0], -w0, f"gap1 = {gaps[0] * 1e3:.0f} µm"))
    lay.cotas = cotas
    return lay


def ressonador_gap(nome, sub: Substrato, w0, w, comp, s_gap, curto=False, raio_via=0.2, l_alim=5.0,
                   info=None):
    """Linha de alimentação (porta 1), gap e ressonador aberto ou em curto (via no fim)."""
    topo = [_ret(0.0, -w0 / 2, l_alim, w0 / 2), _ret(l_alim + s_gap, -w / 2, l_alim + s_gap + comp, w / 2)]
    vias = []
    x_fim = l_alim + s_gap + comp
    if curto:
        lado = max(2 * raio_via + 0.4, w)
        topo.append(_ret(x_fim, -lado / 2, x_fim + lado, lado / 2))
        vias.append((x_fim + lado / 2, 0.0, raio_via))
        x_fim += lado
    total = x_fim + _margem(sub)
    m = _margem(sub) + max(w, w0) / 2
    lay = Layout(nome, sub, topo, (0.0, -m, total, m), portas=[Porta(0.0, 0.0, w0, "xmin")], info=info or {},
                 vias=vias, perfil=(np.array([0.0, l_alim, l_alim + s_gap, x_fim]), np.array([w0, w0, w, w])))
    lay.celula = (l_alim - 2, min(x_fim, l_alim + 0.3 * comp) + 1)
    lay.cotas = [_cota_h(l_alim, l_alim + s_gap, w0, f"gap = {s_gap * 1e3:.0f} µm"),
                 _cota_h(l_alim + s_gap, l_alim + s_gap + comp, -w, f"ℓ = {comp:.3f}")]
    return lay


# ---------------------------------------------------------------------------
# Híbridas e divisores
# ---------------------------------------------------------------------------
def branch_line(nome, sub: Substrato, w0, w_s, l_s, w_p, l_p, l_alim=5.0, info=None, sw_s=None, sw_p=None):
    """Híbrida de quadratura: braços série (horizontais, w_s × l_s) e em derivação (verticais,
    w_p × l_p), medidos entre linhas de centro. sw_s/sw_p: (GeometriaSW, n) para braços slow-wave
    (stubs voltados para fora do quadrado).

    Portas: 1 (sup. esq.), 2 (sup. dir.), 3 (inf. dir.), 4 (inf. esq.), como em dispositivos.branch_line.
    """
    xl, yt = l_alim + (sw_p[0].l_carga + sw_p[0].w_linha if sw_p and sw_p[0].tipo != "Patches" else 0.0), 0.0
    if sw_s:
        l_s = sw_s[1] * sw_s[0].periodo
    if sw_p:
        l_p = sw_p[1] * sw_p[0].periodo
    xr, yb = xl + l_s, yt - l_p
    topo = []
    for y, lado in ((yt, 1), (yb, -1)):
        if sw_s:
            topo += braco_slow_wave((xl, y), "x", sw_s[0], sw_s[1], lado)
        else:
            topo.append(_ret(xl, y - w_s / 2, xr, y + w_s / 2))
    for x, lado in ((xl, -1), (xr, 1)):
        if sw_p:
            topo += braco_slow_wave((x, yt), "-y", sw_p[0], sw_p[1], lado)
        else:
            topo.append(_ret(x - w_p / 2, yb, x + w_p / 2, yt))
    lado = max(w_s if not sw_s else sw_s[0].w_linha, w_p if not sw_p else sw_p[0].w_linha)
    for x in (xl, xr):                                   # cantos
        for y in (yt, yb):
            topo.append(_ret(x - lado / 2, y - lado / 2, x + lado / 2, y + lado / 2))
    total = xr + l_alim
    for y in (yt, yb):
        topo.append(_ret(0.0, y - w0 / 2, xl, y + w0 / 2))
        topo.append(_ret(xr, y - w0 / 2, total, y + w0 / 2))
    placa = _placa(topo, sub, total)
    portas = [Porta(0.0, yt, w0, "xmin"), Porta(total, yt, w0, "xmax"), Porta(total, yb, w0, "xmax"),
              Porta(0.0, yb, w0, "xmin")]
    lay = Layout(nome, sub, topo, placa, portas=portas, info=info or {})
    lay.celula = (l_alim - 1, total - l_alim + 1)
    lay.cotas = [_cota_h(xl, xr, yt + lado, f"ℓ série = {l_s:.3f}"),
                 _cota_v(xl + lado, yb, yt, f"ℓ derivação = {l_p:.3f}")]
    if not sw_s:
        lay.cotas.append(_cota_v(0.5 * (xl + xr), yt - w_s / 2, yt + w_s / 2, f"Ws = {w_s:.3f}"))
    if not sw_p:
        lay.cotas.append(_cota_h(xr - w_p / 2, xr + w_p / 2, 0.5 * (yt + yb), f"Wp = {w_p:.3f}"))
    return lay


def rat_race(nome, sub: Substrato, w0, w_r, raio, l_alim=5.0, info=None, n_seg=144, chanfro=False):
    """Anel de 1,5λg (raio na linha de centro) com portas 1 (150°), 2 (90°), 4 (30°) e 3 (-30°),
    na ordem de dispositivos.rat_race (1-2, 2-4 e 4-3 com 90°; 3-1 com 270°)."""
    r1, r2 = raio - w_r / 2, raio + w_r / 2
    ang = np.linspace(0, 2 * np.pi, n_seg + 1)
    xc = l_alim + r2 + w0
    topo = []
    for a, b in zip(ang[:-1], ang[1:]):
        topo.append(np.array([[xc + r1 * np.cos(a), r1 * np.sin(a)], [xc + r2 * np.cos(a), r2 * np.sin(a)],
                              [xc + r2 * np.cos(b), r2 * np.sin(b)], [xc + r1 * np.cos(b), r1 * np.sin(b)]]))
    d = max(3 * w0, 2.0)
    total = xc + r2 + w0 + l_alim
    p = {k: (xc + raio * np.cos(np.radians(g)), raio * np.sin(np.radians(g)))
         for k, g in ((1, 150), (2, 90), (4, 30), (3, -30))}
    topo += caminho([(0.0, p[1][1]), (p[1][0], p[1][1])], w0)
    topo += caminho([(p[2][0], p[2][1]), (p[2][0], raio + d + w0), (0.0, raio + d + w0)], w0,
                    chanfro_para(w0, sub)[0] if chanfro else None)
    topo += caminho([(p[4][0], p[4][1]), (total, p[4][1])], w0)
    topo += caminho([(p[3][0], p[3][1]), (total, p[3][1])], w0)
    placa = _placa(topo, sub, total)
    portas = [Porta(0.0, p[1][1], w0, "xmin"), Porta(0.0, raio + d + w0, w0, "xmin"),
              Porta(total, p[3][1], w0, "xmax"), Porta(total, p[4][1], w0, "xmax")]
    lay = Layout(nome, sub, topo, placa, portas=portas, info=info or {})
    lay.celula = (xc - r2 - 1, xc + r2 + 1)
    lay.cotas = [_cota_h(xc - raio, xc + raio, 0.0, f"2r = {2 * raio:.3f}"),
                 _cota_v(xc, -r2, -r1, f"Wr = {w_r:.3f}")]
    return lay


def wilkinson(nome, sub: Substrato, w0, w2, l2, w3, l3, wo2=None, wo3=None, g_r=1.0, t2=None, t3=None,
              l_alim=5.0, info=None, chanfro_f0=None):
    """Wilkinson em 'C': braço 2 em cima e 3 embaixo (wk × lk na linha de centro), resistor entre
    as pontas (gap g_r, borda a borda). wo2/wo3: largura das linhas de saída (Z0, ou R2/R3 no
    desigual); t2/t3: (W, ℓ) dos transformadores de saída λ/4 para Z0, ou None.
    chanfro_f0 [Hz]: chanfra as curvas (eq. 3.115) e alonga cada braço de 2·Δb (eq. 3.116), pois
    cada curva chanfrada é Δb mais curta que as linhas de centro."""
    wo2, wo3 = wo2 or w0, wo3 or w0
    xj = l_alim
    la = min(l2, l3) / 3
    xr = xj + la
    topo = [_ret(0.0, -w0 / 2, xj, w0 / 2)]
    cotas, saidas = [], []
    for sinal, w, l, wo, t in ((1, w2, l2, wo2, t2), (-1, w3, l3, wo3, t3)):
        mf, db = chanfro_para(w, sub, chanfro_f0) if chanfro_f0 else (None, 0.0)
        l = l + 2 * db
        yr = sinal * (g_r / 2 + w / 2)
        a = (l - la + abs(yr)) / 2                    # a + la + (a - |yr|) = l
        if a <= abs(yr) + w:
            raise ValueError("Braço curto demais para a topologia em 'C': reduza o gap do resistor.")
        topo += caminho([(xj, 0.0), (xj, sinal * a), (xr, sinal * a), (xr, yr)], w, mf)
        saidas.append((sinal, yr, a, wo, t))
        cotas.append(_cota_v(xj - w, 0.0, sinal * a, f"a{2 if sinal > 0 else 3} = {a:.3f}"))
    # saídas: logo após o resistor (borda interna alinhada ao gap) e afastando-se do eixo
    wo_max = max(wo2, wo3)
    x_t = xr + max(w2, w3) / 2 + wo_max / 2 + max(1.0, 2 * sub.h_mm)
    yv = g_r / 2 + wo_max + max(3 * wo_max, 2.0)
    fins = []
    for sinal, yr, a, wo, t in saidas:
        yo, y = sinal * (g_r / 2 + wo / 2), sinal * yv
        topo += caminho([(xr, yo), (x_t, yo), (x_t, y), (x_t + 1.0, y)], wo,
                        chanfro_para(wo, sub)[0] if chanfro_f0 else None)
        x = x_t + 1.0
        if t:
            topo.append(_ret(x, y - t[0] / 2, x + t[1], y + t[0] / 2))
            cotas.append(_cota_h(x, x + t[1], y + sinal * t[0], f"λ/4 = {t[1]:.3f}"))
            x += t[1]
        fins.append((x, y, w0 if t else wo))
    total = max(x for x, _, _ in fins) + l_alim
    for x, y, w in fins:
        topo.append(_ret(x, y - w / 2, total, y + w / 2))
    placa = _placa(topo, sub, total)
    portas = [Porta(0.0, 0.0, w0, "xmin"), Porta(total, fins[0][1], fins[0][2], "xmax"),
              Porta(total, fins[1][1], fins[1][2], "xmax")]
    lay = Layout(nome, sub, topo, placa, portas=portas, info=info or {})
    lay.celula = (xj - 2, x_t + 3)
    ys = max(s[2] for s in saidas)
    cotas += [_cota_h(xj, xr, ys + max(w2, w3), f"c = {la:.3f}"),
              _cota_v(xr + max(w2, w3), -g_r / 2, g_r / 2, f"R: gap {g_r:.2f}")]
    lay.cotas = cotas
    return lay


def juncao_t(nome, sub: Substrato, w0, w2, w3, t2=None, t3=None, l_alim=5.0, info=None, chanfro=False):
    """Junção T: entrada à esquerda, saídas para cima/baixo e depois para a direita.

    As saídas são linhas de Z2/Z3 (larguras w2, w3); com t2/t3 = (W, ℓ), segue o λ/4 e a linha
    de Z0 (larguras de porta w0)."""
    xj = l_alim
    a = max(3 * max(w0, w2, w3), 3.0)
    topo = [_ret(0.0, -w0 / 2, xj, w0 / 2)]
    fins, cotas = [], []
    for sinal, w, t in ((1, w2, t2), (-1, w3, t3)):
        y = sinal * a
        topo += caminho([(xj, 0.0), (xj, y), (xj + a, y)], w, chanfro_para(w, sub)[0] if chanfro else None)
        x = xj + a
        if t:
            topo.append(_ret(x, y - t[0] / 2, x + t[1], y + t[0] / 2))
            cotas.append(_cota_h(x, x + t[1], y + sinal * t[0], f"λ/4 = {t[1]:.3f}"))
            x += t[1]
        fins.append((x, y, w0 if t else w))
    total = max(x for x, _, _ in fins) + l_alim
    for x, y, w in fins:
        topo.append(_ret(x, y - w / 2, total, y + w / 2))
    placa = _placa(topo, sub, total)
    portas = [Porta(0.0, 0.0, w0, "xmin"), Porta(total, fins[0][1], fins[0][2], "xmax"),
              Porta(total, fins[1][1], fins[1][2], "xmax")]
    lay = Layout(nome, sub, topo, placa, portas=portas, info=info or {})
    lay.celula = (xj - 2, min(total, xj + a + 5))
    lay.cotas = cotas + [_cota_v(xj + max(w2, w3), 0.0, a, f"a = {a:.3f}"),
                         _cota_v(xj + a - 1, a - w2 / 2, a + w2 / 2, f"W2 = {w2:.3f}"),
                         _cota_v(xj + a - 1, -a - w3 / 2, -a + w3 / 2, f"W3 = {w3:.3f}")]
    return lay


# ---------------------------------------------------------------------------
# Garg, Bahl e Bianchi: DGS, slotline, fixture de conectores e balun de Marchand
# ---------------------------------------------------------------------------
def haltere(xc, a, b, g, l_pescoco):
    """Fendas no terra de uma DGS em haltere centrada em (xc, 0): cabeças a (x) × b (y) e pescoço de
    largura g (x) e comprimento l_pescoco (y) sob a linha (Fig. 6.1)."""
    yc = l_pescoco / 2 + b / 2
    return [_ret(xc - a / 2, yc - b / 2, xc + a / 2, yc + b / 2),
            _ret(xc - a / 2, -yc - b / 2, xc + a / 2, -yc + b / 2),
            _ret(xc - g / 2, -l_pescoco / 2 - 1e-6, xc + g / 2, l_pescoco / 2 + 1e-6)]


def linha_dgs(nome, sub: Substrato, w, n, periodo, a, b, g, l_pescoco, l_alim=5.0, stubs=None, info=None):
    """Linha de largura w com n DGS em haltere a cada período; stubs: [(x, comprimento, largura)]
    de pares simétricos de stubs abertos (passa-baixas da seção 6.4.1), x ao longo da linha."""
    comp = 2 * l_alim + n * periodo
    topo = [_ret(0.0, -w / 2, comp, w / 2)]
    fendas = []
    for k in range(int(n)):
        fendas += haltere(l_alim + (k + 0.5) * periodo, a, b, g, l_pescoco)
    y_max = l_pescoco / 2 + b
    for x, ls, ws in (stubs or []):
        xs = l_alim + x
        topo += [_ret(xs - ws / 2, w / 2 - 1e-6, xs + ws / 2, w / 2 + ls), _ret(xs - ws / 2, -w / 2 - ls,
                                                                              xs + ws / 2, -w / 2 + 1e-6)]
        y_max = max(y_max, w / 2 + ls)
    m = y_max + _margem(sub)
    lay = Layout(nome, sub, topo, (0.0, -m, comp, m), portas=[Porta(0.0, 0.0, w, "xmin"), Porta(comp, 0.0, w, "xmax")],
                 info=info or {}, perfil=(np.array([0.0, comp]), np.array([w, w])), fendas_terra=fendas)
    lay.celula = (l_alim, l_alim + periodo)
    xc = l_alim + 0.5 * periodo
    lay.cotas = [_cota_h(xc - a / 2, xc + a / 2, l_pescoco / 2 + b + 0.3, f"a = {a:.3f}"),
                 _cota_v(xc + a / 2 + 0.3, l_pescoco / 2, l_pescoco / 2 + b, f"b = {b:.3f}"),
                 _cota_h(xc - g / 2, xc + g / 2, -l_pescoco / 2 - 0.3, f"g = {g:.3f}"),
                 _cota_v(xc - a / 2 - 0.3, -l_pescoco / 2, l_pescoco / 2, f"pescoço = {l_pescoco:.3f}")]
    if n > 1:
        lay.cotas.append(_cota_h(l_alim, l_alim + periodo, -m + _margem(sub) / 2, f"p = {periodo:.3f}"))
    return lay


def slotline(nome, sub: Substrato, w_fenda, comp, largura_plano=10.0, l_alim=0.0, info=None):
    """Slotline: fenda de largura w_fenda entre dois planos metálicos no topo, sem terra inferior."""
    total = comp + 2 * l_alim
    topo = [_ret(0.0, w_fenda / 2, total, w_fenda / 2 + largura_plano),
            _ret(0.0, -w_fenda / 2 - largura_plano, total, -w_fenda / 2)]
    m = w_fenda / 2 + largura_plano + max(3 * sub.h_mm, 1.0)
    lay = Layout(nome, sub, topo, (0.0, -m, total, m), portas=[Porta(0.0, 0.0, w_fenda, "xmin"),
                                                               Porta(total, 0.0, w_fenda, "xmax")],
                 info=info or {}, sem_terra=True)
    lay.celula = (0.0, min(total, 10 * w_fenda + 5))
    lay.cotas = [_cota_v(1.0, -w_fenda / 2, w_fenda / 2, f"fenda = {w_fenda:.4f}")]
    return lay


def fixture_conector(nome, sub: Substrato, w0, comp, trecho=None, info=None):
    """Linha de teste entre dois conectores de borda: [trecho de compensação (W, ℓ)], linha de Z0
    de comprimento comp, [trecho]. As portas ficam nas bordas, onde pousa o pino do conector."""
    tr = [(trecho[1], trecho[0])] if trecho else []
    x, w = perfil_degraus_local(tr + [(comp, w0)] + tr)
    w_ponta = trecho[0] if trecho else w0
    lay = linha_perfil(nome, sub, x, w, w_ponta, w_ponta, 0.0, info=info,
                       celula=(0.0, min(x[-1], (trecho[1] if trecho else 0) + 4 * w0)))
    if trecho:
        lay.cotas = [_cota_h(0.0, trecho[1], trecho[0] * 2, f"compensação ℓ = {trecho[1]:.3f}"),
                     _cota_v(trecho[1] / 2, -trecho[0] / 2, trecho[0] / 2, f"Wc = {trecho[0]:.3f}"),
                     _cota_v(trecho[1] + w0, -w0 / 2, w0 / 2, f"W0 = {w0:.3f}")]
    return lay


def perfil_degraus_local(trechos):
    xs, ws, x = [], [], 0.0
    for comp, larg in trechos:
        xs += [x, x + comp]
        ws += [larg, larg]
        x += comp
    return np.array(xs), np.array(ws)


def balun_marchand(nome, sub: Substrato, w0, w, s, comp, w_out, raio_via=0.2, l_alim=5.0, info=None, chanfro=False,
                   capacitores=None):
    """Marchand com duas seções acopladas em linha: linha a (entrada à esquerda, ponta aberta à
    direita); linhas b com pontas externas aterradas por vias e saídas nas pontas internas, que
    descem e vão para as bordas (porta 2 à esquerda, porta 3 à direita).

    capacitores: {nó: C [F]} do Marchand miniaturizado (1 meio da linha a, 2 ponta aberta, 3 e 4
    saídas), marcados como cotas nas posições dos capacitores ao terra."""
    ya, yb = s / 2 + w / 2, -(s / 2 + w / 2)
    gb = max(2 * w, 0.5)                        # separação entre as linhas b das duas seções
    lado = max(2 * raio_via + 0.4, w)
    x1 = l_alim + lado
    x2 = x1 + comp
    x3 = x2 + gb + comp
    total = x3 + lado + l_alim
    d = max(4 * w_out, 3 * sub.h_mm, 2.0)
    mf = chanfro_para(w_out, sub)[0] if chanfro else None
    topo = [_ret(0.0, ya - w0 / 2, x1, ya + w0 / 2), _ret(x1, ya - w / 2, x3, ya + w / 2),
            _ret(x1, yb - w / 2, x2, yb + w / 2), _ret(x2 + gb, yb - w / 2, x3, yb + w / 2),
            _ret(x1 - lado, yb - lado / 2, x1, yb + lado / 2), _ret(x3, yb - lado / 2, x3 + lado, yb + lado / 2)]
    vias = [(x1 - lado / 2, yb, raio_via), (x3 + lado / 2, yb, raio_via)]
    y_out = yb - d
    topo += caminho([(x2 - w_out / 2, yb), (x2 - w_out / 2, y_out), (0.0, y_out)], w_out, mf)
    topo += caminho([(x2 + gb + w_out / 2, yb), (x2 + gb + w_out / 2, y_out), (total, y_out)], w_out, mf)
    placa = _placa(topo, sub, total)
    lay = Layout(nome, sub, topo, placa, portas=[Porta(0.0, ya, w0, "xmin"), Porta(0.0, y_out, w_out, "xmin"),
                                                 Porta(total, y_out, w_out, "xmax")], info=info or {}, vias=vias)
    lay.celula = (x1 - lado, x1 + min(comp, 6 * (w + s)))
    lay.cotas = [_cota_h(x1, x2, ya + w, f"ℓ = {comp:.3f}"), _cota_v(x1 + 1, yb + w / 2, ya - w / 2, f"S = {s:.4f}"),
                 _cota_v(x1 + 2, ya - w / 2, ya + w / 2, f"W = {w:.4f}")]
    if capacitores:
        pos = {1: (x2 + gb / 2, ya), 2: (x3, ya), 3: (x2, yb), 4: (x2 + gb, yb)}
        for no, c in sorted(capacitores.items()):
            if c > 0 and no in pos:
                xc, yc = pos[no]
                lay.cotas.append(_cota_v(xc, yc - w / 2, yc + w / 2, f"C{no} = {c * 1e12:.3f} pF ao terra"))
        lay.info["capacitores_terra_pF"] = ", ".join(f"no {k}: {c * 1e12:.4f}" for k, c in sorted(capacitores.items()))
    return lay


def combline(nome, sub: Substrato, ws, ss, comp, w0, caps, cb, raio_via=0.2, l_alim=5.0, info=None):
    """Filtro combline (Steer 10.14): N ressonadores verticais de largura ws[k] e espaçamentos ss[k],
    em curto embaixo (via) e com capacitor ao terra no topo; entrada e saída por capacitores série
    Cb entre as linhas de alimentação (horizontais, acima) e o topo do 1º e do último ressonador."""
    n = len(ws)
    lado = max(2 * raio_via + 0.4, max(ws))
    gap_cb = max(0.5, 2 * sub.h_mm / 3)                 # espaço do capacitor série (pad SMD)
    d = max(3 * w0, 2.0)
    x = l_alim + lado
    xs = []
    for k in range(n):
        xs.append(x)
        x += ws[k] + (ss[k] if k < n - 1 else 0.0)
    total = x + lado + l_alim
    topo, vias = [], []
    for k in range(n):
        xa, xb = xs[k], xs[k] + ws[k]
        xc = 0.5 * (xa + xb)
        topo.append(_ret(xa, 0.0, xb, comp))
        topo.append(_ret(xc - lado / 2, -lado, xc + lado / 2, 0.0))
        vias.append((xc, -lado / 2, raio_via))
    y_f = comp + gap_cb + d
    xc1, xcn = xs[0] + ws[0] / 2, xs[-1] + ws[-1] / 2
    topo += caminho([(0.0, y_f), (xc1, y_f), (xc1, comp + gap_cb)], w0)
    topo += caminho([(total, y_f), (xcn, y_f), (xcn, comp + gap_cb)], w0)
    placa = _placa(topo, sub, total)
    lay = Layout(nome, sub, topo, placa, portas=[Porta(0.0, y_f, w0, "xmin"), Porta(total, y_f, w0, "xmax")],
                 info=info or {}, vias=vias)
    lay.celula = (xs[0], xs[min(1, n - 1)] + ws[min(1, n - 1)])
    lay.cotas = [_cota_v(xs[0] - 0.4, 0.0, comp, f"ℓ = {comp:.3f}"),
                 _cota_v(xc1 + w0, comp, comp + gap_cb, f"Cb = {cb * 1e12:.3f} pF série"),
                 _cota_v(xcn + w0, comp, comp + gap_cb, f"Cb = {cb * 1e12:.3f} pF série")]
    for k in range(n):
        lay.cotas.append(_cota_h(xs[k], xs[k] + ws[k], comp + 0.3 * gap_cb, f"C{k + 1} = {caps[k] * 1e12:.3f} pF"))
        if k < n - 1:
            lay.cotas.append(_cota_h(xs[k] + ws[k], xs[k + 1], 0.5 * comp, f"S = {ss[k]:.3f}"))
    lay.info["capacitores_topo_terra_pF"] = ", ".join(f"{c * 1e12:.4f}" for c in caps)
    lay.info["Cb_serie_pF"] = f"{cb * 1e12:.4f}"
    return lay
