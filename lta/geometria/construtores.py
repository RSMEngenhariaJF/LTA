"""Construção de layouts a partir dos projetos calculados nas abas."""
import numpy as np

from ..core.microstrip import Substrato, sintese
from .layout import (Layout, Porta, margem_lateral, perfil_degraus,
                     poligono_de_perfil, simplificar_perfil)


def _montar(nome, sub, x, w, w_porta, l_alim, furos=(), info=None, meia_larg=None):
    """Acrescenta linhas de alimentação nas pontas e monta o Layout."""
    x = np.asarray(x) + l_alim
    if l_alim <= 0:
        xs, ws = np.asarray(x), np.asarray(w)
    else:
        # degraus verticais entre alimentação e estrutura, se as larguras diferirem
        xs = np.concatenate([[0.0, l_alim], x, [x[-1], x[-1] + l_alim]])
        ws = np.concatenate([[w_porta, w_porta], w, [w_porta, w_porta]])
    xs, ws = simplificar_perfil(xs, ws)
    poli = poligono_de_perfil(xs, ws)
    comp = float(xs[-1])
    if meia_larg is None:
        meia_larg = margem_lateral(sub, float(np.max(ws)))
    placa = (0.0, -meia_larg, comp, meia_larg)
    portas = [Porta(0.0, 0.0, float(ws[0]), "xmin"), Porta(comp, 0.0, float(ws[-1]), "xmax")]
    return Layout(nome, sub, [poli], placa, list(furos), portas, info or {}, perfil=(xs, ws))


def _cota_h(x1, x2, y, texto):
    return (x1, y, x2, y, texto)


def _cota_v(x, w, texto):
    return (x, -w / 2, x, w / 2, texto)


def linha_carga_capacitiva(geo, n_celulas, sub: Substrato, l_alim=5.0, nome="carga_capacitiva", info=None):
    """Linha de alta impedância com patches (trechos largos) centrados em cada célula."""
    if geo.l_patch >= geo.periodo:
        raise ValueError(f"Patch ({geo.l_patch:.3f} mm) maior que o período "
                         f"({geo.periodo:.3f} mm). Reduza a impedância do patch.")
    resto = (geo.periodo - geo.l_patch) / 2
    trechos = []
    for _ in range(int(n_celulas)):
        trechos += [(resto, geo.w_linha), (geo.l_patch, geo.w_patch), (resto, geo.w_linha)]
    x, w = perfil_degraus(trechos)
    lay = _montar(nome, sub, x, w, geo.w_porta, l_alim, info=info)
    c0, d = l_alim, 0.25 * geo.w_patch
    lay.celula = (c0, c0 + geo.periodo)
    lay.cotas = [
        _cota_h(c0, c0 + geo.periodo, geo.w_patch / 2 + d, f"l = {geo.periodo:.3f}"),
        _cota_h(c0 + resto, c0 + resto + geo.l_patch, -geo.w_patch / 2 - d, f"lp = {geo.l_patch:.3f}"),
        _cota_v(c0 + resto + 0.8 * geo.l_patch, geo.w_patch, f"Wp = {geo.w_patch:.3f}"),
        _cota_v(c0 + resto / 2, geo.w_linha, f"W = {geo.w_linha:.3f}"),
    ]
    if l_alim > 0:
        lay.cotas.append(_cota_v(l_alim / 2, geo.w_porta, f"W porta = {geo.w_porta:.3f}"))
    return lay


def linha_largura_modulada(z_mm, w_mm, sub: Substrato, l_alim=5.0, nome="ebg_largura", info=None, periodo=None):
    """EBG com modulação da largura da trilha (perfil contínuo W(z))."""
    lay = _montar(nome, sub, z_mm, w_mm, float(w_mm[0]), l_alim, info=info)
    if periodo:
        z = np.asarray(z_mm)
        na_cel = z <= z[0] + periodo
        zc, wc = z[na_cel] + l_alim, np.asarray(w_mm)[na_cel]
        i_max, i_min = int(np.argmax(wc)), int(np.argmin(wc))
        lay.celula = (l_alim, l_alim + periodo)
        lay.cotas = [
            _cota_h(l_alim, l_alim + periodo, wc.max() / 2 + 0.3 * wc.max(), f"l = {periodo:.3f}"),
            _cota_v(zc[i_max], wc[i_max], f"Wmax = {wc[i_max]:.3f}"),
            _cota_v(zc[i_min], wc[i_min], f"Wmin = {wc[i_min]:.3f}"),
        ]
    else:  # perfil sem período: detalhe em torno do trecho mais largo (±10% do comprimento)
        z = np.asarray(z_mm) + l_alim
        w = np.asarray(w_mm)
        i_max = int(np.argmax(w))
        meia = 0.1 * (z[-1] - z[0])
        x0, x1 = max(z[0], z[i_max] - meia), min(z[-1], z[i_max] + meia)
        janela = (z >= x0) & (z <= x1)
        i_min = int(np.flatnonzero(janela)[np.argmin(w[janela])])
        lay.celula = (float(x0), float(x1))
        lay.cotas = [_cota_v(z[i_max], w[i_max], f"Wmax = {w[i_max]:.3f}"),
                     _cota_v(z[i_min], w[i_min], f"W local mín = {w[i_min]:.3f}"),
                     _cota_h(x0, x1, w.max() / 2 + 0.3 * w.max(), f"trecho = {x1 - x0:.2f}")]
    return lay


def linha_furos_terra(periodo, raio, n_furos, z_ref, sub: Substrato, l_alim=5.0,
                      nome="ebg_furos", info=None):
    """Linha uniforme com furos circulares no plano de terra (Fig. 2.15, versão 1D)."""
    w = sintese(z_ref, sub)
    comp = periodo * n_furos
    x = np.array([0.0, comp])
    ws = np.array([w, w])
    furos = [(l_alim + periodo * (i + 0.5), 0.0, raio) for i in range(int(n_furos))]
    meia = max(margem_lateral(sub, w), raio + 2 * sub.h_mm)
    lay = _montar(nome, sub, x, ws, w, l_alim, furos=furos, info=info, meia_larg=meia)
    c0 = furos[0][0]
    lay.celula = (l_alim, l_alim + periodo)
    lay.cotas = [
        _cota_v(c0, 2 * raio, f"Ø = {2 * raio:.3f}"),
        _cota_v(l_alim + 0.08 * periodo, w, f"W = {w:.3f}"),
    ]
    if n_furos > 1:
        lay.cotas.append(_cota_h(c0, c0 + periodo, raio * 1.25, f"l = {periodo:.3f}"))
    return lay


def _ret(x0, y0, x1, y1):
    return np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], dtype=float)


def _interdigital(xc, n_dedos, comp_dedo, larg_dedo):
    """Capacitor interdigital centrado em xc: duas barras + n dedos alternados (largura = espaço).

    Retorna (polígonos, x inicial, x final).
    """
    larg_tot = (2 * n_dedos - 1) * larg_dedo
    ext = comp_dedo + 3 * larg_dedo                     # dedos + folga na ponta + duas barras
    xa, xb = xc - ext / 2, xc + ext / 2
    y0 = -larg_tot / 2
    polis = [_ret(xa, y0, xa + larg_dedo, -y0), _ret(xb - larg_dedo, y0, xb, -y0)]
    for k in range(n_dedos):
        ya = y0 + 2 * k * larg_dedo
        if k % 2 == 0:
            polis.append(_ret(xa + larg_dedo, ya, xa + larg_dedo + comp_dedo, ya + larg_dedo))
        else:
            polis.append(_ret(xb - larg_dedo - comp_dedo, ya, xb - larg_dedo, ya + larg_dedo))
    return polis, xa, xb


def linha_csrr(n_celulas, periodo, w_linha, geo, sub: Substrato, serie="Gap simples", gap=0.2,
               n_dedos=6, comp_dedo=1.0, larg_dedo=0.1, l_alim=5.0, nome="linha_csrr", info=None):
    """Microstrip com elemento série (gap ou interdigital) e CSRR no terra no centro de cada célula.

    geo: GeometriaCSRR (lta.core.csrr). Fig. 3.37 / 3.42 do livro (CSRR + gap / interdigital).
    """
    from ..core.csrr import fendas_csrr
    if periodo < 2 * geo.r_ext:
        raise ValueError(f"Período ({periodo:.3f} mm) menor que o diâmetro do CSRR ({2 * geo.r_ext:.3f} mm): "
                         "os anéis de células vizinhas se sobrepõem.")
    comp = 2 * l_alim + n_celulas * periodo
    topo, fendas, x_esq = [], [], 0.0
    if serie == "Interdigital":
        if n_dedos < 2:
            raise ValueError("Use ao menos 2 dedos.")
        extensao = comp_dedo + 3 * larg_dedo
    else:
        extensao = gap
    if extensao >= periodo:
        raise ValueError("O elemento série não cabe no período.")
    for i in range(n_celulas):
        xc = l_alim + (i + 0.5) * periodo
        fendas += fendas_csrr(xc, 0.0, geo)
        xa, xb = xc - extensao / 2, xc + extensao / 2
        topo.append(_ret(x_esq, -w_linha / 2, xa, w_linha / 2))
        if serie == "Interdigital":
            topo += _interdigital(xc, n_dedos, comp_dedo, larg_dedo)[0]
        x_esq = xb
    topo.append(_ret(x_esq, -w_linha / 2, comp, w_linha / 2))
    meia = max(margem_lateral(sub, w_linha), geo.r_ext + 2 * sub.h_mm)
    portas = [Porta(0.0, 0.0, w_linha, "xmin"), Porta(comp, 0.0, w_linha, "xmax")]
    lay = Layout(nome, sub, topo, (0.0, -meia, comp, meia), portas=portas, info=info or {},
                 perfil=(np.array([0.0, comp]), np.array([w_linha, w_linha])), fendas_terra=fendas)
    xc = l_alim + 0.5 * periodo
    lay.celula = (l_alim, l_alim + periodo)
    lay.cotas = [_cota_h(l_alim, l_alim + periodo, geo.r_ext * 1.15, f"l = {periodo:.3f}"),
                 _cota_v(l_alim + 0.06 * periodo, w_linha, f"W = {w_linha:.3f}"),
                 _cota_v(xc + 0.0, 2 * geo.r_ext, f"2·rext = {2 * geo.r_ext:.3f}"),
                 _cota_h(xc - extensao / 2, xc + extensao / 2, -geo.r_ext * 1.1,
                         f"{'gap' if serie != 'Interdigital' else 'interdigital'} = {extensao:.3f}")]
    return lay


def linha_cl_loaded(geo, n_celulas, sub: Substrato, l_alim=5.0, nome="crlh_cl_loaded", info=None):
    """CRLH microstrip CL-loaded (Fig. 3.27): interdigitais em série e stubs com via em derivação.

    geo: GeometriaCL (lta.core.cl_loaded). Células T: capacitores nas fronteiras (2CL nas
    pontas, CL entre células) e stub no centro de cada célula.
    """
    p, w = geo.periodo, geo.w_linha
    if geo.extensao_idc >= p:
        raise ValueError(f"Interdigital ({geo.extensao_idc:.2f} mm) não cabe no período ({p:.2f} mm).")
    lado_pad = 2 * geo.raio_via + 0.4
    x0 = l_alim + (geo.comp_dedo_ext + 3 * geo.larg_dedo) / 2   # capacitores das pontas (2CL) após a alimentação
    comp = 2 * x0 + n_celulas * p
    topo, vias, x_esq = [], [], 0.0
    for i in range(n_celulas + 1):
        xb_c = x0 + i * p
        comp_d = geo.comp_dedo_ext if i in (0, n_celulas) else geo.comp_dedo
        polis, xa, xb = _interdigital(xb_c, geo.n_dedos, comp_d, geo.larg_dedo)
        topo.append(_ret(x_esq, -w / 2, xa, w / 2))
        topo += polis
        x_esq = xb
        if i < n_celulas:                               # stub + pad + via no centro da célula
            xc = x0 + (i + 0.5) * p
            y1 = w / 2 + geo.comp_stub
            topo.append(_ret(xc - geo.w_stub / 2, 0.0, xc + geo.w_stub / 2, y1))
            topo.append(_ret(xc - lado_pad / 2, y1, xc + lado_pad / 2, y1 + lado_pad))
            vias.append((xc, y1 + lado_pad / 2, geo.raio_via))
    topo.append(_ret(x_esq, -w / 2, comp, w / 2))
    margem = max(5 * sub.h_mm, 2.0)
    y_max = w / 2 + geo.comp_stub + lado_pad + margem
    y_min = -max(margem_lateral(sub, w), (2 * geo.n_dedos - 1) * geo.larg_dedo / 2 + margem)
    portas = [Porta(0.0, 0.0, w, "xmin"), Porta(comp, 0.0, w, "xmax")]
    lay = Layout(nome, sub, topo, (0.0, y_min, comp, y_max), portas=portas, info=info or {},
                 perfil=(np.array([0.0, comp]), np.array([w, w])), vias=vias)
    xc = x0 + 0.5 * p
    lay.celula = (x0, x0 + p)
    lay.cotas = [_cota_h(x0, x0 + p, -w * 1.5, f"p = {p:.3f}"),
                 _cota_v(xc + geo.w_stub, 2 * (w / 2 + geo.comp_stub), f"stub = {geo.comp_stub:.3f}"),
                 _cota_h(x0 + p - geo.extensao_idc / 2, x0 + p + geo.extensao_idc / 2, w * 1.2,
                         f"IDC = {geo.extensao_idc:.3f}"),
                 _cota_v(x0 + 0.3 * p, w, f"W = {w:.3f}")]
    return lay


def linha_siw(p, comprimento, sub: Substrato, w_porta, w_cunha, l_cunha, l_alim=5.0, nome="siw", info=None):
    """SIW (fileiras de vias, Fig. 6.33) com transições em cunha para microstrip nas duas pontas.

    p: ProjetoSIW (lta.core.siw). Medidas em mm, exceto p (em m). A cunha vai de w_porta
    (linha de 50 Ω) a w_cunha na entrada do guia; suas dimensões são ponto de partida e
    devem ser otimizadas em EM (livro, p. 465).
    """
    w = p.w * 1e3
    d, s = p.d * 1e3, p.s * 1e3
    if w_cunha > w or w_porta > w_cunha:
        raise ValueError("A cunha deve ir da largura da porta até, no máximo, a largura do guia.")
    x1 = l_alim + l_cunha                 # início do guia
    x2 = x1 + comprimento
    total = x2 + l_cunha + l_alim
    meia_placa = w / 2 + d / 2 + 0.1      # cobre cobre as vias
    topo = [_ret(0.0, -w_porta / 2, l_alim, w_porta / 2),
            np.array([[l_alim, -w_porta / 2], [x1, -w_cunha / 2], [x1, w_cunha / 2], [l_alim, w_porta / 2]]),
            _ret(x1, -meia_placa, x2, meia_placa),
            np.array([[x2, -w_cunha / 2], [x2 + l_cunha, -w_porta / 2], [x2 + l_cunha, w_porta / 2],
                      [x2, w_cunha / 2]]),
            _ret(x2 + l_cunha, -w_porta / 2, total, w_porta / 2)]
    n_vias = int(np.floor((comprimento - d) / s)) + 1
    xs = x1 + d / 2 + np.arange(n_vias) * s
    vias = [(x, y, d / 2) for x in xs for y in (-w / 2, w / 2)]
    margem = max(5 * sub.h_mm, 2.0)
    portas = [Porta(0.0, 0.0, w_porta, "xmin"), Porta(total, 0.0, w_porta, "xmax")]
    lay = Layout(nome, sub, topo, (0.0, -meia_placa - margem, total, meia_placa + margem), portas=portas,
                 info=info or {}, perfil=(np.array([0.0, total]), np.array([w_porta, w_porta])), vias=vias)
    lay.celula = (x1 - 0.3 * l_cunha, x1 + 4 * s)
    lay.cotas = [_cota_v(x1 + 1.5 * s, w, f"W = {w:.3f}"),
                 _cota_h(xs[0], xs[1] if n_vias > 1 else xs[0] + s, w / 2 + d, f"s = {s:.3f}"),
                 _cota_h(l_alim, x1, -w_cunha, f"cunha = {l_cunha:.2f}"),
                 _cota_v(x1 - 0.05, w_cunha, f"Wc = {w_cunha:.3f}")]
    return lay


def linha_convencional(tipo, elemento, d, sub: Substrato, l_alim=5.0, nome="linha", info=None):
    """Linhas planares convencionais (cap. 0 / cap. 1 do livro).

    tipo: "Microstrip", "CPW", "CPW com terra (GCPW)" ou "Coplanar strips (CPS)".
    elemento: "Trecho de linha", "Stub em aberto", "Stub em curto" ou "Transformador λ/4"
    (stubs e transformador só em microstrip). d: dimensões em mm (chaves w, g, comp, w_stub,
    comp_stub, w1, w2, wt, raio_via, passo_vias, largura_terra).
    """
    info = info or {}
    if tipo == "Microstrip":
        if elemento == "Transformador λ/4":
            x, w = perfil_degraus([(l_alim, d["w1"]), (d["comp"], d["wt"]), (l_alim, d["w2"])])
            lay = _montar(nome, sub, x, w, d["w1"], 0.0, info=info)
            lay.portas[-1] = Porta(lay.comprimento, 0.0, d["w2"], "xmax")
            lay.celula = (l_alim, l_alim + d["comp"])
            lay.cotas = [_cota_h(l_alim, l_alim + d["comp"], d["wt"], f"λ/4 = {d['comp']:.3f}"),
                         _cota_v(l_alim + d["comp"] / 2, d["wt"], f"Wt = {d['wt']:.3f}")]
            return lay
        if elemento.startswith("Stub"):
            w, ws, ls = d["w"], d["w_stub"], d["comp_stub"]
            comp = 2 * l_alim + ws
            xc = comp / 2
            topo = [_ret(0.0, -w / 2, comp, w / 2), _ret(xc - ws / 2, w / 2, xc + ws / 2, w / 2 + ls)]
            vias = []
            topo_y = w / 2 + ls
            if elemento == "Stub em curto":
                lado = max(2 * d["raio_via"] + 0.4, ws)
                topo.append(_ret(xc - lado / 2, topo_y, xc + lado / 2, topo_y + lado))
                vias.append((xc, topo_y + lado / 2, d["raio_via"]))
                topo_y += lado
            margem = max(5 * sub.h_mm, 2.0)
            portas = [Porta(0.0, 0.0, w, "xmin"), Porta(comp, 0.0, w, "xmax")]
            # sem perfil: o corte da pré-visualização passa pelo centro (stub)
            lay = Layout(nome, sub, topo, (0.0, -w / 2 - margem, comp, topo_y + margem), portas=portas,
                         info=info, vias=vias)
            lay.celula = (0.0, comp)
            lay.cotas = [(xc + ws, w / 2, xc + ws, w / 2 + ls, f"stub = {ls:.3f}"),
                         _cota_h(xc - ws / 2, xc + ws / 2, -w, f"Ws = {ws:.3f}"),
                         _cota_v(0.5, w, f"W = {w:.3f}")]
            return lay
        x = np.array([0.0, d["comp"]])
        ws = np.array([d["w"], d["w"]])
        lay = _montar(nome, sub, x, ws, d["w"], l_alim, info=info)
        lay.celula = (l_alim, l_alim + min(d["comp"], 6 * d["w"] + 5))
        lay.cotas = [_cota_v(l_alim + 1, d["w"], f"W = {d['w']:.3f}")]
        return lay
    if elemento != "Trecho de linha":
        raise ValueError("Stubs e transformador estão disponíveis só em microstrip.")
    s, g = d["w"], d["g"]
    comp = 2 * l_alim + d["comp"]
    wt = d["largura_terra"]
    if tipo == "Coplanar strips (CPS)":
        topo = [_ret(0.0, g / 2, comp, g / 2 + s), _ret(0.0, -g / 2 - s, comp, -g / 2)]
        borda = g / 2 + s
        portas = [Porta(0.0, 0.0, g + 2 * s, "xmin"), Porta(comp, 0.0, g + 2 * s, "xmax")]
        perfil = (np.array([0.0, comp]), np.array([g + 2 * s, g + 2 * s]))
        cotas = [_cota_v(l_alim + 1, g, f"S = {g:.3f}"), _cota_h(l_alim, l_alim + 2, g / 2 + s + 0.3,
                                                                f"W = {s:.3f}")]
        vias = []
    else:
        topo = [_ret(0.0, -s / 2, comp, s / 2), _ret(0.0, s / 2 + g, comp, s / 2 + g + wt),
                _ret(0.0, -s / 2 - g - wt, comp, -s / 2 - g)]
        borda = s / 2 + g + wt
        portas = [Porta(0.0, 0.0, s, "xmin"), Porta(comp, 0.0, s, "xmax")]
        perfil = (np.array([0.0, comp]), np.array([s, s]))
        cotas = [_cota_v(l_alim + 1, s, f"S = {s:.3f}"),
                 _cota_h(l_alim + 2, l_alim + 2 + g, s / 2 + g + 0.3, f"G = {g:.3f}")]
        vias = []
        if tipo == "CPW com terra (GCPW)" and d.get("passo_vias", 0) > 0:
            r = d["raio_via"]
            y = s / 2 + g + max(r + 0.2, 0.5)
            if y + r > borda:
                raise ValueError("Plano lateral estreito demais para a cerca de vias.")
            xs = np.arange(r + 0.2, comp - r - 0.2, d["passo_vias"])
            vias = [(xv, sy * y, r) for xv in xs for sy in (-1, 1)]
    margem = max(3 * sub.h_mm, 1.0)
    lay = Layout(nome, sub, topo, (0.0, -borda - margem, comp, borda + margem), portas=portas, info=info,
                 perfil=perfil, vias=vias, sem_terra=tipo in ("CPW", "Coplanar strips (CPS)"))
    lay.celula = (l_alim, l_alim + min(d["comp"], 8 * (s + 2 * g) + 3))
    lay.cotas = cotas
    return lay
