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
