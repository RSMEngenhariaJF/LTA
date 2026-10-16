"""Pré-visualização da estrutura: vistas com cotas, corte e comparação com referência."""
import logging
import tkinter as tk
from tkinter import ttk

import numpy as np
from matplotlib.patches import Circle, Polygon, Rectangle

from ..geometria.layout import intervalos_em_x
from .comum import Grafico

COR_COBRE = "#d98c1f"
COR_BORDA = "#8a5208"
COR_SUBSTRATO = "#dfe8d8"
COR_FURO = "#1f4e9a"
COR_REF = "#7a7a7a"
COR_COTA = "#b0122b"

# Com aspecto igual e "datalim", o matplotlib amplia os limites para preencher o
# painel (desejado aqui) e registra um aviso informativo a cada desenho.
logging.getLogger("matplotlib.axes._base").setLevel(logging.ERROR)


def _desenhar_cota(ax, cota, cor=COR_COTA, etiqueta=None):
    """Seta de cota; com etiqueta, escreve só a letra (valor vai para o quadro de cotas)."""
    x1, y1, x2, y2, texto = cota
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), annotation_clip=False,
                arrowprops=dict(arrowstyle="<->", color=cor, lw=0.9, shrinkA=0, shrinkB=0))
    vertical = abs(x2 - x1) < abs(y2 - y1)
    ax.annotate(etiqueta or texto, ((x1 + x2) / 2, (y1 + y2) / 2), xytext=(3, 0) if vertical else (0, 3),
                textcoords="offset points", fontsize=8 if etiqueta else 7.5, color=cor,
                weight="bold" if etiqueta else "normal", annotation_clip=False,
                ha="left" if vertical else "center", va="center" if vertical else "bottom",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))


def _copper(ax, layout, ref=False):
    for p in layout.topo:
        if ref:
            ax.add_patch(Polygon(p, closed=True, fill=False, ec=COR_REF, ls="--", lw=1.0, zorder=5))
        else:
            ax.add_patch(Polygon(p, closed=True, fc=COR_COBRE, ec=COR_BORDA, lw=0.6, alpha=0.9, zorder=3))
    for cx, cy, r in layout.furos_terra:
        if ref:
            ax.add_patch(Circle((cx, cy), r, fill=False, ec=COR_REF, ls=":", lw=1.0, zorder=5))
        else:
            ax.add_patch(Circle((cx, cy), r, fc="white", ec=COR_FURO, ls="--", lw=1, zorder=2))
    for cx, cy, r in layout.vias:
        ax.add_patch(Circle((cx, cy), r, fc="none" if ref else "#333", ec=COR_REF if ref else "k",
                            lw=0.8, zorder=6))
    for p in layout.fendas_terra:
        if ref:
            ax.add_patch(Polygon(p, closed=True, fill=False, ec=COR_REF, ls=":", lw=1.0, zorder=5))
        else:
            ax.add_patch(Polygon(p, closed=True, fc="#dbe8f7", ec=COR_FURO, ls="--", lw=0.8, zorder=4, alpha=0.8))


def _mesmo_tipo(ref, atual):
    return ref is not None and ref.info.get("tipo") == atual.info.get("tipo")


def _diferencas(ref, atual):
    if ref is None:
        return []
    out = []
    for k in sorted(set(ref.info) | set(atual.info)):
        a, b = ref.info.get(k, "—"), atual.info.get(k, "—")
        if str(a) != str(b):
            out.append(f"{k}: {a} → {b}")
    for nome, a, b in (("comprimento_mm", ref.comprimento, atual.comprimento),
                       ("largura_placa_mm", ref.largura, atual.largura)):
        if abs(a - b) > 1e-6:
            out.append(f"{nome}: {a:.3f} → {b:.3f}")
    return out


class PreviaEstrutura(ttk.Frame):
    """Vista geral, zoom da célula com cotas e corte transversal.

    A referência (desenhada tracejada) é o estado anterior, ou um estado
    fixado pelo usuário com o botão 'Fixar referência'.
    """

    def __init__(self, mestre):
        super().__init__(mestre)
        barra = ttk.Frame(self)
        barra.pack(fill="x", padx=4, pady=(4, 0))
        self.ver_ref = tk.BooleanVar(value=True)
        ttk.Checkbutton(barra, text="Mostrar referência (tracejado)", variable=self.ver_ref,
                        command=self._redesenhar).pack(side="left")
        ttk.Button(barra, text="Fixar referência", command=self.fixar).pack(side="left", padx=6)
        ttk.Button(barra, text="Soltar referência", command=self.soltar).pack(side="left")
        self.rotulo_ref = ttk.Label(barra, text="Referência: estado anterior", foreground="#555")
        self.rotulo_ref.pack(side="left", padx=10)

        # o texto é empacotado antes do gráfico para não ser espremido por ele
        self.texto = tk.Text(self, height=4, wrap="word", font=("Consolas", 9), state="disabled")
        self.texto.pack(side="bottom", fill="x", padx=4, pady=(0, 4))
        self.grafico = Grafico(self)
        self.grafico.pack(fill="both", expand=True)

        self.atual = None
        self.anterior = None
        self.fixada = None
        self.erro = None
        self._redesenhar()

    # -- referência -------------------------------------------------------
    def fixar(self):
        if self.atual:
            self.fixada = self.atual
            self.rotulo_ref.configure(text=f"Referência fixada: {self.atual.nome}")
            self._redesenhar()

    def soltar(self):
        self.fixada = None
        self.rotulo_ref.configure(text="Referência: estado anterior")
        self._redesenhar()

    @property
    def referencia(self):
        return self.fixada if self.fixada is not None else self.anterior

    # -- entrada ----------------------------------------------------------
    def mostrar(self, layout):
        if self.atual is not None and (_diferencas(self.atual, layout) or self.atual.nome != layout.nome):
            self.anterior = self.atual
        self.atual = layout
        self.erro = None
        self._redesenhar()

    def mostrar_erro(self, msg):
        self.erro = msg
        self._redesenhar()

    def _escrever(self, txt):
        self.texto.configure(state="normal")
        self.texto.delete("1.0", "end")
        self.texto.insert("end", txt)
        self.texto.configure(state="disabled")

    # -- desenho ----------------------------------------------------------
    def _redesenhar(self):
        fig = self.grafico.fig
        fig.clear()
        lay = self.atual
        if lay is None:
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5, f"⚠ {self.erro}" if self.erro else "Edite os parâmetros para ver a estrutura",
                    ha="center", va="center", color="#b0122b" if self.erro else "k", wrap=True)
            ax.set_axis_off()
            self._escrever("")
            self.grafico.desenhar()
            return
        ref = self.referencia if self.ver_ref.get() else None
        if ref is lay or not _mesmo_tipo(ref, lay):
            ref = None   # não sobrepõe estruturas de tipos diferentes
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.35])
        ax_geral = fig.add_subplot(gs[0, :])
        ax_cel = fig.add_subplot(gs[1, 0])
        ax_corte = fig.add_subplot(gs[1, 1])
        self._vista_geral(ax_geral, lay, ref)
        x_corte = self._vista_celula(ax_cel, lay, ref)
        self._corte(ax_corte, lay, x_corte)
        if self.erro:
            fig.text(0.5, 0.5, f"⚠ {self.erro}\n(mostrando a última estrutura válida)", ha="center",
                     va="center", color="#b0122b", fontsize=11,
                     bbox=dict(fc="#fff3f3", ec="#b0122b", boxstyle="round"))
        self.grafico.desenhar()

        s = lay.substrato
        linhas = [f"{lay.nome} | {lay.comprimento:.3f} x {lay.largura:.3f} mm | "
                  f"{s.nome}: εr={s.er:g}, h={s.h_mm:g} mm, t={s.t_um:g} µm"]
        r = self.referencia
        if r is not None and not _mesmo_tipo(r, lay):
            linhas.append(f"Referência é de outro tipo ({r.info.get('tipo', r.nome)}): sem comparação.")
        elif r is not None:
            dif = _diferencas(r, lay)
            linhas.append("Alterações em relação à referência: " + " | ".join(dif) if dif
                          else "Sem alterações em relação à referência.")
        if lay.perfil is not None and np.min(lay.perfil[1]) < 0.1:
            linhas.append(f"⚠ Largura mínima de trilha {np.min(lay.perfil[1]):.3f} mm: abaixo de ~0,1 mm "
                          "costuma exceder a capacidade de fabricação. Verificar com o fornecedor.")
        self._escrever("\n".join(linhas))

    def _vista_geral(self, ax, lay, ref):
        x0, y0, x1, y1 = lay.placa
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=COR_SUBSTRATO, ec="#556b2f", lw=1))
        _copper(ax, lay)
        if ref is not None:
            rx0, ry0, rx1, ry1 = ref.placa
            ax.add_patch(Rectangle((rx0, ry0), rx1 - rx0, ry1 - ry0, fill=False, ec=COR_REF, ls="--", lw=0.8))
            _copper(ax, ref, ref=True)
        if lay.celula:
            c0, c1 = lay.celula
            ax.add_patch(Rectangle((c0, y0), c1 - c0, y1 - y0, fill=False, ec="#1f77b4", ls="-.", lw=1))
            ax.annotate("detalhe (zoom)", (c0, y1), xytext=(0, 2), textcoords="offset points",
                        fontsize=7, color="#1f77b4", va="bottom")
        for i, p in enumerate(lay.portas, start=1):
            ax.annotate(f"P{i}", (p.x, p.y), xytext=(-14 if p.lado == "xmin" else 4, 4),
                        textcoords="offset points", fontsize=8, weight="bold")
        _desenhar_cota(ax, (x0, y0 - 0.06 * (y1 - y0), x1, y0 - 0.06 * (y1 - y0),
                            f"{lay.comprimento:.2f} mm"), cor="#333")
        _desenhar_cota(ax, (x1 + 0.01 * (x1 - x0), y0, x1 + 0.01 * (x1 - x0), y1,
                            f"{lay.largura:.2f}"), cor="#333")
        xr = (min(x0, ref.placa[0]) if ref else x0, max(x1, ref.placa[2]) if ref else x1)
        yr = (min(y0, ref.placa[1]) if ref else y0, max(y1, ref.placa[3]) if ref else y1)
        ax.set_xlim(xr[0] - 0.02 * (xr[1] - xr[0]), xr[1] + 0.09 * (xr[1] - xr[0]))
        ax.set_ylim(yr[0] - 0.25 * (yr[1] - yr[0]), yr[1] + 0.12 * (yr[1] - yr[0]))
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_title("Vista superior (cobre do topo; furos e fendas do terra tracejados)", fontsize=9)
        ax.set_xlabel("x [mm]", fontsize=8)
        ax.tick_params(labelsize=7)

    def _vista_celula(self, ax, lay, ref):
        if lay.celula:
            c0, c1 = lay.celula
        else:
            c0, c1 = lay.placa[0], lay.placa[2]
        m = 0.2 * (c1 - c0)
        xa, xb = max(lay.placa[0], c0 - m), min(lay.placa[2], c1 + m)
        x0, y0, x1, y1 = lay.placa
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=COR_SUBSTRATO, ec="none"))
        _copper(ax, lay)
        if ref is not None:
            _copper(ax, ref, ref=True)
        letras = "abcdefghij"
        quadro = []
        for letra, cota in zip(letras, lay.cotas):
            _desenhar_cota(ax, cota, etiqueta=letra)
            quadro.append(f"{letra}: {cota[4]}")
        if quadro:
            ax.text(0.01, 0.99, "\n".join(quadro), transform=ax.transAxes, va="top", ha="left",
                    fontsize=7.5, family="monospace", color=COR_COTA,
                    bbox=dict(fc="white", ec=COR_COTA, alpha=0.85, boxstyle="round,pad=0.3"))
        # extensão: cobre, furos e cotas da célula
        xa = min([xa] + [min(c[0], c[2]) for c in lay.cotas])
        xb = max([xb] + [max(c[0], c[2]) for c in lay.cotas])
        ys = [0.0]
        if lay.perfil is not None:
            xs, ws = lay.perfil
            sel = (xs >= xa) & (xs <= xb)
            ys += list(ws[sel] / 2) + list(-ws[sel] / 2)
        ys += [cy + s * r for cx, cy, r in lay.furos_terra if xa - r <= cx <= xb + r for s in (-1, 1)]
        for p in lay.fendas_terra:
            p = np.asarray(p)
            if p[:, 0].max() >= xa and p[:, 0].min() <= xb:
                ys += [p[:, 1].min(), p[:, 1].max()]
        ys += [c[1] for c in lay.cotas] + [c[3] for c in lay.cotas]
        ya, yb = min(ys), max(ys)
        pad = 0.25 * max(yb - ya, 1e-3)
        mx = 0.05 * (xb - xa)
        ax.set_xlim(xa - mx, xb + mx)
        ax.set_ylim(ya - pad, yb + 2.2 * pad)   # espaço extra no topo para o quadro de cotas
        ax.set_aspect("equal", adjustable="datalim")
        # posição do corte: furo da célula, ou onde a largura mais difere da porta
        x_corte = (c0 + c1) / 2
        furos = [cx for cx, _, _ in lay.furos_terra if c0 <= cx <= c1]
        fendas = [float(np.mean(np.asarray(p)[:, 0])) for p in lay.fendas_terra
                  if c0 <= np.mean(np.asarray(p)[:, 0]) <= c1]
        if furos:
            x_corte = furos[0]
        elif fendas:
            x_corte = fendas[0]
        elif lay.perfil is not None:
            xx = np.linspace(c0, c1, 400)
            wp = lay.portas[0].largura
            x_corte = float(xx[np.argmax([abs(lay.largura_em(x) - wp) for x in xx])])
        ax.axvline(x_corte, color="#1f77b4", ls="-.", lw=0.8)
        ax.annotate("A", (x_corte, yb + pad * 0.8), color="#1f77b4", fontsize=8, ha="center")
        detalhe = "personalizado" in str(lay.info.get("tipo", ""))
        ax.set_title("Trecho de maior largura, com cotas [mm]" if detalhe else "Célula unitária com cotas [mm]",
                     fontsize=9)
        ax.set_xlabel("x [mm]", fontsize=8)
        ax.tick_params(labelsize=7)
        return x_corte

    def _corte(self, ax, lay, x):
        s = lay.substrato
        _, y0, _, y1 = lay.placa
        h, t = s.h_mm, s.t_um / 1000
        t_vis = max(t, 0.04 * h)  # cobre exagerado para ficar visível
        w = lay.largura_em(x)
        ax.add_patch(Rectangle((y0, 0), y1 - y0, h, fc=COR_SUBSTRATO, ec="#556b2f", lw=0.8))
        # terra com eventual abertura do furo neste x
        aberturas = [(cy - np.sqrt(r**2 - (x - cx) ** 2), cy + np.sqrt(r**2 - (x - cx) ** 2))
                     for cx, cy, r in lay.furos_terra if abs(x - cx) < r]
        aberturas += [iv for p in lay.fendas_terra for iv in intervalos_em_x(p, x)]
        limites = [y0] + [v for a in sorted(aberturas) for v in a] + [y1]
        for a, b in zip(limites[::2], limites[1::2]):
            ax.add_patch(Rectangle((a, -t_vis), b - a, t_vis, fc=COR_COBRE, ec=COR_BORDA, lw=0.5))
        if w > 0:
            ax.add_patch(Rectangle((-w / 2, h), w, t_vis, fc=COR_COBRE, ec=COR_BORDA, lw=0.5))
            _desenhar_cota(ax, (-w / 2, h + t_vis + 0.25 * h, w / 2, h + t_vis + 0.25 * h, f"W = {w:.3f}"))
        _desenhar_cota(ax, (y1 * 0.92, 0, y1 * 0.92, h, f"h = {h:g}"))
        ax.text(y0 + 0.03 * (y1 - y0), h / 2, f"εr = {s.er:g}", fontsize=8, va="center")
        for a, b in aberturas:
            ax.annotate("fenda" if lay.fendas_terra else "furo", ((a + b) / 2, -t_vis), xytext=(0, -10),
                        textcoords="offset points", fontsize=7, ha="center", color=COR_FURO)
        ax.set_xlim(y0, y1)
        ax.set_ylim(-0.6 * h, 1.9 * h)
        ax.set_title(f"Corte A (x = {x:.2f} mm), cobre fora de escala", fontsize=9)
        ax.set_xlabel("y [mm]", fontsize=8)
        ax.set_ylabel("z [mm]", fontsize=8)
        ax.tick_params(labelsize=7)


# ---------------------------------------------------------------------------
# Esquemático da célula concentrada (aba Célula ABCD)
# ---------------------------------------------------------------------------
def _segmento(ax, p0, p1, tipo, rotulo, lado=1):
    """Desenha indutor ('L') ou capacitor ('C') entre p0 e p1 (com terminais retos)."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    comp = np.hypot(*d)
    u = d / comp
    v = np.array([-u[1], u[0]])
    a, b = 0.3 * comp, 0.7 * comp
    pa, pb = p0 + u * a, p0 + u * b
    ax.plot(*zip(p0, pa), color="k", lw=1.2)
    ax.plot(*zip(pb, p1), color="k", lw=1.2)
    if tipo == "L":
        s = np.linspace(0, 1, 200)
        n = 4
        fase = s * n * np.pi
        ponto = pa[None, :] + (b - a) * s[:, None] * u[None, :] + 0.12 * np.abs(np.sin(fase))[:, None] * v[None, :] * lado
        ax.plot(ponto[:, 0], ponto[:, 1], color="k", lw=1.2)
    else:
        meio = (pa + pb) / 2
        g = 0.05
        ax.plot(*zip(pa, meio - u * g), color="k", lw=1.2)
        ax.plot(*zip(meio + u * g, pb), color="k", lw=1.2)
        for sgn in (-1, 1):
            c = meio + sgn * u * g
            ax.plot(*zip(c - 0.16 * v, c + 0.16 * v), color="k", lw=2)
    horizontal = abs(u[0]) > abs(u[1])
    m = (p0 + p1) / 2 + v * (0.32 if horizontal else 0.24) * lado
    if horizontal:
        ha = "center"
    else:  # ramo vertical: texto para fora do símbolo
        ha = "left" if (v * lado)[0] > 0 else "right"
    ax.text(*m, rotulo, fontsize=8, ha=ha, va="center", color="#1f4e9a")


def _fmt_l(val):
    return f"{val * 1e9:.3g} nH"


def _fmt_c(val):
    return f"{val * 1e12:.3g} pF"


def _caixa(ax, p0, p1, rotulo, lado=1):
    """Impedância genérica (retângulo) entre p0 e p1, para ramos com L e C juntos."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    u = (p1 - p0) / np.hypot(*(p1 - p0))
    v = np.array([-u[1], u[0]])
    pa, pb = p0 + (p1 - p0) * 0.3, p0 + (p1 - p0) * 0.7
    ax.plot(*zip(p0, pa), color="k", lw=1.2)
    ax.plot(*zip(pb, p1), color="k", lw=1.2)
    cantos = [pa + 0.12 * v, pb + 0.12 * v, pb - 0.12 * v, pa - 0.12 * v, pa + 0.12 * v]
    ax.plot(*zip(*cantos), color="k", lw=1.2)
    m = (p0 + p1) / 2 + v * 0.32 * lado
    ax.text(*m, rotulo, fontsize=8, ha="center", va="center", color="#1f4e9a")


def _ramo(ax, p0, p1, elems, lado=1):
    """Desenha os elementos de um ramo em série ao longo de p0 -> p1 (ou caixa se forem paralelos)."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    if not elems:
        ax.plot(*zip(p0, p1), color="k", lw=1.2)
        return
    passos = np.linspace(0, 1, len(elems) + 1)
    for (tipo, rot), a, b in zip(elems, passos[:-1], passos[1:]):
        _segmento(ax, p0 + (p1 - p0) * a, p0 + (p1 - p0) * b, tipo, rot, lado)


def _elementos_serie(ls, cs, fator=1.0, nome_l="Ls", nome_c="Cs"):
    out = []
    if ls > 0:
        out.append(("L", f"{nome_l} = {_fmt_l(ls * fator)}"))
    if cs > 0:
        out.append(("C", f"{nome_c} = {_fmt_c(cs / fator)}"))
    return out


def _desenhar_trelica(ax, ls, cs, lp, cp, n):
    x0, x1, topo, base = 0.6, 4.2, 2.0, 0.0
    for y in (topo, base):
        ax.plot([0, x0], [y, y], color="k", lw=1.2)
        ax.plot([x1, x1 + 0.6], [y, y], color="k", lw=1.2)
        _ramo(ax, (x0, y), (x1, y), _elementos_serie(ls, cs), lado=1 if y == topo else -1)
    if lp > 0 and cp > 0:
        rot = f"Cp ‖ Lp ({_fmt_c(cp)}, {_fmt_l(lp)})"
        diag = None
    else:
        diag = ("C", f"Cp = {_fmt_c(cp)}") if cp > 0 else ("L", f"Lp = {_fmt_l(lp)}")
    for (a, b), (t0, t1), lado in ((((x0, topo), (x1, base)), (0.12, 0.42), 1),
                                   (((x0, base), (x1, topo)), (0.58, 0.88), -1)):
        a, b = np.array(a), np.array(b)
        pa, pb = a + (b - a) * t0, a + (b - a) * t1
        ax.plot(*zip(a, pa), color="k", lw=1.2)
        ax.plot(*zip(pb, b), color="k", lw=1.2)
        if diag is None:
            _caixa(ax, pa, pb, rot, lado)
        else:
            _segmento(ax, pa, pb, diag[0], diag[1], lado)
    for x in (x0, x1):
        for y in (topo, base):
            ax.plot(x, y, "ko", ms=3)
    return x1 + 0.6, topo, base, "Treliça (X) – rede balanceada"


def _desenhar_ponte(ax, ls, cs, lp, cp, n):
    topo, base, ponte = 1.6, 0.0, 2.7
    xa, xm, xb, xf = 0.5, 2.3, 4.1, 4.6
    ax.plot([0, xa], [topo, topo], color="k", lw=1.2)
    ax.plot([xb, xf], [topo, topo], color="k", lw=1.2)
    _ramo(ax, (xa, topo), (xm, topo), _elementos_serie(ls, cs, 1, "Z1: Ls", "Z1: Cs"))
    _ramo(ax, (xm, topo), (xb, topo), _elementos_serie(ls, cs, 1, "Z1: Ls", "Z1: Cs"))
    # derivação Z2 e ponte Z3 = 4*Z2
    if lp > 0 and cp > 0:
        _caixa(ax, (xm, topo), (xm, base), f"Z2: Cp ‖ Lp", lado=1)
        z3 = None
    elif cp > 0:
        _segmento(ax, (xm, topo), (xm, base), "C", f"Z2: Cp = {_fmt_c(cp)}", lado=1)
        z3 = ("C", f"Z3 = 4·Z2: Cp/4 = {_fmt_c(cp / 4)}")
    else:
        _segmento(ax, (xm, topo), (xm, base), "L", f"Z2: Lp = {_fmt_l(lp)}", lado=1)
        z3 = ("L", f"Z3 = 4·Z2: 4Lp = {_fmt_l(4 * lp)}")
    ax.plot([xa, xa], [topo, ponte], color="k", lw=1.2)
    ax.plot([xb, xb], [topo, ponte], color="k", lw=1.2)
    if z3 is None:
        _caixa(ax, (xa, ponte), (xb, ponte), "Z3 = 4·Z2")
    else:
        _segmento(ax, (xa, ponte), (xb, ponte), z3[0], z3[1])
    for x, y in ((xa, topo), (xm, topo), (xb, topo), (xm, base)):
        ax.plot(x, y, "ko", ms=3)
    return xf, topo, base, "T em ponte (Z3 = 4·Z2)"


def desenhar_celula(ax, topologia, ls=0.0, cs=0.0, lp=0.0, cp=0.0, n=1):
    """Esquemático com símbolos de L e C e valores de cada elemento."""
    ax.clear()
    ax.set_axis_off()
    if topologia in ("Treliça (X)", "T em ponte"):
        func = _desenhar_trelica if topologia == "Treliça (X)" else _desenhar_ponte
        x, topo, base, titulo = func(ax, ls, cs, lp, cp, n)
        if topologia == "T em ponte":
            ax.plot([0, x], [base, base], color="k", lw=1.2)
        for xp, nome in ((0, "porta 1"), (x, "porta 2")):
            ax.plot(xp, topo, "o", mfc="white", mec="k", ms=6)
            ax.plot(xp, base, "o", mfc="white", mec="k", ms=6)
            ax.text(xp, topo + 0.25, nome, ha="center", fontsize=8)
        ax.text(x / 2, -0.6, f"{titulo}  ×  {n} em cascata", ha="center", fontsize=10, weight="bold")
        ax.set_xlim(-0.6, x + 0.6)
        ax.set_ylim(-1.0, 3.5)
        ax.set_aspect("equal")
        return
    ops = {"T": [("S", 0.5), ("P", 1.0), ("S", 0.5)], "Pi": [("P", 0.5), ("S", 1.0), ("P", 0.5)],
           "L": [("S", 1.0), ("P", 1.0)], "L invertida": [("P", 1.0), ("S", 1.0)]}[topologia]
    topo, base = 1.6, 0.0
    x = 0.0
    ax.plot([x, x + 0.4], [topo, topo], color="k", lw=1.2)
    x += 0.4
    for tipo, fator in ops:
        meio = fator != 1
        if tipo == "S":
            elems = [("L", f"{'Ls/2' if meio else 'Ls'} = {_fmt_l(ls * fator)}")] if ls > 0 else []
            if cs > 0:
                elems.append(("C", f"{'2Cs' if meio else 'Cs'} = {_fmt_c(cs / fator)}"))
            if not elems:
                ax.plot([x, x + 1.6], [topo, topo], color="k", lw=1.2)
                x += 1.6
            for e, rot in elems:
                _segmento(ax, (x, topo), (x + 1.6, topo), e, rot)
                x += 1.6
        else:
            ax.plot(x, topo, "ko", ms=3)
            elems = [("C", f"{'Cp/2' if meio else 'Cp'} = {_fmt_c(cp * fator)}")] if cp > 0 else []
            if lp > 0:
                elems.append(("L", f"{'2Lp' if meio else 'Lp'} = {_fmt_l(lp / fator)}"))
            if len(elems) == 2:
                for dx, (e, rot) in zip((-0.35, 0.35), elems):
                    ax.plot([x - 0.35, x + 0.35], [topo - 0.2, topo - 0.2], color="k", lw=1.2)
                    ax.plot([x - 0.35, x + 0.35], [base + 0.2, base + 0.2], color="k", lw=1.2)
                    _segmento(ax, (x + dx, topo - 0.2), (x + dx, base + 0.2), e, rot, lado=-1 if dx < 0 else 1)
                ax.plot([x, x], [topo, topo - 0.2], color="k", lw=1.2)
                ax.plot([x, x], [base, base + 0.2], color="k", lw=1.2)
            elif elems:
                _segmento(ax, (x, topo), (x, base), elems[0][0], elems[0][1], lado=1)
            ax.plot(x, base, "ko", ms=3)
            ax.plot([x, x + 0.6], [topo, topo], color="k", lw=1.2)
            x += 0.6
    ax.plot([0, x], [base, base], color="k", lw=1.2)
    for xp, nome in ((0, "porta 1"), (x, "porta 2")):
        ax.plot([xp, xp], [base, topo], color="#999", lw=0.6, ls=":")
        ax.plot(xp, topo, "o", mfc="white", mec="k", ms=6)
        ax.plot(xp, base, "o", mfc="white", mec="k", ms=6)
        ax.text(xp, topo + 0.25, nome, ha="center", fontsize=8)
    ax.text(x / 2, -0.55, f"Célula {topologia}  ×  {n} em cascata", ha="center", fontsize=10, weight="bold")
    ax.set_xlim(-0.6, x + 0.6)
    ax.set_ylim(-0.9, topo + 0.6)
    ax.set_aspect("equal")


class PreviaCelula(ttk.Frame):
    def __init__(self, mestre):
        super().__init__(mestre)
        self.grafico = Grafico(self)
        self.grafico.pack(fill="both", expand=True)

    def mostrar(self, topologia, n, **el):
        ax = self.grafico.eixos(1, 1)[0]
        desenhar_celula(ax, topologia, n=n, **el)
        self.grafico.desenhar()

    def mostrar_erro(self, msg):
        ax = self.grafico.eixos(1, 1)[0]
        ax.set_axis_off()
        ax.text(0.5, 0.5, f"⚠ {msg}", ha="center", va="center", color="#b0122b")
        self.grafico.desenhar()
