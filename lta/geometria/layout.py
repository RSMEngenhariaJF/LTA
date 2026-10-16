"""Modelo de layout planar (microstrip) independente de formato de saída.

Sistema de coordenadas: x ao longo da propagação (porta 1 em x=0), y
transversal centrado na trilha, unidades em mm. Camadas:
  - topo: polígonos de cobre (trilha)
  - terra: retângulo do plano de terra, com furos circulares e fendas
    poligonais opcionais (ex.: CSRR gravado no terra)
  - vias: furos metalizados ligando o cobre do topo ao terra
  - placa: contorno do substrato
"""
from dataclasses import dataclass, field

import numpy as np

from ..core.microstrip import Substrato


@dataclass
class Porta:
    x: float          # mm
    y: float          # mm
    largura: float    # mm (largura da trilha na porta)
    lado: str         # "xmin" ou "xmax"


@dataclass
class Layout:
    nome: str
    substrato: Substrato
    topo: list                      # lista de arrays (N,2) [mm], polígonos fechados (sem repetir o 1o ponto)
    placa: tuple                    # (x0, y0, x1, y1) [mm]
    furos_terra: list = field(default_factory=list)   # (cx, cy, r) [mm]
    portas: list = field(default_factory=list)
    info: dict = field(default_factory=dict)          # parâmetros de projeto (texto)
    # usados só na pré-visualização:
    cotas: list = field(default_factory=list)         # (x1, y1, x2, y2, texto) [mm]
    celula: tuple = None                              # (x0, x1) da 1a célula unitária
    perfil: tuple = None                              # (x, w) contorno da trilha
    fendas_terra: list = field(default_factory=list)  # polígonos (N,2) [mm] removidos do terra
    vias: list = field(default_factory=list)          # (cx, cy, r) [mm], furos metalizados topo-terra
    sem_terra: bool = False                           # CPW/CPS: sem plano de terra na face inferior

    def largura_em(self, x):
        """Largura da trilha na posição x (0 se fora da trilha)."""
        if self.perfil is None:
            return 0.0
        xs, ws = self.perfil
        if x < xs[0] or x > xs[-1]:
            return 0.0
        return float(np.interp(x, xs, ws))

    @property
    def terra(self):
        return self.placa

    @property
    def comprimento(self):
        return self.placa[2] - self.placa[0]

    @property
    def largura(self):
        return self.placa[3] - self.placa[1]


def intervalos_em_x(poligono, x):
    """Trechos [(y_a, y_b)] em que a vertical x atravessa o interior do polígono (par-ímpar)."""
    p = np.asarray(poligono, dtype=float)
    q = np.roll(p, -1, axis=0)
    ys = []
    for (x1, y1), (x2, y2) in zip(p, q):
        if (x1 <= x < x2) or (x2 <= x < x1):
            ys.append(y1 + (x - x1) * (y2 - y1) / (x2 - x1))
    ys.sort()
    return list(zip(ys[::2], ys[1::2]))


def poligono_de_perfil(x, w):
    """Contorno de uma trilha simétrica em y a partir do perfil de largura w(x)."""
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    sup = np.column_stack([x, w / 2])
    inf = np.column_stack([x[::-1], -w[::-1] / 2])
    pts = np.vstack([sup, inf])
    # remove pontos consecutivos repetidos
    d = np.any(np.abs(np.diff(pts, axis=0, append=pts[:1])) > 1e-9, axis=1)
    return pts[d]


def perfil_degraus(trechos):
    """Converte [(comprimento, largura), ...] em perfil (x, w) com degraus verticais."""
    xs, ws = [], []
    x = 0.0
    for comp, larg in trechos:
        xs += [x, x + comp]
        ws += [larg, larg]
        x += comp
    return np.array(xs), np.array(ws)


def simplificar_perfil(x, w, tol=1e-4):
    """Remove pontos colineares do perfil (reduz tamanho dos arquivos)."""
    keep = [0]
    for i in range(1, len(x) - 1):
        x0, w0 = x[keep[-1]], w[keep[-1]]
        if abs(x[i + 1] - x0) < 1e-12:
            keep.append(i)
            continue
        w_interp = w0 + (w[i + 1] - w0) * (x[i] - x0) / (x[i + 1] - x0)
        if abs(w_interp - w[i]) > tol or abs(x[i] - x[i - 1]) < 1e-12:
            keep.append(i)
    keep.append(len(x) - 1)
    return np.asarray(x)[keep], np.asarray(w)[keep]


def margem_lateral(sub: Substrato, w_max):
    """Meia-largura da placa: borda do cobre + max(5h, 2 mm) de dielétrico/terra."""
    return w_max / 2 + max(5 * sub.h_mm, 2.0)
