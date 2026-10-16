"""Linha CRLH em microstrip pela abordagem CL-loaded (Fig. 3.27, seção 3.5.1).

Célula T (Fig. 3.20): braços LR/2 + 2CL; derivação LL ‖ CR. Implementação:
  - LR, CR: trecho da linha hospedeira de impedância √(LR/CR) e período p;
  - CL: capacitor interdigital em série (fórmula de Bahl);
  - LL: stub em curto (via) em derivação no centro da célula.
O livro não traz síntese geométrica para esta linha (remete a Caloz e Itoh): as
fórmulas aqui são estimativas iniciais (stub curto, via de Goldfarb e Pucel,
capacitor de Bahl), sem parasitas do interdigital e do T de junção. Feche o
projeto com simulação EM e extração (modelo "CRLH (CL-loaded)" do Apêndice G).
"""
from dataclasses import dataclass

import numpy as np

from .crlh import CelulaCRLH
from .csrr import capacitancia_interdigital
from .microstrip import C0, Substrato, analise, sintese

MU0 = 4e-7 * np.pi


def indutancia_via(raio_mm, sub: Substrato):
    """Via cilíndrica até o terra (Goldfarb e Pucel, IEEE MGWL 1991) [H]."""
    h, r = sub.h_mm * 1e-3, raio_mm * 1e-3
    return MU0 / (2 * np.pi) * (h * np.log((h + np.hypot(r, h)) / r) + 1.5 * (r - np.hypot(r, h)))


@dataclass
class GeometriaCL:
    periodo: float       # mm
    w_linha: float       # mm
    n_dedos: int
    comp_dedo: float     # mm, dedos dos capacitores internos (CL)
    comp_dedo_ext: float  # mm, dedos dos capacitores das pontas (2CL)
    larg_dedo: float     # mm (largura = espaçamento)
    w_stub: float        # mm
    comp_stub: float     # mm
    raio_via: float      # mm
    l_via: float         # H, indutância estimada da via
    l_stub: float        # H

    @property
    def extensao_idc(self):
        return self.comp_dedo + 3 * self.larg_dedo


def projeto(c: CelulaCRLH, sub: Substrato, n_dedos=8, larg_dedo=0.1, z_stub=100.0, raio_via=0.2):
    """Geometria inicial a partir dos elementos LR, CR, LL, CL da célula."""
    zh = np.sqrt(c.LR / c.CR)
    w = sintese(zh, sub)
    _, ee = analise(w, sub)
    periodo = c.LR * C0 / (zh * np.sqrt(ee)) * 1e3            # LR = Zh·√εeff·p/c
    c_mm = capacitancia_interdigital(n_dedos, 1.0, larg_dedo, sub)   # Bahl: C ∝ comprimento
    comp = c.CL / c_mm
    ws = sintese(z_stub, sub)
    _, ee_s = analise(ws, sub)
    lv = indutancia_via(raio_via, sub)
    if c.LL <= lv:
        raise ValueError(f"LL = {c.LL * 1e9:.3g} nH é menor que a indutância da via ({lv * 1e9:.3g} nH): "
                         "aumente LL ou o raio da via.")
    ls = (c.LL - lv) * C0 / (z_stub * np.sqrt(ee_s)) * 1e3   # stub curto: L ≈ Zs·√εeff·ℓ/c
    geo = GeometriaCL(float(periodo), float(w), n_dedos, float(comp), float(2 * comp), larg_dedo, float(ws),
                      float(ls), raio_via, float(lv), float(c.LL - lv))
    avisos = []
    if geo.extensao_idc > 0.6 * periodo:
        avisos.append(f"Interdigital ({geo.extensao_idc:.2f} mm) ocupa mais de 60% do período ({periodo:.2f} mm): "
                      "aumente o nº de dedos ou reduza LR (o próprio capacitor tem indutância série).")
    if (2 * n_dedos - 1) * larg_dedo > 3 * w:
        avisos.append("Interdigital muito mais largo que a linha: reduza o nº de dedos.")
    lambda_mm = C0 / (c.wR / (2 * np.pi) * np.sqrt(ee_s)) * 1e3      # λ no stub em fR
    if ls > 0.15 * lambda_mm:
        avisos.append(f"Stub longo ({ls:.2f} mm > 0,15 λ em fR): a aproximação de stub curto perde validade.")
    return geo, avisos
