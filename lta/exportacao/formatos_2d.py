"""Exportação 2D do layout: DXF (R12), Gerber RS-274X, KiCad (.kicad_pcb) e GDSII.

Todos em mm. Uso típico:
  - Altium: importar DXF (File > Import > DXF/DWG) ou Gerber (CAMtastic),
            ou abrir o .kicad_pcb pelo Import Wizard (versões recentes).
  - KiCad:  abrir o .kicad_pcb diretamente (pressionar B para preencher o terra).
  - ADS:    importar GDSII ou DXF no Layout.
  - HFSS/CST: também aceitam DXF/GDSII, mas prefira os scripts 3D gerados.
"""
import datetime
import struct
import uuid

import numpy as np


def _circulo(cx, cy, r, n=64):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([cx + r * np.cos(t), cy + r * np.sin(t)])


def _retangulo(x0, y0, x1, y1):
    return np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])


# ---------------------------------------------------------------------------
# DXF R12
# ---------------------------------------------------------------------------
CAMADAS_DXF = {"TOPO_COBRE": 1, "TERRA": 5, "TERRA_FUROS": 3, "PLACA": 7}


def dxf(layout, caminho):
    linhas = ["0", "SECTION", "2", "HEADER", "9", "$ACADVER", "1", "AC1009", "0", "ENDSEC",
              "0", "SECTION", "2", "TABLES", "0", "TABLE", "2", "LAYER", "70", str(len(CAMADAS_DXF))]
    for nome, cor in CAMADAS_DXF.items():
        linhas += ["0", "LAYER", "2", nome, "70", "0", "62", str(cor), "6", "CONTINUOUS"]
    linhas += ["0", "ENDTAB", "0", "ENDSEC", "0", "SECTION", "2", "ENTITIES"]

    def poli(pts, camada):
        out = ["0", "POLYLINE", "8", camada, "66", "1", "10", "0.0", "20", "0.0", "30", "0.0", "70", "1"]
        for x, y in pts:
            out += ["0", "VERTEX", "8", camada, "10", f"{x:.6f}", "20", f"{y:.6f}", "30", "0.0"]
        return out + ["0", "SEQEND", "8", camada]

    for p in layout.topo:
        linhas += poli(p, "TOPO_COBRE")
    linhas += poli(_retangulo(*layout.terra), "TERRA")
    linhas += poli(_retangulo(*layout.placa), "PLACA")
    for p in layout.fendas_terra:
        linhas += poli(p, "TERRA_FUROS")
    for cx, cy, r in layout.furos_terra:
        linhas += ["0", "CIRCLE", "8", "TERRA_FUROS", "10", f"{cx:.6f}", "20", f"{cy:.6f}",
                   "30", "0.0", "40", f"{r:.6f}"]
    linhas += ["0", "ENDSEC", "0", "EOF"]
    caminho.write_text("\n".join(linhas) + "\n", encoding="ascii")
    return [caminho]


# ---------------------------------------------------------------------------
# Gerber RS-274X (com atributos X2 básicos)
# ---------------------------------------------------------------------------
def _gcoord(v):
    return f"{int(round(v * 1e6))}"


def _gerber_cab(funcao, polaridade="Positive"):
    return ["G04 Gerado por LTA - linhas periodicas*", "%FSLAX46Y46*%", "%MOMM*%",
            f"%TF.FileFunction,{funcao}*%", f"%TF.FilePolarity,{polaridade}*%", "%LPD*%", "G01*"]


def _gerber_regiao(pts):
    out = ["G36*", f"X{_gcoord(pts[0][0])}Y{_gcoord(pts[0][1])}D02*"]
    for x, y in pts[1:]:
        out.append(f"X{_gcoord(x)}Y{_gcoord(y)}D01*")
    out += [f"X{_gcoord(pts[0][0])}Y{_gcoord(pts[0][1])}D01*", "G37*"]
    return out


def gerber(layout, caminho_base):
    """Gera <base>_topo.gtl, <base>_terra.gbl e <base>_contorno.gko."""
    base = caminho_base.with_suffix("")
    arqs = []

    topo = _gerber_cab("Copper,L1,Top")
    for p in layout.topo:
        topo += _gerber_regiao(p)
    topo.append("M02*")
    arqs.append(base.parent / f"{base.name}_topo.gtl")
    arqs[-1].write_text("\n".join(topo) + "\n", encoding="ascii")

    terra = _gerber_cab("Copper,L2,Bot")
    diams = sorted({round(2 * r, 6) for _, _, r in layout.furos_terra})
    ap = {d: 10 + i for i, d in enumerate(diams)}
    terra += [f"%ADD{ap[d]}C,{d:.6f}*%" for d in diams]
    terra += _gerber_regiao(_retangulo(*layout.terra))
    if layout.furos_terra or layout.fendas_terra:
        # polaridade "clear": círculos e fendas removem cobre do plano de terra
        terra.append("%LPC*%")
        for cx, cy, r in layout.furos_terra:
            terra += [f"D{ap[round(2 * r, 6)]}*", f"X{_gcoord(cx)}Y{_gcoord(cy)}D03*"]
        for p in layout.fendas_terra:
            terra += _gerber_regiao(p)
    terra.append("M02*")
    arqs.append(base.parent / f"{base.name}_terra.gbl")
    arqs[-1].write_text("\n".join(terra) + "\n", encoding="ascii")

    cont = _gerber_cab("Profile,NP")
    cont += ["%ADD10C,0.100000*%", "D10*"]
    pts = _retangulo(*layout.placa)
    cont.append(f"X{_gcoord(pts[0][0])}Y{_gcoord(pts[0][1])}D02*")
    for x, y in list(pts[1:]) + [pts[0]]:
        cont.append(f"X{_gcoord(x)}Y{_gcoord(y)}D01*")
    cont.append("M02*")
    arqs.append(base.parent / f"{base.name}_contorno.gko")
    arqs[-1].write_text("\n".join(cont) + "\n", encoding="ascii")
    return arqs


# ---------------------------------------------------------------------------
# KiCad 6+ (.kicad_pcb). KiCad 7/8 abrem e convertem automaticamente.
# ---------------------------------------------------------------------------
_CAMADAS_KICAD = [(0, "F.Cu", "signal"), (31, "B.Cu", "signal"), (32, "B.Adhes", "user"),
                  (33, "F.Adhes", "user"), (34, "B.Paste", "user"), (35, "F.Paste", "user"),
                  (36, "B.SilkS", "user"), (37, "F.SilkS", "user"), (38, "B.Mask", "user"),
                  (39, "F.Mask", "user"), (40, "Dwgs.User", "user"), (41, "Cmts.User", "user"),
                  (42, "Eco1.User", "user"), (43, "Eco2.User", "user"), (44, "Edge.Cuts", "user"),
                  (45, "Margin", "user"), (46, "B.CrtYd", "user"), (47, "F.CrtYd", "user"),
                  (48, "B.Fab", "user"), (49, "F.Fab", "user")]
_ORIGEM_KICAD = (50.0, 100.0)  # desloca o desenho para dentro da folha


def _kxy(x, y):
    # KiCad: eixo y para baixo
    return f"(xy {x + _ORIGEM_KICAD[0]:.6f} {-y + _ORIGEM_KICAD[1]:.6f})"


def _kpts(pts):
    return "(pts " + " ".join(_kxy(x, y) for x, y in pts) + ")"


def kicad(layout, caminho):
    sub = layout.substrato
    t_mm = sub.t_um / 1000
    s = ["(kicad_pcb (version 20211014) (generator lta_linhas_periodicas)",
         f"  (general (thickness {sub.h_mm + 2 * t_mm:.4f}))",
         '  (paper "A4")',
         "  (layers"]
    s += [f'    ({n} "{nome}" {tipo})' for n, nome, tipo in _CAMADAS_KICAD]
    s += ["  )",
          "  (setup",
          "    (stackup",
          f'      (layer "F.Cu" (type "copper") (thickness {t_mm:.4f}))',
          f'      (layer "dielectric 1" (type "core") (thickness {sub.h_mm:.4f}) (material "{sub.nome}")'
          f" (epsilon_r {sub.er:g}) (loss_tangent {sub.tand:g}))",
          f'      (layer "B.Cu" (type "copper") (thickness {t_mm:.4f}))',
          '      (copper_finish "None")',
          "      (dielectric_constraints no)",
          "    )",
          "    (pad_to_mask_clearance 0)",
          "  )",
          '  (net 0 "")']
    for p in layout.topo:
        s.append(f'  (gr_poly {_kpts(p)} (layer "F.Cu") (width 0) (fill solid) (tstamp {uuid.uuid4()}))')
    x0, y0, x1, y1 = layout.placa
    s.append(f'  (gr_rect (start {x0 + _ORIGEM_KICAD[0]:.6f} {-y1 + _ORIGEM_KICAD[1]:.6f}) '
             f'(end {x1 + _ORIGEM_KICAD[0]:.6f} {-y0 + _ORIGEM_KICAD[1]:.6f}) '
             f'(layer "Edge.Cuts") (width 0.1) (fill none) (tstamp {uuid.uuid4()}))')
    # Plano de terra como zona (preenchida pelo KiCad: tecla B)
    s.append(f'  (zone (net 0) (net_name "") (layer "B.Cu") (tstamp {uuid.uuid4()}) (hatch edge 0.508)'
             " (connect_pads yes (clearance 0)) (min_thickness 0.1) (filled_areas_thickness no)"
             " (fill yes (thermal_gap 0.508) (thermal_bridge_width 0.508))"
             f" (polygon {_kpts(_retangulo(*layout.terra))}))")
    # Furos no terra: áreas de exclusão de preenchimento
    for pts in [_circulo(cx, cy, r) for cx, cy, r in layout.furos_terra] + list(layout.fendas_terra):
        s.append(f'  (zone (net 0) (net_name "") (layer "B.Cu") (tstamp {uuid.uuid4()}) (hatch edge 0.508)'
                 " (connect_pads (clearance 0)) (min_thickness 0.1) (filled_areas_thickness no)"
                 " (keepout (tracks allowed) (vias allowed) (pads allowed) (copperpour not_allowed)"
                 " (footprints allowed))"
                 " (fill (thermal_gap 0.508) (thermal_bridge_width 0.508))"
                 f" (polygon {_kpts(pts)}))")
    s.append(")")
    caminho.write_text("\n".join(s) + "\n", encoding="utf-8")
    return [caminho]


# ---------------------------------------------------------------------------
# GDSII (binário). Unidade de banco de dados = 1 nm; unidade de usuário = 1 um.
# ---------------------------------------------------------------------------
CAMADAS_GDS = {"TOPO_COBRE": 1, "TERRA": 2, "TERRA_FUROS": 3, "PLACA": 10}


def _real8(v):
    if v == 0:
        return bytes(8)
    sinal = 0x80 if v < 0 else 0
    v = abs(v)
    e = 64
    while v >= 1:
        v /= 16
        e += 1
    while v < 1 / 16:
        v *= 16
        e -= 1
    m = int(round(v * 2**56))
    if m >= 2**56:
        m //= 16
        e += 1
    return bytes([sinal | e]) + m.to_bytes(7, "big")


def _rec(tipo, dados=b""):
    return struct.pack(">HH", 4 + len(dados), tipo) + dados


def _str(s):
    b = s.encode("ascii")
    return b + (b"\0" if len(b) % 2 else b"")


def gds(layout, caminho, nome_celula=None):
    agora = datetime.datetime.now()
    data = struct.pack(">6h", agora.year, agora.month, agora.day, agora.hour, agora.minute, agora.second) * 2
    out = bytearray()
    out += _rec(0x0002, struct.pack(">h", 600))
    out += _rec(0x0102, data)
    out += _rec(0x0206, _str("LTA"))
    out += _rec(0x0305, _real8(1e-3) + _real8(1e-9))
    out += _rec(0x0502, data)
    out += _rec(0x0606, _str((nome_celula or layout.nome).upper()[:32]))

    def boundary(pts, camada):
        pts = np.vstack([pts, pts[:1]])
        if len(pts) > 8191:
            raise ValueError("Polígono com mais de 8190 vértices não é suportado em GDSII.")
        coords = np.round(pts * 1e6).astype(">i4").ravel().tobytes()
        return (_rec(0x0800) + _rec(0x0D02, struct.pack(">h", camada))
                + _rec(0x0E02, struct.pack(">h", 0)) + _rec(0x1003, coords) + _rec(0x1100))

    for p in layout.topo:
        out += boundary(p, CAMADAS_GDS["TOPO_COBRE"])
    out += boundary(_retangulo(*layout.terra), CAMADAS_GDS["TERRA"])
    out += boundary(_retangulo(*layout.placa), CAMADAS_GDS["PLACA"])
    for cx, cy, r in layout.furos_terra:
        out += boundary(_circulo(cx, cy, r), CAMADAS_GDS["TERRA_FUROS"])
    for p in layout.fendas_terra:
        out += boundary(np.asarray(p), CAMADAS_GDS["TERRA_FUROS"])
    out += _rec(0x0700) + _rec(0x0400)
    caminho.write_bytes(bytes(out))
    return [caminho]
