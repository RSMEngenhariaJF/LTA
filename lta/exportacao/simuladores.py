"""Geração de modelos 3D editáveis para simuladores EM.

- HFSS (Ansys Electronics Desktop): script IronPython/CPython executado em
  Tools > Run Script. Cria projeto, variáveis de projeto, substrato, trilha,
  terra (com furos), caixa de ar, portas de onda, setup e varredura.
- CST Studio Suite: macro VBA (.bas). Executar em Home > Macros > Edit/Run
  (ou colar no History List). Cria parâmetros, materiais, sólidos, portas,
  fronteiras e faixa de frequência.

Trilha e terra no HFSS são folhas PEC (espessura nula); no CST são sólidos PEC
com espessura t. Ajuste para cobre/rugosidade conforme a necessidade.
"""
from pathlib import Path

import numpy as np


def _fmt(v):
    return f"{v:.6f}".rstrip("0").rstrip(".") if v != 0 else "0"


# ---------------------------------------------------------------------------
# HFSS
# ---------------------------------------------------------------------------
def hfss(layout, caminho, f_min_ghz, f_max_ghz, f_sol_ghz=None, n_pontos=401):
    sub = layout.substrato
    x0, y0, x1, y1 = layout.placa
    h = sub.h_mm
    h_ar = max(10 * h, 5.0)
    f_sol = f_sol_ghz or f_max_ghz
    nome = layout.nome

    L = []
    a = L.append
    a("# -*- coding: utf-8 -*-")
    a(f"# Modelo HFSS gerado por LTA (linhas periodicas) - estrutura '{nome}'")
    a("# Executar no Ansys Electronics Desktop: Tools > Run Script...")
    a("# Unidades: mm. Trilha e terra como folhas PerfectE (sem perdas de condutor).")
    a("import ScriptEnv")
    a('ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")')
    a("oDesktop.RestoreWindow()")
    a("oProject = oDesktop.NewProject()")
    a(f'oProject.InsertDesign("HFSS", "{nome}", "DrivenModal", "")')
    a(f'oDesign = oProject.SetActiveDesign("{nome}")')
    a('oEditor = oDesign.SetActiveEditor("3D Modeler")')
    a('oEditor.SetModelUnits(["NAME:Units Parameter", "Units:=", "mm", "Rescale:=", False])')
    a("")
    a("def var(nome, valor):")
    a('    oDesign.ChangeProperty(["NAME:AllTabs", ["NAME:LocalVariableTab",')
    a('        ["NAME:PropServers", "LocalVariables"],')
    a('        ["NAME:NewProps", ["NAME:" + nome, "PropType:=", "VariableProp", "UserDef:=", True, "Value:=", valor]]]])')
    a("")
    a(f'var("h_sub", "{_fmt(h)}mm")')
    a(f'var("h_ar", "{_fmt(h_ar)}mm")')
    a("")
    a("def attr(nome, material, cor, transp=0.0, solve=True):")
    a('    return ["NAME:Attributes", "Name:=", nome, "Flags:=", "", "Color:=", cor,')
    a('            "Transparency:=", transp, "PartCoordinateSystem:=", "Global", "UDMId:=", "",')
    a('            "MaterialValue:=", "\\"" + material + "\\"", "SolveInside:=", solve]')
    a("")
    a("oDefMgr = oProject.GetDefinitionManager()")
    a('oDefMgr.AddMaterial(["NAME:substrato_lta", "CoordinateSystemType:=", "Cartesian",')
    a('    "BulkOrSurfaceType:=", 1, ["NAME:PhysicsTypes", "set:=", ["Electromagnetic"]],')
    a(f'    "permittivity:=", "{sub.er:g}", "dielectric_loss_tangent:=", "{sub.tand:g}"])')
    a("")
    a("# Substrato")
    a('oEditor.CreateBox(["NAME:BoxParameters",')
    a(f'    "XPosition:=", "{_fmt(x0)}mm", "YPosition:=", "{_fmt(y0)}mm", "ZPosition:=", "0mm",')
    a(f'    "XSize:=", "{_fmt(x1 - x0)}mm", "YSize:=", "{_fmt(y1 - y0)}mm", "ZSize:=", "h_sub"],')
    a('    attr("Substrato", "substrato_lta", "(143 175 143)", 0.6))')
    a("")
    a("# Caixa de ar (contém o substrato)")
    a('oEditor.CreateBox(["NAME:BoxParameters",')
    a(f'    "XPosition:=", "{_fmt(x0)}mm", "YPosition:=", "{_fmt(y0)}mm", "ZPosition:=", "0mm",')
    a(f'    "XSize:=", "{_fmt(x1 - x0)}mm", "YSize:=", "{_fmt(y1 - y0)}mm", "ZSize:=", "h_sub+h_ar"],')
    a('    attr("Ar", "vacuum", "(128 255 255)", 0.9))')
    a("")
    a('def poligono(nome, pts, z="h_sub"):')
    a('    p = [["NAME:PLPoint", "X:=", "%.6fmm" % x, "Y:=", "%.6fmm" % y, "Z:=", z] for x, y in pts]')
    a("    p.append(p[0])")
    a('    seg = [["NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", i, "NoOfPoints:=", 2]')
    a("           for i in range(len(pts))]")
    a('    oEditor.CreatePolyline(["NAME:PolylineParameters", "IsPolylineCovered:=", True,')
    a('        "IsPolylineClosed:=", True, ["NAME:PolylinePoints"] + p, ["NAME:PolylineSegments"] + seg],')
    a('        attr(nome, "vacuum", "(255 128 0)"))')
    a("")
    trilhas = []
    for i, p in enumerate(layout.topo):
        nm = f"Trilha{i + 1}" if len(layout.topo) > 1 else "Trilha"
        trilhas.append(nm)
        pts = ", ".join(f"({x:.6f}, {y:.6f})" for x, y in p)
        a(f'poligono("{nm}", [{pts}])')
    a("")
    a("# Plano de terra")
    a('oEditor.CreateRectangle(["NAME:RectangleParameters", "IsCovered:=", True,')
    a(f'    "XStart:=", "{_fmt(x0)}mm", "YStart:=", "{_fmt(y0)}mm", "ZStart:=", "0mm",')
    a(f'    "Width:=", "{_fmt(x1 - x0)}mm", "Height:=", "{_fmt(y1 - y0)}mm", "WhichAxis:=", "Z"],')
    a('    attr("Terra", "vacuum", "(255 200 100)"))')
    if layout.furos_terra or layout.fendas_terra:
        furos = []
        for i, (cx, cy, r) in enumerate(layout.furos_terra):
            nm = f"Furo{i + 1}"
            furos.append(nm)
            a('oEditor.CreateCircle(["NAME:CircleParameters", "IsCovered:=", True,')
            a(f'    "XCenter:=", "{_fmt(cx)}mm", "YCenter:=", "{_fmt(cy)}mm", "ZCenter:=", "0mm",')
            a(f'    "Radius:=", "{_fmt(r)}mm", "WhichAxis:=", "Z", "NumSegments:=", "0"],')
            a(f'    attr("{nm}", "vacuum", "(0 0 255)"))')
        for i, p in enumerate(layout.fendas_terra):
            nm = f"Fenda{i + 1}"
            furos.append(nm)
            pts = ", ".join(f"({x:.6f}, {y:.6f})" for x, y in p)
            a(f'poligono("{nm}", [{pts}], "0mm")')
        a(f'oEditor.Subtract(["NAME:Selections", "Blank Parts:=", "Terra", "Tool Parts:=", "{",".join(furos)}"],')
        a('    ["NAME:SubtractParameters", "KeepOriginals:=", False])')
    a("")
    vias = []
    for i, (cx, cy, r) in enumerate(layout.vias):
        nm = f"Via{i + 1}"
        vias.append(nm)
        a('oEditor.CreateCylinder(["NAME:CylinderParameters",')
        a(f'    "XCenter:=", "{_fmt(cx)}mm", "YCenter:=", "{_fmt(cy)}mm", "ZCenter:=", "0mm",')
        a(f'    "Radius:=", "{_fmt(r)}mm", "Height:=", "h_sub", "WhichAxis:=", "Z", "NumSides:=", "0"],')
        a(f'    attr("{nm}", "vacuum", "(160 160 160)"))')
    a("")
    a('oBnd = oDesign.GetModule("BoundarySetup")')
    a(f'oBnd.AssignPerfectE(["NAME:PEC_Cobre", "Objects:=", {trilhas + ["Terra"] + vias!r}, '
      '"InfGroundPlane:=", False])')
    a("")
    a("# Radiacao nas faces laterais (y) e superior (z) da caixa de ar")
    a("def face(x, y, z):")
    a('    return oEditor.GetFaceByPosition(["NAME:FaceParameters", "BodyName:=", "Ar",')
    a('        "XPosition:=", "%smm" % x, "YPosition:=", "%smm" % y, "ZPosition:=", z])')
    xm, ym = _fmt((x0 + x1) / 2), _fmt((y0 + y1) / 2)
    a(f'faces_rad = [face("{xm}", "{_fmt(y0)}", "(h_sub+h_ar)/2"), face("{xm}", "{_fmt(y1)}", "(h_sub+h_ar)/2"),')
    a(f'             face("{xm}", "{ym}", "h_sub+h_ar")]')
    a('oBnd.AssignRadiation(["NAME:Rad1", "Faces:=", faces_rad, "IsFssReference:=", False, "IsForPML:=", False])')
    a("")
    a("# Portas de onda nas faces x da caixa")
    for k, porta in enumerate(layout.portas, start=1):
        wp = min(y1 - y0, max(10 * porta.largura, 10 * h))
        hp = min(h + h_ar, 10 * h)
        a('oEditor.CreateRectangle(["NAME:RectangleParameters", "IsCovered:=", True,')
        a(f'    "XStart:=", "{_fmt(porta.x)}mm", "YStart:=", "{_fmt(porta.y - wp / 2)}mm", "ZStart:=", "0mm",')
        a(f'    "Width:=", "{_fmt(wp)}mm", "Height:=", "{_fmt(hp)}mm", "WhichAxis:=", "X"],')
        a(f'    attr("Porta{k}", "vacuum", "(255 0 0)"))')
        a(f'oBnd.AssignWavePort(["NAME:{k}", "Objects:=", ["Porta{k}"], "NumModes:=", 1,')
        a('    "RenormalizeAllTerminals:=", True, "UseLineModeAlignment:=", False, "DoDeembed:=", False,')
        a(f'    ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", True,')
        a(f'        ["NAME:IntLine", "Start:=", ["{_fmt(porta.x)}mm", "{_fmt(porta.y)}mm", "0mm"],')
        a(f'                         "End:=", ["{_fmt(porta.x)}mm", "{_fmt(porta.y)}mm", "h_sub"]],')
        a('        "AlignmentGroup:=", 0, "CharImp:=", "Zpi", "RenormImp:=", "50ohm"]],')
        a('    "ShowReporterFilter:=", False, "ReporterFilter:=", [True], "UseAnalyticAlignment:=", False])')
    a("")
    a('oAna = oDesign.GetModule("AnalysisSetup")')
    a(f'oAna.InsertSetup("HfssDriven", ["NAME:Setup1", "Frequency:=", "{f_sol:g}GHz",')
    a('    "MaxDeltaS:=", 0.02, "MaximumPasses:=", 15, "MinimumPasses:=", 2, "MinimumConvergedPasses:=", 2,')
    a('    "PercentRefinement:=", 30, "IsEnabled:=", True, "BasisOrder:=", 1])')
    a('oAna.InsertFrequencySweep("Setup1", ["NAME:Varredura", "IsEnabled:=", True,')
    a(f'    "RangeType:=", "LinearCount", "RangeStart:=", "{f_min_ghz:g}GHz", "RangeEnd:=", "{f_max_ghz:g}GHz",')
    a(f'    "RangeCount:=", {int(n_pontos)}, "Type:=", "Interpolating", "SaveFields:=", False,')
    a('    "SaveRadFields:=", False, "InterpTolerance:=", 0.5, "InterpMaxSolns:=", 250,')
    a('    "InterpMinSolns:=", 0, "InterpMinSubranges:=", 1, "ExtrapToDC:=", False])')
    a('oEditor.FitAll()')
    caminho.write_text("\n".join(L) + "\n", encoding="utf-8")
    return [caminho]


# ---------------------------------------------------------------------------
# CST Studio Suite
# ---------------------------------------------------------------------------
def cst(layout, caminho, f_min_ghz, f_max_ghz):
    sub = layout.substrato
    x0, y0, x1, y1 = layout.placa
    h = sub.h_mm
    t = sub.t_um / 1000

    L = []
    a = L.append
    a(f"' Macro CST gerada por LTA (linhas periodicas) - estrutura '{layout.nome}'")
    a("' Executar em Home > Macros > Edit/Run VBA Macro (ou colar no History List).")
    a("' Unidades: mm / GHz. Cobre modelado como PEC com espessura t.")
    a("Sub Main ()")
    a("With Units")
    a('    .Geometry "mm"')
    a('    .Frequency "GHz"')
    a('    .Time "ns"')
    a("End With")
    a(f'StoreDoubleParameter "h", {_fmt(h)}')
    a(f'StoreDoubleParameter "t", {_fmt(t)}')
    a(f'StoreDoubleParameter "er", {sub.er:g}')
    a(f'StoreDoubleParameter "tand", {sub.tand:g}')
    a("RebuildOnParametricChange False, False")
    a("")
    a("With Material")
    a("    .Reset")
    a('    .Name "substrato_lta"')
    a('    .Folder ""')
    a('    .FrqType "all"')
    a('    .Type "Normal"')
    a('    .Epsilon "er"')
    a('    .Mu "1"')
    a('    .TanD "tand"')
    a(f'    .TanDFreq "{(f_min_ghz + f_max_ghz) / 2:g}"')
    a('    .TanDGiven "True"')
    a('    .TanDModel "ConstTanD"')
    a('    .Colour "0.5", "0.8", "0.5"')
    a("    .Create")
    a("End With")
    a("")
    a("With Brick")
    a('    .Reset : .Name "Substrato" : .Component "lta" : .Material "substrato_lta"')
    a(f'    .Xrange "{_fmt(x0)}", "{_fmt(x1)}" : .Yrange "{_fmt(y0)}", "{_fmt(y1)}" : .Zrange "0", "h"')
    a("    .Create")
    a("End With")
    a("With Brick")
    a('    .Reset : .Name "Terra" : .Component "lta" : .Material "PEC"')
    a(f'    .Xrange "{_fmt(x0)}", "{_fmt(x1)}" : .Yrange "{_fmt(y0)}", "{_fmt(y1)}" : .Zrange "-t", "0"')
    a("    .Create")
    a("End With")
    for i, (cx, cy, r) in enumerate(layout.furos_terra, start=1):
        a("With Cylinder")
        a(f'    .Reset : .Name "Furo{i}" : .Component "lta" : .Material "PEC"')
        a(f'    .OuterRadius "{_fmt(r)}" : .InnerRadius "0" : .Axis "z" : .Zrange "-t", "0"')
        a(f'    .Xcenter "{_fmt(cx)}" : .Ycenter "{_fmt(cy)}" : .Segments "0"')
        a("    .Create")
        a("End With")
        a(f'Solid.Subtract "lta:Terra", "lta:Furo{i}"')
    a("")
    a("' Trilha: curva poligonal no plano z = h, extrudada com espessura t")
    a("WCS.ActivateWCS \"local\"")
    a('WCS.SetNormal "0", "0", "1"')
    a('WCS.SetOrigin "0", "0", "h"')
    a('WCS.SetUVector "1", "0", "0"')
    for i, p in enumerate(layout.topo, start=1):
        a("With Polygon")
        a(f'    .Reset : .Name "contorno{i}" : .Curve "curvas"')
        a(f"    .Point {p[0][0]:.6f}, {p[0][1]:.6f}")
        for x, y in p[1:]:
            a(f"    .LineTo {x:.6f}, {y:.6f}")
        a(f"    .LineTo {p[0][0]:.6f}, {p[0][1]:.6f}")
        a("    .Create")
        a("End With")
        a("With ExtrudeCurve")
        a(f'    .Reset : .Name "Trilha{i}" : .Component "lta" : .Material "PEC"')
        a('    .Thickness "t" : .Twistangle "0" : .Taperangle "0" : .DeleteProfile "True"')
        a(f'    .Curve "curvas:contorno{i}"')
        a("    .Create")
        a("End With")
    if layout.fendas_terra:
        a("' Fendas no terra (ex.: CSRR): curvas em z = -t extrudadas com t e subtraídas")
        a('WCS.SetOrigin "0", "0", "-t"')
        for i, p in enumerate(layout.fendas_terra, start=1):
            a("With Polygon")
            a(f'    .Reset : .Name "fenda{i}" : .Curve "curvas"')
            a(f"    .Point {p[0][0]:.6f}, {p[0][1]:.6f}")
            for x, y in p[1:]:
                a(f"    .LineTo {x:.6f}, {y:.6f}")
            a(f"    .LineTo {p[0][0]:.6f}, {p[0][1]:.6f}")
            a("    .Create")
            a("End With")
            a("With ExtrudeCurve")
            a(f'    .Reset : .Name "Fenda{i}" : .Component "lta" : .Material "PEC"')
            a('    .Thickness "t" : .Twistangle "0" : .Taperangle "0" : .DeleteProfile "True"')
            a(f'    .Curve "curvas:fenda{i}"')
            a("    .Create")
            a("End With")
            a(f'Solid.Subtract "lta:Terra", "lta:Fenda{i}"')
    a('WCS.ActivateWCS "global"')
    for i, (cx, cy, r) in enumerate(layout.vias, start=1):
        a("With Cylinder")
        a(f'    .Reset : .Name "Via{i}" : .Component "lta" : .Material "PEC"')
        a(f'    .OuterRadius "{_fmt(r)}" : .InnerRadius "0" : .Axis "z" : .Zrange "0", "h"')
        a(f'    .Xcenter "{_fmt(cx)}" : .Ycenter "{_fmt(cy)}" : .Segments "0"')
        a("    .Create")
        a("End With")
    a("")
    a("With Background")
    a('    .Type "Normal"')
    a('    .XminSpace "0" : .XmaxSpace "0" : .YminSpace "0" : .YmaxSpace "0"')
    a('    .ZminSpace "0" : .ZmaxSpace "10*h"')
    a("End With")
    a("With Boundary")
    a('    .Xmin "expanded open" : .Xmax "expanded open"')
    a('    .Ymin "expanded open" : .Ymax "expanded open"')
    a('    .Zmin "electric" : .Zmax "expanded open"')
    a("End With")
    for k, porta in enumerate(layout.portas, start=1):
        wp = min(y1 - y0, max(10 * porta.largura, 10 * h)) / 2
        a("With Port")
        a(f'    .Reset : .PortNumber "{k}" : .Label "" : .NumberOfModes "1"')
        a('    .AdjustPolarization "False" : .PolarizationAngle "0.0" : .ReferencePlaneDistance "0"')
        a(f'    .TextSize "50" : .Coordinates "Free" : .Orientation "{porta.lado}"')
        a('    .PortOnBound "True" : .ClipPickedPortToBound "False"')
        a(f'    .Xrange "{_fmt(porta.x)}", "{_fmt(porta.x)}"')
        a(f'    .Yrange "{_fmt(porta.y - wp)}", "{_fmt(porta.y + wp)}"')
        a('    .Zrange "-t", "8*h"')
        a("    .Create")
        a("End With")
    a(f'Solver.FrequencyRange "{f_min_ghz:g}", "{f_max_ghz:g}"')
    a("End Sub")
    caminho.write_text("\r\n".join(L) + "\r\n", encoding="utf-8")
    return [caminho]


# ---------------------------------------------------------------------------
# Touchstone
# ---------------------------------------------------------------------------
def touchstone(caminho, f_hz, s, z0=50.0, comentario=""):
    """Grava .s2p (formato RI, frequência em GHz). s: (Nf, 2, 2)."""
    s = np.asarray(s)
    linhas = [f"! {c}" for c in comentario.splitlines() if c] + ["! Gerado por LTA - linhas periodicas",
                                                                 f"# GHz S RI R {z0:g}"]
    for fi, si in zip(f_hz, s):
        # ordem Touchstone v1 para 2 portas: S11 S21 S12 S22
        vals = [si[0, 0], si[1, 0], si[0, 1], si[1, 1]]
        linhas.append(f"{fi / 1e9:.9g} " + " ".join(f"{v.real:.9e} {v.imag:.9e}" for v in vals))
    caminho.write_text("\n".join(linhas) + "\n", encoding="ascii")
    return [caminho]


def ler_touchstone(caminho):
    """Lê .s2p (Touchstone v1): retorna (f [Hz], S (Nf, 2, 2), z0).

    Aceita unidades Hz/kHz/MHz/GHz, formatos RI, MA e DB, comentários '!' e
    dados quebrados em várias linhas. Só parâmetros S de 2 portas.
    """
    escala = {"HZ": 1.0, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}
    unid, fmt, z0 = 1e9, "MA", 50.0
    nums = []
    for linha in Path(caminho).read_text(encoding="latin-1").splitlines():
        linha = linha.split("!", 1)[0].strip()
        if not linha:
            continue
        if linha.startswith("#"):
            toks = linha[1:].upper().split()
            for i, t in enumerate(toks):
                if t in escala:
                    unid = escala[t]
                elif t in ("RI", "MA", "DB"):
                    fmt = t
                elif t == "R" and i + 1 < len(toks):
                    z0 = float(toks[i + 1])
                elif t in ("Y", "Z", "H", "G"):
                    raise ValueError(f"Parâmetros {t} não suportados; exporte parâmetros S.")
            continue
        if linha.startswith("["):
            raise ValueError("Touchstone v2 com palavras-chave não suportado; exporte em v1 (.s2p).")
        nums += [float(x) for x in linha.split()]
    if not nums or len(nums) % 9:
        raise ValueError("Arquivo não parece um .s2p de 2 portas (9 números por frequência).")
    d = np.array(nums).reshape(-1, 9)
    a, b = d[:, 1::2], d[:, 2::2]
    if fmt == "RI":
        v = a + 1j * b
    elif fmt == "MA":
        v = a * np.exp(1j * np.radians(b))
    else:
        v = 10 ** (a / 20) * np.exp(1j * np.radians(b))
    s = np.empty((len(d), 2, 2), dtype=complex)
    s[:, 0, 0], s[:, 1, 0], s[:, 0, 1], s[:, 1, 1] = v[:, 0], v[:, 1], v[:, 2], v[:, 3]
    return d[:, 0] * unid, s, z0


def csv_s(caminho, f_hz, s):
    """CSV com módulo [dB] e fase [graus] de S11, S21, S12, S22 (separador ';')."""
    s = np.asarray(s)
    nomes = (("S11", 0, 0), ("S21", 1, 0), ("S12", 0, 1), ("S22", 1, 1))
    cab = ["f_GHz"] + [f"{n}_{q}" for n, _, _ in nomes for q in ("dB", "graus")]
    linhas = [";".join(cab)]
    for fi, si in zip(f_hz, s):
        vals = [f"{fi / 1e9:.9g}"]
        for _, i, j in nomes:
            v = si[i, j]
            vals += [f"{20 * np.log10(max(abs(v), 1e-15)):.6f}", f"{np.degrees(np.angle(v)):.4f}"]
        linhas.append(";".join(vals))
    caminho.write_text("\n".join(linhas) + "\n", encoding="ascii")
    return [caminho]
