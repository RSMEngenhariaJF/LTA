"""Registro dos formatos de exportação.

Para adicionar um formato novo: escreva a função em um módulo deste pacote e
inclua uma entrada em FORMATOS (rótulo -> (extensão, função, destinos)).
"""
from pathlib import Path

from . import formatos_2d, simuladores

# rótulo: (extensão, função(layout, caminho, opcoes) -> [Path], ferramentas de destino)
FORMATOS = {
    "KiCad (.kicad_pcb)": (".kicad_pcb", lambda lay, c, o: formatos_2d.kicad(lay, c), "KiCad, Altium (Import Wizard)"),
    "DXF R12 (.dxf)": (".dxf", lambda lay, c, o: formatos_2d.dxf(lay, c), "Altium, ADS, HFSS, CST, KiCad"),
    "Gerber RS-274X": (".gbr", lambda lay, c, o: formatos_2d.gerber(lay, c), "Altium (CAMtastic), fabricação"),
    "GDSII (.gds)": (".gds", lambda lay, c, o: formatos_2d.gds(lay, c), "ADS, HFSS, CST"),
    "HFSS - script (.py)": ("_hfss.py", lambda lay, c, o: simuladores.hfss(
        lay, c, o["f_min_ghz"], o["f_max_ghz"], o.get("f_sol_ghz")), "Ansys HFSS / AEDT"),
    "CST - macro VBA (.bas)": ("_cst.bas", lambda lay, c, o: simuladores.cst(
        lay, c, o["f_min_ghz"], o["f_max_ghz"]), "CST Studio Suite"),
}

INSTRUCOES = """Arquivos gerados por LTA - linhas periodicas
Estrutura: {nome}
Substrato: {sub.nome}  er = {sub.er:g}  h = {sub.h_mm:g} mm  t = {sub.t_um:g} um  tan(d) = {sub.tand:g}
Placa: {comp:.3f} mm x {larg:.3f} mm (x = propagacao, porta 1 em x = 0). Unidades: mm.

Camadas (DXF / GDSII):
  TOPO_COBRE  (GDS 1)  trilha
  TERRA       (GDS 2)  plano de terra (face inferior)
  TERRA_FUROS (GDS 3)  furos a subtrair do plano de terra
  PLACA       (GDS 10) contorno do substrato

Como abrir:
  KiCad   : abrir o .kicad_pcb no Editor de PCB; pressionar B para preencher a zona de terra.
  Altium  : (a) File > Import > DXF/DWG no PCB, mapeando TOPO_COBRE -> Top Layer e TERRA -> Bottom Layer;
            (b) File > Import Wizard > KiCad (abrir o .kicad_pcb); ou
            (c) importar os Gerber (.gtl/.gbl/.gko) no CAMtastic e exportar para PCB.
  ADS     : Layout > File > Import > GDSII (ou DXF), unidade 1 nm / um; mapear camadas 1 e 2 para
            cond e o plano de terra do substrato; o .s2p (se exportado) entra em um bloco S2P.
  HFSS    : Tools > Run Script > *_hfss.py (cria o projeto completo, com variaveis h_sub e h_ar).
  CST     : Home > Macros > Edit/Run VBA Macro > *_cst.bas (parametros h, t, er, tand editaveis).

Parametros de projeto:
{info}

ATENCAO: geometria e uma estimativa inicial a partir de modelos de circuito/quase-estaticos.
Validar por simulacao EM e medicao, com revisao tecnica, antes de qualquer uso em processo controlado.
"""


def exportar(layout, pasta, formatos, opcoes):
    pasta = Path(pasta)
    pasta.mkdir(parents=True, exist_ok=True)
    gerados = []
    for rotulo in formatos:
        ext, func, _ = FORMATOS[rotulo]
        gerados += func(layout, pasta / f"{layout.nome}{ext}", opcoes)
    info = "\n".join(f"  {k}: {v}" for k, v in layout.info.items()) or "  (nenhum)"
    leia = pasta / f"{layout.nome}_LEIA-ME.txt"
    leia.write_text(INSTRUCOES.format(nome=layout.nome, sub=layout.substrato, comp=layout.comprimento,
                                      larg=layout.largura, info=info), encoding="utf-8")
    return gerados + [leia]
