"""Captura telas da interface para o manual (doc/fontes/img/tela_*.png).

Uso (na raiz do projeto):  python doc/fontes/capturar_telas.py
A janela do programa abre, executa os exemplos do livro e fecha sozinha.
"""
import ctypes
import ctypes.wintypes
import sys
import time
import tkinter as tk
from pathlib import Path
from tkinter import messagebox

RAIZ = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RAIZ))
IMG = Path(__file__).resolve().parent / "img"

from PIL import ImageGrab  # noqa: E402

from lta.gui import app  # noqa: E402


def main():
    try:  # coordenadas em pixels físicos (telas com escala > 100%)
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except (AttributeError, OSError):
        pass
    IMG.mkdir(parents=True, exist_ok=True)
    messagebox.showerror = lambda t, m: print("ERRO:", t, m)

    raiz = tk.Tk()
    escala = raiz.winfo_fpixels("1i") / 96
    raiz.geometry(f"{int(1400 * escala)}x{int(800 * escala)}+20+10")
    raiz.attributes("-topmost", True)
    aplic = app.Aplicacao(raiz)
    sub, cel, car, ebg, exp = (aplic.ferramenta(c) for c in (app.AbaSubstrato, app.AbaCelula, app.AbaCarga,
                                                               app.AbaEBG, app.AbaExportar))

    class _Nav:  # mantém a sintaxe nb.select(ferramenta) usada abaixo
        @staticmethod
        def select(ferramenta):
            aplic.mostrar(type(ferramenta))

    nb = _Nav()

    def esperar(t=1.2):
        t0 = time.time()
        while time.time() - t0 < t:
            raiz.update()
            time.sleep(0.02)

    def foto(nome):
        esperar()
        # janela inteira (inclui título e barra de menu)
        ret = ctypes.wintypes.RECT()
        ctypes.windll.user32.GetWindowRect(int(raiz.wm_frame(), 16), ctypes.byref(ret))
        x, y, w, h = ret.left, ret.top, ret.right - ret.left, ret.bottom - ret.top
        ImageGrab.grab(bbox=(x, y, x + w, y + h), all_screens=True).save(IMG / f"tela_{nome}.png")
        print("capturada:", nome)

    esperar()
    nb.select(sub)
    sub.calc.set({"z0": 50})
    sub.z_para_w()
    foto("01_substrato")

    nb.select(cel)
    cel.form.set(dict(topologia="T", comparar="L", ls=2.5, cs=0, lp=0, cp=1.0, n=20, zp=50,
                      fmin=0.01, fmax=8, npts=1001))
    cel.atualizar_previa()
    cel.ver_estrutura()
    foto("02_celula_estrutura")
    cel.calcular()
    foto("03_celula_resultados")
    cel.exemplo_trelica()
    foto("13_trelica_resultados")
    cel.exemplo_ponte()
    cel.ver_estrutura()
    foto("15_ponte_estrutura")

    nb.select(car)
    car.exemplo()
    car.ver_estrutura()
    esperar(0.3)
    car.fis.vars["zpatch"].set("15")   # edição simulada: mostra referência e diferenças
    esperar(0.8)
    foto("04_carga_estrutura")
    car.fis.set({"zpatch": 10})
    car.calcular()
    foto("05_carga_resultados")

    nb.select(ebg)
    ebg.ex_sen()
    foto("06_ebg_sen_resultados")
    ebg.ex_quad()
    ebg.ver_estrutura()
    esperar(0.3)
    ebg.form.vars["k"].set("0.6")      # edição simulada
    esperar(0.8)
    foto("07_ebg_quad_estrutura")
    ebg.form.set({"k": 1.0})
    ebg.calcular()
    foto("08_ebg_quad_resultados")
    ebg.ex_furos()
    foto("09_ebg_furos_estrutura")

    ebg.modelo.set("Duas bandas: 3 e 5 GHz")
    ebg.carregar_modelo()
    ebg.form.set({"fmin": 0.5, "fmax": 8, "npts": 1000})
    esperar(0.3)
    ebg.expr.delete("1.0", "end")   # edição simulada: 2a banda de 5 para 5,5 GHz
    ebg.expr.insert("1.0", "0.4*cos(2*pi*z/bragg(3)) + 0.4*cos(2*pi*z/bragg(5.5))")
    ebg._agendar_previa()
    esperar(0.8)
    foto("11_ebg_pers_estrutura")
    ebg.calcular()
    foto("12_ebg_pers_resultados")

    # ---- perdas e dispersão (EBG da Fig. 2.11 com o modelo físico completo)
    ebg.ex_sen()
    ebg.modelo_fis.set(True, True)
    ebg.calcular()
    foto("37_ebg_perdas")
    ebg.modelo_fis.set(False, False)

    capturar_novas(aplic, foto, esperar)

    ebg.ex_sen()
    ebg.publicar()
    nb.select(exp)
    exp.pasta.set("saidas/exemplos_manual")
    exp.exportar()
    exp.sel_s.set("ebg")
    exp.exportar_s()
    exp.exportar()
    foto("10_exportar")

    raiz.destroy()


def capturar_novas(aplic, foto, esperar):
    """Telas das ferramentas acrescentadas a partir da v1.2 (caps. 0 e 2 a 6 do livro)."""
    from lta.core import linhas_planares as lp
    from lta.gui import aba_dispositivos, aba_extracao, aba_filtros, aba_sintonia
    from lta.gui.aba_substrato import PREDEFINIDOS

    def ir(cls):
        aplic.mostrar(cls)
        aplic.estado.msg("")   # não deixa a mensagem de outra ferramenta na barra de status
        return aplic.ferramenta(cls)

    # Cap. 0: linhas convencionais (RO4003C)
    aplic.estado.definir_substrato(PREDEFINIDOS["Rogers RO4003C 0.813 mm"])
    lin = ir(app.AbaLinhas)
    lin.form.set(dict(tipo=lp.TIPOS[3], elemento="Trecho de linha", modo="Síntese: Z0 → dimensões", z0=50,
                      g=0.2, lterra=3, f0=2.45, fase=90, pvia=1.5, rvia=0.2))
    lin.calcular()
    lin.ver_estrutura()
    foto("17_linhas_gcpw_estrutura")
    lin.form.set(dict(tipo=lp.TIPOS[0], elemento="Transformador λ/4", z0=50, z2=100, fmax=7.5))
    lin.calcular()
    foto("18_linhas_transformador_resultados")

    # Cap. 2: Floquet e aplicações EBG
    flo = ir(app.AbaFloquet)
    flo.ex_ebg()
    foto("19_floquet")
    apl = ir(app.AbaAplicacoesEBG)
    apl.ex_cdl()
    foto("20_cdl")
    apl.ex_res()
    foto("21_ressonador")

    # Cap. 3: metamateriais
    cr = ir(app.AbaCRLH)
    cr.ex_desbalanceada()
    foto("22_crlh")
    cl = ir(app.AbaCLLoaded)
    cl.exemplo()
    cl.ver_estrutura()
    foto("23_cl_loaded_estrutura")
    cs = ir(app.AbaCSRR)
    cs.ex_balanceada()
    cs.ver_estrutura()
    foto("24_csrr_estrutura")
    cs.calcular()
    foto("25_csrr_resultados")
    ex = ir(app.AbaExtracao)
    ex.fonte.set(next(iter(aba_extracao.EXEMPLOS)))
    ex._fonte_mudou()
    ex.extrair()
    foto("26_extracao")

    # Cap. 4: aplicações
    mb = ir(app.AbaMultibanda)
    mb.ex_quad()
    foto("27_multibanda_quad")
    dv = ir(app.AbaDispositivos)
    dv.ex_ratrace()
    foto("28_ratrace")
    ft = ir(app.AbaFiltros)
    ft.ex_hibrido()
    foto("29_filtro_hibrido")

    # Cap. 5: sintonizáveis e não lineares
    si = ir(app.AbaSintonia)
    si.ex_vlsrr()
    foto("30_vlsrr")
    si.ex_crlh()
    foto("31_defasador")
    nlt = ir(app.AbaNLTL)
    nlt.ex_pulso()
    foto("32_nltl_pulso")

    # Cap. 6: outras linhas avançadas
    di = ir(app.AbaDiferencial)
    di.form.set({"preset": "DS-CSRR 3 GHz"})
    di._carregar()
    di.calcular()
    foto("33_diferencial")
    mi = ir(app.AbaMIW)
    mi.exemplo()
    foto("34_miw")
    sw = ir(app.AbaSIW)
    sw.calcular()
    sw.ver_estrutura()
    foto("35_siw_estrutura")
    nf = ir(app.AbaNaoFoster)
    nf.exemplo()
    foto("36_naofoster")
    aplic.estado.definir_substrato(PREDEFINIDOS["Rogers RO3010 (livro)"])
    _ = (aba_dispositivos, aba_filtros, aba_sintonia, esperar)


if __name__ == "__main__":
    main()
