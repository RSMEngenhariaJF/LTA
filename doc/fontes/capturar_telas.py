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


if __name__ == "__main__":
    main()
