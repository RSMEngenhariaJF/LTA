"""Aba EBG: linhas não uniformes analisadas por modos acoplados (seção 2.4)."""
import tkinter as tk
from tkinter import ttk

import numpy as np

from ..core import abcd
from ..core import modos_acoplados as cm
from ..core import perfil_personalizado as pp
from ..core.microstrip import analise, sintese
from ..geometria import construtores
from .aba_substrato import PREDEFINIDOS
from .comum import AbaBase, AreaTexto, Formulario, protegido
from .previa import PreviaEstrutura

PERTURBACOES = ("Senoidal (largura)", "Quadrada (largura)", "Furos no terra", "Personalizado (função)")
PERSONALIZADO = PERTURBACOES[3]
AJUDA = ("z em mm (0 a L). Variáveis: z, L, Zref. Funções: sin, cos, exp, log, sqrt, abs, tanh, "
         "sign, where, minimum, maximum, clip, pi.\n"
         "bragg(f) = período [mm] da banda em f GHz.  janela(z) = Hamming em [0, L].  "
         "hann(z), ret(z, a, b), degrau(z, a).\n"
         "Em K(z) [1/cm], cada termo A·cos(2·pi·z/bragg(f)) cria uma banda em f "
         "(|S21|min ≈ sech(A·L/2)). A média de K é removida para a linha terminar em Zref.")


class AbaEBG(AbaBase):
    TITULO = "EBG / Modos acoplados"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("tipo", "Perturbação", PERTURBACOES, ""),
            ("periodo", "Período l", 10.0, "mm"),
            ("n", "Nº de períodos", 10, ""),
            ("k", "|K| (amplitude)", 0.64, "1/cm"),
            ("zref", "Z referência", 50.0, "Ω"),
            ("r_l", "Raio/período (furos)", 0.25, ""),
            ("lalim", "Alimentação", 5.0, "mm"),
            ("fmin", "f mín", 0.1, "GHz"),
            ("fmax", "f máx", 18.0, "GHz"),
            ("npts", "Pontos", 800, ""),
        ], titulo="Estrutura periódica")
        self.form.pack(fill="x")
        self._painel_personalizado()
        self.botoes(("Calcular", self.calcular), ("Ex. Fig. 2.11 (sen)", self.ex_sen),
                    ("Ex. Fig. 2.14 (quad)", self.ex_quad), ("Ex. Fig. 2.16 (furos)", self.ex_furos),
                    ("Enviar p/ exportação", self.publicar), ("Exportar S2P", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 16)
        self.saida.pack(fill="both", expand=True)
        self._ultimo = None
        self._lay = None
        self.observar_previa(self.form, self.form_pers)
        self.expr.bind("<KeyRelease>", lambda e: self._agendar_previa())
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    def _painel_personalizado(self):
        caixa = ttk.LabelFrame(self.painel, text="Perfil personalizado (função de z)", padding=6)
        caixa.pack(fill="x", pady=(8, 0))
        linha = ttk.Frame(caixa)
        linha.pack(fill="x")
        ttk.Label(linha, text="Modelo").pack(side="left")
        self.modelo = tk.StringVar(value=next(iter(pp.MODELOS)))
        ttk.Combobox(linha, textvariable=self.modelo, values=list(pp.MODELOS), state="readonly",
                     width=18).pack(side="left", padx=4, fill="x", expand=True)
        ttk.Button(linha, text="Carregar", command=self.carregar_modelo).pack(side="left")
        self.form_pers = Formulario(caixa, [
            ("func", "A função define", pp.TIPOS_FUNCAO, ""),
            ("comp", "Comprimento L", 150.0, "mm"),
            ("passo", "Passo dz", 0.2, "mm"),
        ])
        self.form_pers.pack(fill="x", pady=(4, 2))
        ttk.Label(caixa, text="Expressão:").pack(anchor="w")
        self.expr = tk.Text(caixa, height=4, width=36, wrap="word", font=("Consolas", 9))
        self.expr.pack(fill="x")
        self.expr.insert("1.0", pp.MODELOS[self.modelo.get()][1])
        self.desc_modelo = ttk.Label(caixa, text=pp.MODELOS[self.modelo.get()][4], foreground="#1f4e9a",
                                     wraplength=300)
        self.desc_modelo.pack(anchor="w", pady=(2, 0))
        ttk.Label(caixa, text=AJUDA, foreground="#555", wraplength=300, font=("Segoe UI", 8)).pack(anchor="w")

    def carregar_modelo(self):
        tipo, expr, comp, passo, desc = pp.MODELOS[self.modelo.get()]
        self.form.set({"tipo": PERSONALIZADO})
        self.form_pers.set({"func": tipo, "comp": comp, "passo": passo})
        self.expr.delete("1.0", "end")
        self.expr.insert("1.0", expr)
        self.desc_modelo.configure(text=desc)
        self.ver_estrutura()
        self.atualizar_previa()

    def _ro3010(self):
        self.estado.definir_substrato(PREDEFINIDOS["Rogers RO3010 (livro)"])

    def ex_sen(self):
        self._ro3010()
        self.form.set(dict(tipo=PERTURBACOES[0], periodo=10, n=10, k=0.64, zref=50, fmin=0.1, fmax=18))
        self.calcular()

    def ex_quad(self):
        self._ro3010()
        self.form.set(dict(tipo=PERTURBACOES[1], periodo=10, n=10, k=1.0, zref=50, fmin=0.1, fmax=18))
        self.calcular()

    def ex_furos(self):
        self._ro3010()
        self.form.set(dict(tipo=PERTURBACOES[2], periodo=18.9, n=7, r_l=0.25, zref=50, fmin=0.1, fmax=6))
        self.calcular()

    def _gerar(self):
        """Lê os campos e gera o layout (e o perfil, nos casos de largura modulada)."""
        v = self.form.get()
        if v["n"] < 1 or v["periodo"] <= 0:
            raise ValueError("Período e número de períodos devem ser positivos.")
        sub = self.estado.substrato
        if v["tipo"] == PERTURBACOES[2]:
            if not 0 < v["r_l"] < 0.5:
                raise ValueError("Raio/período deve estar entre 0 e 0,5 (furos não podem se sobrepor).")
            w = sintese(v["zref"], sub)
            _, ee = analise(w, sub)
            r = v["r_l"] * v["periodo"]
            fb = cm.frequencia_bragg(v["periodo"] * 1e-3, ee)
            info = {"tipo": "EBG com furos no plano de terra", "periodo_mm": v["periodo"], "furos": v["n"],
                    "raio_mm": f"{r:.4f}", "W_mm": f"{w:.4f}", "f_Bragg_estimada_GHz": f"{fb / 1e9:.3f}"}
            lay = construtores.linha_furos_terra(v["periodo"], r, v["n"], v["zref"], sub, v["lalim"],
                                                 nome=f"ebg_furos_{v['n']}f", info=info)
            return v, lay, dict(w=w, ee=ee, r=r, fb=fb)
        if v["tipo"] == PERSONALIZADO:
            return self._gerar_personalizado(v, sub)
        tipo = "Senoidal" if v["tipo"] == PERTURBACOES[0] else "Quadrada"
        p = cm.perfil(tipo, v["k"] * 100, v["periodo"] * 1e-3, v["n"], v["zref"])
        w, ee = cm.larguras(p, sub)
        ee_med = cm.eeff_medio(p.z, ee)
        f_b = cm.frequencia_bragg(p.periodo, ee_med)
        info = {"tipo": f"EBG largura modulada ({tipo})", "periodo_mm": v["periodo"], "periodos": v["n"],
                "K_1_cm": v["k"], "Zref_ohm": v["zref"], "Z0_min_max_ohm": f"{p.z0.min():.1f}/{p.z0.max():.1f}",
                "eeff_medio": f"{ee_med:.4f}", "f_Bragg_GHz": f"{f_b / 1e9:.3f}"}
        lay = construtores.linha_largura_modulada(p.z * 1e3, w, sub, v["lalim"],
                                                  nome=f"ebg_{tipo.lower()}_{v['n']}p", info=info,
                                                  periodo=v["periodo"])
        return v, lay, dict(tipo=tipo, p=p, w=w, ee=ee, ee_med=ee_med, avisos=[])

    def _gerar_personalizado(self, v, sub):
        vp = self.form_pers.get()
        expr = self.expr.get("1.0", "end").strip()
        if not expr:
            raise ValueError("Digite a expressão do perfil personalizado.")
        p, avisos = pp.perfil(vp["func"], expr, vp["comp"], vp["passo"], sub, v["zref"])
        w, ee = cm.larguras(p, sub)
        ee_med = cm.eeff_medio(p.z, ee)
        info = {"tipo": "EBG perfil personalizado", "funcao": vp["func"], "expressao": expr,
                "comprimento_mm": vp["comp"], "Zref_ohm": v["zref"],
                "Z0_min_max_ohm": f"{p.z0.min():.1f}/{p.z0.max():.1f}", "eeff_medio": f"{ee_med:.4f}"}
        lay = construtores.linha_largura_modulada(p.z * 1e3, w, sub, v["lalim"], nome="ebg_personalizado",
                                                  info=info)
        return v, lay, dict(tipo="Personalizado", p=p, w=w, ee=ee, ee_med=ee_med, avisos=avisos)

    def atualizar_previa(self):
        try:
            _, lay, _ = self._gerar()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        self._lay = lay
        self.previa.mostrar(lay)

    @protegido
    def calcular(self):
        v, lay, d = self._gerar()
        self._lay = lay
        self.previa.mostrar(lay)
        if v["tipo"] == PERTURBACOES[2]:
            return self._furos(v, d)
        p, w, ee, ee_med, tipo = d["p"], d["w"], d["ee"], d["ee_med"], d["tipo"]
        f = np.linspace(v["fmin"], v["fmax"], v["npts"]) * 1e9
        s11, s21 = cm.resolver_cmt(f, p.z, p.z0, ee)
        a11, a21 = cm.resolver_abcd(f, p.z, p.z0, ee, v["zref"])
        bandas = cm.bandas_analiticas(p, ee_med)

        # S22: mesma estrutura vista da porta 2 (perfil invertido); recíproca -> S12 = S21
        s22, _ = cm.resolver_cmt(f, p.z, p.z0[::-1].copy(), ee[::-1].copy())
        s = np.empty((f.size, 2, 2), dtype=complex)
        s[:, 0, 0], s[:, 1, 0], s[:, 0, 1], s[:, 1, 1] = s11, s21, s21, s22
        self._ultimo = (f, s, v["zref"], v)
        self.estado.publicar_resultado("ebg", f, s, v["zref"])

        ax = self.grafico.eixos(2, 2)
        zmm = p.z * 1e3
        ax[0].plot(zmm, p.z0, label="Z0(z)")
        ax[0].set_ylabel("Z0 [Ω]")
        ax[0].set_xlabel("z [mm]")
        axw = ax[0].twinx()
        axw.plot(zmm, w, color="C1", label="W(z)")
        axw.set_ylabel("W [mm]", color="C1")
        ax[0].set_title("Perfil de impedância e largura", fontsize=9)
        for i, (nome, sa, sb) in enumerate((("CMT numérico", s21, s11), ("Cascata ABCD", a21, a11))):
            ls = "-" if i == 0 else "--"
            ax[1].plot(f / 1e9, abcd.db(sa), ls, label=nome)
            ax[2].plot(f / 1e9, abcd.db(sb), ls, label=nome)
        for b in bandas:
            if v["fmin"] * 1e9 <= b.f_max <= v["fmax"] * 1e9:
                fa = f[np.abs(f - b.f_max) < 1.5 * b.banda]
                a11n, a21n = cm.s_analitico(fa, b.n, b.k, p.periodo, p.comprimento, ee_med)
                ax[1].plot(fa / 1e9, abcd.db(a21n), ":", color="k", label=f"Analítico n={b.n}")
                ax[2].plot(fa / 1e9, abcd.db(a11n), ":", color="k", label=f"Analítico n={b.n}")
        personalizado = tipo == "Personalizado"
        if personalizado:
            est = pp.estimativa_fourier(f, p.z, p.z0, ee_med)
            ax[2].plot(f / 1e9, abcd.db(est), ":", color="k", label="Fourier de K(z) (eq. 2.95)")
        ax[1].set_ylabel("|S21| [dB]")
        ax[2].set_ylabel("|S11| [dB]")
        for a in ax[1:3]:
            a.set_xlabel("f [GHz]")
            a.set_ylim(-70, 3)
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        ax[0].grid(True, alpha=0.3)
        k = cm.coeficiente_acoplamento(p) / 100
        zm = 0.5 * (zmm[1:] + zmm[:-1])
        um = np.ones_like(zm, dtype=bool) if personalizado else zm <= v["periodo"] * 2
        ax[3].plot(zm[um], k[um])
        ax[3].set_xlabel("z [mm]" if personalizado else "z [mm] (2 primeiros períodos)")
        ax[3].set_ylabel("K(z) [1/cm]")
        ax[3].set_title("Coeficiente de acoplamento (eq. 2.50)", fontsize=9)
        ax[3].grid(True, alpha=0.3)
        self.grafico.desenhar()
        self.ver_resultados()

        self._banda = abcd.descrever_banda(f, s)
        if personalizado:
            return self._texto_personalizado(v, p, w, ee_med, f, s21, d["avisos"])
        linhas = [f"Perfil {tipo}: l = {v['periodo']:g} mm, {v['n']} períodos, L = {p.comprimento * 1e3:g} mm",
                  f"Z0: {p.z0.min():.1f} a {p.z0.max():.1f} Ω   W: {w.min():.3f} a {w.max():.3f} mm",
                  f"εeff médio (eq. 2.57) = {ee_med:.4f}", "",
                  "Bandas rejeitadas (analítico, eqs. 2.84-2.89):",
                  " n   fmax[GHz]  |S21|min[dB]  BW[GHz]"]
        for b in bandas:
            linhas.append(f" {b.n}   {b.f_max / 1e9:8.3f}   {20 * np.log10(b.s21_min):10.2f}   {b.banda / 1e9:7.3f}")
        linhas += ["", self._banda, "", "CMT e ABCD desprezam descontinuidades e radiação;",
                   "εeff vem de fórmula fechada (o livro usa EM)."]
        self.saida.escrever("\n".join(linhas))

    def _texto_personalizado(self, v, p, w, ee_med, f, s21, avisos):
        vp = self.form_pers.get()
        linhas = [f"Perfil personalizado ({vp['func']}), L = {vp['comp']:g} mm",
                  f"Z0: {p.z0.min():.1f} a {p.z0.max():.1f} Ω   W: {w.min():.3f} a {w.max():.3f} mm",
                  f"εeff médio (eq. 2.57) = {ee_med:.4f}", "",
                  "Bandas detectadas (|S21| < -10 dB, CMT):",
                  " f_mín[GHz]  |S21|[dB]  faixa -10 dB [GHz]"]
        bandas = pp.bandas_detectadas(f, s21)
        for a, b, c, dbm in bandas:
            linhas.append(f" {c / 1e9:9.3f}  {dbm:9.1f}   {a / 1e9:.3f} a {b / 1e9:.3f}")
        if not bandas:
            linhas.append(" nenhuma na faixa calculada")
        if w.min() < 0.1:
            linhas.append(f"\n⚠ Largura mínima {w.min():.3f} mm: verificar fabricabilidade.")
        for a in avisos:
            linhas.append("⚠ " + a)
        linhas += ["", self._banda, "", "Curva pontilhada em |S11|: estimativa de Fourier,", "válida para baixa refletividade."]
        self.saida.escrever("\n".join(linhas))

    def _furos(self, v, d):
        self._ultimo = None
        self.ver_estrutura()
        self.saida.escrever(
            f"Linha {v['zref']:g} Ω: W = {d['w']:.4f} mm, εeff = {d['ee']:.4f}\n"
            f"Furos: raio {d['r']:.3f} mm, período {v['periodo']:g} mm, {v['n']} furos\n\n"
            f"Frequência de Bragg estimada (eq. 2.9):\n  f = c/(2 l √εeff) = {d['fb'] / 1e9:.3f} GHz\n\n"
            "O Z0(z) com furos no terra não tem fórmula fechada\n"
            "(eq. 2.51 exige solução EM). Use a exportação para\n"
            "HFSS/CST para obter a resposta em frequência.")

    @protegido
    def publicar(self):
        _, lay, _ = self._gerar()
        self._lay = lay
        self.previa.mostrar(lay)
        self.estado.publicar_layout(lay)
        self.estado.msg(f"Layout '{lay.nome}' enviado para Exportar (Ctrl+E).")

    def exportar_s2p(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s, z, v = self._ultimo
            self.salvar_s2p(f"ebg_{v['n']}p", f, s, z, f"EBG {v['tipo']} l={v['periodo']}mm N={v['n']} K={v['k']}/cm (CMT)")
