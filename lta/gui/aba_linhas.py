"""Aba Linhas de transmissão convencionais (capítulo 0: base para as linhas artificiais)."""
import tkinter as tk
from tkinter import ttk

import numpy as np

from ..core import abcd
from ..core import linhas_planares as lp
from ..core.descontinuidades import extensao_ponta_aberta, texto_limites, via
from ..core.microstrip import C0, NP_DB
from ..geometria import circuitos, construtores
from .comum import AbaBase, AreaTexto, Formulario, protegido
from .previa import PreviaEstrutura

ELEMENTOS = ("Trecho de linha", "Stub em aberto", "Stub em curto", "Transformador λ/4")
MODOS = ("Síntese: Z0 → dimensões", "Análise: dimensões → Z0")


class AbaLinhas(AbaBase):
    TITULO = "Linhas de transmissão convencionais"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("tipo", "Linha", lp.TIPOS, ""),
            ("elemento", "Elemento", ELEMENTOS, ""),
            ("modo", "Modo", MODOS, ""),
            ("z0", "Z0 desejada", 50.0, "Ω"),
            ("w", "W ou S (análise)", 1.8, "mm"),
            ("g", "Fenda G / separação S (CPW, CPS)", 0.2, "mm"),
            ("b", "Distância entre planos b (stripline)", 1.6, "mm"),
            ("lterra", "Largura do terra coplanar", 3.0, "mm"),
            ("f0", "Frequência de projeto f0", 2.45, "GHz"),
            ("fase", "Comprimento elétrico em f0", 90.0, "°"),
            ("z2", "Carga do transformador Z2", 100.0, "Ω"),
            ("zp", "Z das portas", 50.0, "Ω"),
            ("rvia", "Raio das vias (curto, GCPW)", 0.2, "mm"),
            ("pvia", "Passo da cerca de vias (GCPW; 0 = sem)", 1.5, "mm"),
            ("fmax", "f máx", 7.5, "GHz"),
        ], titulo="Linha planar")
        self.form.pack(fill="x")
        self.perdas = tk.BooleanVar(value=True)
        ttk.Checkbutton(self.painel, text="Incluir perdas (dielétrico; condutor só em microstrip)",
                        variable=self.perdas).pack(anchor="w", pady=(6, 0))
        self.corrigir = tk.BooleanVar(value=True)
        ttk.Checkbutton(self.painel, text="Descontar ponta aberta / via do stub (Kirschning-Jansen, Goldfarb-Pucel)",
                        variable=self.corrigir, command=self.atualizar_previa).pack(anchor="w")
        self.botoes(("Calcular", self.calcular), ("Enviar p/ exportação", self.publicar),
                    ("Exportar S2P", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 20)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Linhas planares quase-TEM (livro, cap. 1):\n"
            "• Microstrip: Hammerstad-Jensen (+ dispersão/perdas)\n"
            "• Stripline centrada: Cohn (exata, t = 0)\n"
            "• CPW, CPW com terra e CPS: Simons (integrais\n"
            "  elípticas, substrato de espessura h)\n\n"
            "Elementos: trecho, stub em aberto/curto em derivação\n"
            "e transformador λ/4 entre Z0 e Z2 (√(Z0·Z2)).")
        self._ultimo = None
        self.observar_previa(self.form)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    # ------------------------------------------------------------------ projeto
    def _geometria(self, v, sub, z):
        tipo = v["tipo"]
        f0 = v["f0"] * 1e9
        if v["modo"] == MODOS[0]:
            w = lp.sintetizar(tipo, z, sub, g=v["g"], b=v["b"], f=f0)
        else:
            w = v["w"]
        geo = lp.Geometria(tipo, w, v["g"], v["b"])
        z0, ee = lp.analisar(geo, sub, f0)
        return geo, z0, ee

    def _gerar(self):
        v = self.form.get()
        sub = self.estado.substrato
        if v["f0"] <= 0 or v["fase"] <= 0:
            raise ValueError("f0 e comprimento elétrico devem ser positivos.")
        if v["elemento"] != ELEMENTOS[0] and v["tipo"] != lp.TIPOS[0]:
            raise ValueError("Stubs e transformador estão implementados só em microstrip.")
        f0 = v["f0"] * 1e9
        geo, z0, ee = self._geometria(v, sub, v["z0"])
        extra = {}
        if v["elemento"] == ELEMENTOS[3]:
            zt = np.sqrt(v["z0"] * v["z2"])
            geo_t, zt_r, ee_t = self._geometria(dict(v, modo=MODOS[0]), sub, zt)
            geo_2, z2_r, _ = self._geometria(dict(v, modo=MODOS[0]), sub, v["z2"])
            comp = lp.comprimento_eletrico(geo_t, sub, f0, 90.0)
            extra = dict(geo_t=geo_t, zt=zt_r, ee_t=ee_t, geo_2=geo_2, z2=z2_r)
        else:
            comp = lp.comprimento_eletrico(geo, sub, f0, v["fase"])
        # ponta aberta: o campo de franja alonga o stub de Δl; o layout desconta Δl se marcado.
        # stub em curto: a via (Goldfarb-Pucel) também alonga o stub; desconta-se o equivalente em linha.
        extra["dl"] = extensao_ponta_aberta(geo.w, sub) if v["elemento"] == ELEMENTOS[1] else 0.0
        extra["dl_via"] = (lp.encurtamento_via(v["rvia"], geo, sub, f0)
                           if v["elemento"] == ELEMENTOS[2] and v["rvia"] > 0 else 0.0)
        extra["comp_fisico"] = comp - ((extra["dl"] + extra["dl_via"]) if self.corrigir.get() else 0.0)
        if extra["comp_fisico"] <= 0:
            raise ValueError("Stub mais curto que a correção de ponta aberta.")
        info = {"tipo": f"{v['tipo']} - {v['elemento']}", "Z0_ohm": f"{z0:.2f}", "eeff": f"{ee:.4f}",
                "W_mm": f"{geo.w:.4f}", "comprimento_mm": f"{comp:.4f}", "f0_GHz": v["f0"]}
        lay = None
        if v["tipo"] == lp.TIPOS[5]:
            lay = circuitos.slotline("trecho_slotline", sub, geo.w, comp, max(v["lterra"], 5 * geo.w), 5.0,
                                     info=info)
        elif v["tipo"] != lp.TIPOS[1]:
            d = dict(w=geo.w, g=v["g"], comp=comp, largura_terra=v["lterra"], raio_via=v["rvia"],
                     passo_vias=v["pvia"], w_stub=geo.w, comp_stub=extra["comp_fisico"])
            if v["elemento"] == ELEMENTOS[3]:
                d.update(w1=geo.w, wt=extra["geo_t"].w, w2=extra["geo_2"].w)
            nome = {ELEMENTOS[0]: "trecho", ELEMENTOS[1]: "stub_aberto", ELEMENTOS[2]: "stub_curto",
                    ELEMENTOS[3]: "transformador"}[v["elemento"]]
            lay = construtores.linha_convencional(v["tipo"], v["elemento"], d, sub, 5.0,
                                                  nome=f"{nome}_{v['tipo'].split()[0].lower()}", info=info)
        return v, sub, geo, z0, ee, comp, extra, lay

    def atualizar_previa(self):
        try:
            *_, lay = self._gerar()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        if lay is None:
            self.previa.mostrar_erro("Stripline (três camadas): só análise; layout não gerado.")
        else:
            self.previa.mostrar(lay)

    @protegido
    def calcular(self):
        v, sub, geo, z0, ee, comp, extra, lay = self._gerar()
        if lay is not None:
            self.previa.mostrar(lay)
        f = np.linspace(v["fmax"] * 1e9 / 1500, v["fmax"] * 1e9, 1500)
        perdas = self.perdas.get()
        el = v["elemento"]
        z1p, z2p = v["zp"], v["zp"]
        if el == ELEMENTOS[0]:
            m = lp.abcd_trecho(f, geo, sub, comp, perdas)
        elif el in ELEMENTOS[1:3]:
            m = lp.abcd_stub(f, geo, sub, extra["comp_fisico"] + extra["dl"], curto=el == ELEMENTOS[2],
                             perdas=perdas, raio_via=v["rvia"] if v["rvia"] > 0 else None)
        else:
            m = lp.abcd_trecho(f, extra["geo_t"], sub, comp, perdas)
            z1p, z2p = v["z0"], v["z2"]
        s = abcd.abcd_para_s2(m, z1p, z2p)
        self._ultimo = (f, s, z1p, v)
        self.estado.publicar_resultado("linha_convencional", f, s, z1p)
        fg = f / 1e9
        ax = self.grafico.eixos(1, 2)
        ax[0].plot(fg, abcd.db(s[:, 1, 0]), label="|S21|")
        ax[0].plot(fg, abcd.db(s[:, 0, 0]), label="|S11|")
        ax[0].axvline(v["f0"], color="gray", ls=":", lw=0.8)
        ax[0].set_ylim(-50, 2)
        ax[0].set_ylabel("|S| [dB]")
        ax[0].set_title(f"{el}, portas {z1p:g}/{z2p:g} Ω", fontsize=9)
        ax[1].plot(fg, np.degrees(np.angle(s[:, 1, 0])), label="∠S21")
        ax[1].set_ylabel("Fase [°]")
        ax[1].set_title("Fase de transmissão", fontsize=9)
        for a in ax:
            a.set_xlabel("f [GHz]")
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        f0 = v["f0"] * 1e9
        g0, _ = lp.gama(np.array([f0]), geo, sub, True)
        lg = C0 / (f0 * np.sqrt(ee)) * 1e3
        nomes = {lp.TIPOS[0]: "W", lp.TIPOS[1]: "W", lp.TIPOS[4]: "W (fitas)", lp.TIPOS[5]: "fenda"}
        dims = f"{nomes.get(geo.tipo, 'S')} = {geo.w:.4f} mm"
        if geo.tipo in lp.TIPOS[2:5]:
            dims += f", {'S' if geo.tipo == lp.TIPOS[4] else 'G'} = {geo.g:g} mm"
        if geo.tipo == lp.TIPOS[1]:
            dims += f", b = {geo.b:g} mm"
        linhas = [f"{geo.tipo} em {sub.nome} (εr = {sub.er:g}, h = {sub.h_mm:g} mm)",
                  f"{dims}  →  Z0 = {z0:.2f} Ω, εeff = {ee:.4f}",
                  f"λg em {v['f0']:g} GHz = {lg:.3f} mm;  vp = {1 / np.sqrt(ee):.3f}·c",
                  f"α em f0 = {g0.real[0] * NP_DB / 100:.4f} dB/cm"
                  + (" (dielétrico + condutor)" if geo.tipo in lp.TIPOS[:1] + lp.TIPOS[2:5]
                     else " (só dielétrico)")]
        if geo.tipo == lp.TIPOS[5]:
            linhas.append("Slotline dispersiva: Z0 e εeff em f0 (Janaswamy-Schaubert / Garg-Gupta)")
        if sub.incluir_t:
            linhas.append(f"Z0 e εeff incluem a espessura t = {sub.t_um:g} µm (Garg-Bahl)")
        if el == ELEMENTOS[3]:
            linhas += ["", f"Transformador λ/4: Zt = √(Z0·Z2) = {extra['zt']:.2f} Ω, W = {extra['geo_t'].w:.4f} mm, "
                       f"εeff = {extra['ee_t']:.4f}", f"Comprimento λ/4 em f0: {comp:.3f} mm",
                       f"Z2 = {v['z2']:g} Ω → W = {extra['geo_2'].w:.4f} mm"]
            ok = abcd.db(s[:, 0, 0]) < -20
            i = np.argmin(abs(f - f0))
            if ok[i]:
                a = b = i
                while a > 0 and ok[a - 1]:
                    a -= 1
                while b < len(f) - 1 and ok[b + 1]:
                    b += 1
                linhas.append(f"|S11| < -20 dB de {f[a] / 1e9:.3f} a {f[b] / 1e9:.3f} GHz "
                              f"({100 * (f[b] - f[a]) / f0:.1f}%)")
        else:
            linhas += ["", f"{el}: {v['fase']:g}° em f0 → {comp:.3f} mm"]
            if el != ELEMENTOS[0]:
                tipo_ress = "rejeição" if (el == ELEMENTOS[1]) == (abs(v["fase"] - 90) < 1) else "passagem"
                linhas.append(f"Stub de {v['fase']:g}° em f0: {tipo_ress} em f0 (zero de transmissão "
                              "nos múltiplos ímpares de 90° para o aberto e pares para o curto)")
            if el == ELEMENTOS[2] and extra["dl_via"]:
                lv, rv = via(v["rvia"], sub, f0)
                linhas.append(f"Via (Goldfarb-Pucel): L = {lv * 1e12:.1f} pH, R = {float(rv) * 1e3:.2f} mΩ em f0; "
                              f"equivale a {extra['dl_via']:.3f} mm de linha"
                              + (" (descontado)" if self.corrigir.get() else ""))
            if el == ELEMENTOS[1]:
                linhas.append(f"Ponta aberta: Δl = {extra['dl']:.3f} mm (Kirschning-Jansen); stub físico = "
                              f"{extra['comp_fisico']:.3f} mm" + ("" if self.corrigir.get() else
                                                                   " (sem desconto: ressonância abaixo de f0)"))
        if lay is None:
            linhas.append("\nStripline: layout em três camadas não gerado.")
        elif lay.sem_terra:
            linhas.append("\nSem plano de terra inferior: a exportação omite a camada TERRA.")
        if geo.tipo == lp.TIPOS[0]:
            linhas += ["", "Limites do substrato:"] + texto_limites(sub, geo.w, v["fmax"] * 1e9)
        linhas += ["", "Modelos quase-estáticos com metal de espessura nula;",
                   "stubs sem efeito de junção T e via (ponta aberta incluída)."]
        self.saida.escrever("\n".join(linhas))

    @protegido
    def publicar(self):
        *_, lay = self._gerar()
        if lay is None:
            raise ValueError("Stripline: layout não disponível.")
        self.previa.mostrar(lay)
        self.estado.publicar_layout(lay)
        self.estado.msg(f"Layout '{lay.nome}' enviado para Exportar (Ctrl+E).")

    def exportar_s2p(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s, z, v = self._ultimo
            self.salvar_s2p("linha_convencional", f, s, z, f"{v['tipo']} - {v['elemento']}")
