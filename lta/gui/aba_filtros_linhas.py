"""Aba Filtros com linhas convencionais (Pozar, cap. 8), com passa-faixa de linhas acopladas LTA e
combline com capacitores (Steer, seção 10.14)."""
import tkinter as tk
from tkinter import ttk

import numpy as np

from ..core import abcd
from ..core import combline as cl
from ..core import filtros
from ..core import filtros_convencionais as fc
from ..core import linhas_acopladas as la
from ..core import linhas_planares as lp
from ..core.descontinuidades import extensao_ponta_aberta, gap, gap_para_capacitancia, texto_limites
from ..core.microstrip import C0, analise, sintese
from ..geometria import circuitos
from .comum import AbaBase, AreaTexto, Formulario, protegido
from .previa import PreviaEstrutura

TIPOS = ("Passa-baixas: degraus de impedância", "Passa-baixas: stubs (Richards + Kuroda)",
         "Passa-faixa: linhas acopladas em paralelo", "Rejeita-faixa: stubs λ/4 abertos",
         "Passa-faixa: stubs λ/4 em curto", "Passa-faixa: ressonadores série com gaps",
         "Passa-faixa: stubs em curto acoplados por capacitores",
         "Passa-faixa: combline com capacitores (Steer 10.14)")
APROX = ("Chebyshev", "Butterworth")
IMPL = ("Convencional", "LTA: linhas acopladas com modo ímpar compensado")


class AbaFiltrosLinhas(AbaBase):
    TITULO = "Filtros com linhas (Pozar, cap. 8)"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("tipo", "Filtro", TIPOS, ""),
            ("aprox", "Resposta", APROX, ""),
            ("n", "Ordem N", 3, ""),
            ("ripple", "Ripple (Chebyshev)", 0.5, "dB"),
            ("f0", "fc (passa-baixas) ou f0", 2.0, "GHz"),
            ("fbw", "Banda fracionária Δ", 10.0, "%"),
            ("z0", "Z0", 50.0, "Ω"),
            ("zl", "Degraus: Zℓ (baixa)", 20.0, "Ω"),
            ("zh", "Degraus: Zh (alta)", 120.0, "Ω"),
            ("impl", "Linhas acopladas", IMPL, ""),
            ("ncel", "LTA: capacitores por seção", 6, ""),
            ("modelo", "Modelo das linhas acopladas", la.MODELOS, ""),
            ("rvia", "Raio das vias (stubs em curto)", 0.2, "mm"),
            ("zres", "Combline: Z do ressonador central", 80.0, "Ω"),
            ("fmax", "f máx", 7.0, "GHz"),
        ], titulo="Filtro")
        self.form.pack(fill="x")
        self.corrigir = tk.BooleanVar(value=True)
        self.perdas = tk.BooleanVar(value=False)
        self.dispersao = tk.BooleanVar(value=False)
        ttk.Checkbutton(self.painel, text="Corrigir ponta aberta, via e degraus (Garg-Bahl)",
                        variable=self.corrigir).pack(anchor="w", pady=(6, 0))
        ttk.Checkbutton(self.painel, text="Perdas na microstrip", variable=self.perdas).pack(anchor="w")
        ttk.Checkbutton(self.painel, text="Dispersão das linhas acopladas (Kirschning-Jansen)",
                        variable=self.dispersao).pack(anchor="w")
        self.botoes(("Calcular", self.calcular), ("Enviar p/ exportação", self.publicar),
                    ("Ex. 8.6: degraus", self.ex_86), ("Ex. 8.5: Kuroda", self.ex_85),
                    ("Ex. 8.7: linhas acopladas", self.ex_87), ("LTA: sem espúria em 2f0", self.ex_lta),
                    ("Steer: combline 1 GHz", self.ex_combline),
                    ("Exportar S2P", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 22)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Protótipos de Pozar (Tabelas 8.3-8.4) → linhas.\n\n"
            "Curvas: rede ideal (tracejado) e microstrip no substrato atual\n"
            "(contínuo). No passa-faixa de linhas acopladas em microstrip,\n"
            "εeff par ≠ ímpar cria a banda espúria em 2f0; a opção LTA\n"
            "carrega o modo ímpar com capacitores e iguala as velocidades\n"
            "(Martín, seção 4.2.3.2: supressão de espúrias).")
        self._ultimo = None
        self.observar_previa(self.form)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    # ------------------------------------------------------------------ exemplos
    def ex_86(self):
        from ..core.microstrip import Substrato
        self.estado.definir_substrato(Substrato(4.2, 1.58, 12.7, 0.02, "Pozar Ex. 8.6 (er = 4,2)"))
        self.form.set(dict(tipo=TIPOS[0], aprox=APROX[1], n=6, f0=2.5, z0=50, zl=20, zh=120, fmax=8))
        self.calcular()

    def ex_85(self):
        self.form.set(dict(tipo=TIPOS[1], aprox=APROX[0], n=3, ripple=3.0, f0=4, z0=50, fmax=16))
        self.calcular()

    def ex_87(self):
        self.form.set(dict(tipo=TIPOS[2], aprox=APROX[0], n=3, ripple=0.5, f0=2, fbw=10, z0=50, impl=IMPL[0],
                           fmax=7))
        self.calcular()

    def ex_lta(self):
        self.form.set(dict(tipo=TIPOS[2], aprox=APROX[0], n=3, ripple=0.5, f0=2, fbw=10, z0=50, impl=IMPL[1],
                           ncel=6, fmax=7))
        self.calcular()

    # ------------------------------------------------------------------ projeto
    def ex_combline(self):
        r = 10 * np.log10(1 + 0.1**2)
        self.form.set(dict(tipo=TIPOS[7], aprox=APROX[0], n=3, ripple=f"{r:.5f}", f0=1, fbw=10, z0=50, zres=80,
                           fmax=4))
        self.calcular()

    def _projeto(self):
        v = self.form.get()
        sub = self.estado.substrato
        if v["f0"] <= 0 or v["n"] < 1 or v["z0"] <= 0:
            raise ValueError("f0, N e Z0 devem ser positivos.")
        g = filtros.prototipo(v["aprox"], v["n"], v["ripple"])
        f0, z0, fbw = v["f0"] * 1e9, v["z0"], v["fbw"] / 100
        tipo = v["tipo"]
        w0 = float(sintese(z0, sub))
        p = dict(v=v, sub=sub, g=g, f0=f0, z0=z0, fbw=fbw, tipo=tipo, w0=w0, lay=None, texto=[])
        p["texto"].append("g = " + ", ".join(f"{x:.4f}" for x in g))
        corr = self.corrigir.get()

        def largura(z):
            try:
                return float(lp.sintetizar("Microstrip", z, sub))
            except ValueError as e:
                raise ValueError(f"{e} Escolha impedâncias realizáveis no substrato atual.") from None

        def lg_de(w):
            return C0 / (f0 * np.sqrt(float(analise(w, sub)[1]))) * 1e3

        info = {"tipo": f"filtro - {tipo}", "N": v["n"], "f0_GHz": v["f0"], "resposta": v["aprox"]}
        if tipo == TIPOS[0]:
            tr = fc.degraus(g, z0, v["zl"], v["zh"], comeca_com_c=True)
            segs = []
            p["texto"].append(f"Degraus (eq. 8.86), fc = {v['f0']:g} GHz, Zℓ = {v['zl']:g} Ω, Zh = {v['zh']:g} Ω:")
            for i, t in enumerate(tr, 1):
                w = largura(t.z)
                comp = t.bl / (2 * np.pi) * lg_de(w)
                segs.append((comp, w))
                p["texto"].append(f"  {i}: {t.elemento}  Z = {t.z:g} Ω  βℓ = {np.degrees(t.bl):.1f}°  W = {w:.3f} mm"
                                  f"  ℓ = {comp:.3f} mm")
            p["segs_projeto"] = list(segs)
            if corr and len(segs) > 1:
                segs, dls = fc.correcao_degraus(segs, sub, w0, w0)
                p["texto"].append("  degraus (seção 3.4.3): Δl = " + ", ".join(f"{d:.3f}" for d in dls)
                                  + " mm; trechos largos encurtados e estreitos alongados:")
                p["texto"].append("  ℓ corrigidos = " + ", ".join(f"{c:.3f}" for c, _ in segs) + " mm")
            p["tr"], p["segs"] = tr, segs
            x, ws = [], []
            pos = 0.0
            for comp, w in segs:
                x += [pos, pos + comp]
                ws += [w, w]
                pos += comp
            p["lay"] = circuitos.linha_perfil("filtro_degraus", sub, x, ws, w0, w0, 5.0, info=info)
        elif tipo == TIPOS[1]:
            els = fc.richards_kuroda(g)
            p["els"] = els
            p["texto"].append(f"Richards + Kuroda (Ex. 8.5): todos os elementos λ/8 em fc = {v['f0']:g} GHz")
            segs, stubs, pos = [], [], 0.0
            for e in els:
                z = e.z * z0
                w = largura(z)
                l8 = lg_de(w) / 8
                nome = {"UE": "linha (UE)", "P": "stub aberto em derivação"}[e.tipo]
                p["texto"].append(f"  {nome}: Z = {z:.1f} Ω  W = {w:.3f} mm  λ/8 = {l8:.3f} mm")
                if e.tipo == "UE":
                    segs.append((l8, w))
                    pos += l8
                else:
                    dl = extensao_ponta_aberta(w, sub) if corr else 0.0
                    stubs.append(dict(x=pos, w=w, comp=l8 - dl, curto=False, z=z))
            if segs and stubs[0]["x"] == 0.0:
                segs.insert(0, (stubs[0]["w"] / 2 + 0.2, w0))
                for st in stubs:
                    st["x"] += segs[0][0]
            if stubs[-1]["x"] >= sum(c for c, _ in segs) - 1e-9:
                segs.append((stubs[-1]["w"] / 2 + 0.2, w0))
            p["segs"], p["stubs"] = segs, stubs
            p["lay"] = circuitos.linha_com_stubs("filtro_stubs_kuroda", sub, segs, stubs, w0, w0, 5.0, info=info)
        elif tipo == TIPOS[2]:
            lta = v["impl"] == IMPL[1]
            pa = fc.projeto_linhas_acopladas(g, f0, fbw, z0, sub, la.TIPOS[0], lta=lta, n_cel=max(1, v["ncel"]),
                                             corrigir_ponta=corr, modelo=v["modelo"],
                                             dispersivo=self.dispersao.get())
            p["pa"] = pa
            p["texto"] += [f"Linhas acopladas, Δ = {100 * fbw:g}% (eqs. 8.108 e 8.121)"
                           + (" – LTA compensado" if lta else "") + ":"] + pa.texto
            info["implementacao"] = v["impl"]
            for i, s in enumerate(pa.secoes, 1):
                info[f"secao{i}_W_S_l_mm"] = f"{s.w:.4f} / {s.s:.4f} / {s.comp:.3f}"
                if s.n_cel:
                    info[f"secao{i}_Cc_fF"] = f"{s.cc * 1e15:.2f} x {s.n_cel}"
            p["lay"] = circuitos.filtro_linhas_acopladas("filtro_acopladas" + ("_lta" if lta else ""), sub,
                                                         [(s.comp, s.w, s.s) for s in pa.secoes], w0, 5.0, info)
        elif tipo in TIPOS[3:5]:
            rej = tipo == TIPOS[3]
            zs = fc.stubs_quarto(g, fbw, z0, rejeita=rej)
            p["zs"] = zs
            p["texto"].append(f"Stubs λ/4 {'abertos (eq. 8.130)' if rej else 'em curto (eq. 8.131)'}, linhas de "
                              f"ligação λ/4 de Z0:")
            l4 = lg_de(w0) / 4
            segs, stubs = [], []
            erro = None
            for i, z in enumerate(zs):
                p["texto"].append(f"  stub {i + 1}: Z = {z:.1f} Ω")
                try:
                    w = largura(z)
                except ValueError as e:
                    erro = str(e)
                    continue
                if not corr:
                    dl = 0.0
                elif rej:
                    dl = extensao_ponta_aberta(w, sub)
                else:
                    dl = lp.encurtamento_via(v["rvia"], lp.Geometria("Microstrip", w), sub, f0)
                stubs.append(dict(x=i * l4 + w0, w=w, comp=lg_de(w) / 4 - dl, curto=not rej))
            if erro:
                p["texto"].append(f"⚠ {erro} (sem layout; resposta só ideal)")
            else:
                segs = [(w0, w0)] + [(l4, w0)] * (len(zs) - 1) + [(w0, w0)]
                p["segs"], p["stubs"] = segs, stubs
                p["lay"] = circuitos.linha_com_stubs("filtro_stubs_quarto", sub, segs, stubs, w0, w0, 5.0,
                                                 raio_via=v["rvia"], info=info)
        elif tipo == TIPOS[5]:
            pg = fc.ressonadores_serie_gap(g, f0, fbw, z0)
            p["pg"] = pg
            p["texto"].append("Ressonadores série acoplados por gaps (eqs. 8.134-8.135):")
            lg0 = lg_de(w0)
            gaps, cps, erro = [], [], None
            for i, (j, b, c) in enumerate(zip(pg.jz, pg.b, pg.c), 1):
                txt = f"  gap {i}: Z0·J = {j:.4f}  B = {b:.3e} S  C = {c * 1e12:.4f} pF"
                try:
                    s = gap_para_capacitancia(w0, c, sub)
                    gaps.append(s)
                    cps.append(gap(w0, s, sub)[1])
                    txt += f"  s = {s * 1e3:.1f} µm"
                except ValueError as e:
                    erro = str(e)
                p["texto"].append(txt)
            comps = []
            for i, th in enumerate(pg.theta, 1):
                comp = th / (2 * np.pi) * lg0
                if cps and len(cps) == len(pg.c) and corr:      # Cp dos gaps encurta o ressonador
                    comp -= (cps[i - 1] + cps[i]) * z0 * C0 / np.sqrt(float(analise(w0, sub)[1])) * 1e3
                comps.append(comp)
                p["texto"].append(f"  ressonador {i}: θ = {np.degrees(th):.1f}°  ℓ = {comp:.3f} mm")
            if erro:
                p["texto"].append(f"⚠ {erro} Sem layout (use capacitores SMD ou gaps interdigitais).")
            else:
                p["gaps"], p["cps"], p["comps"] = gaps, cps, comps
                p["lay"] = circuitos.filtro_gaps("filtro_gaps", sub, w0, gaps, comps, 5.0, info=info)
        elif tipo == TIPOS[7]:
            pc = cl.projeto(f0, fbw, v["n"], v["aprox"], v["ripple"], z0, v["zres"])
            p["pc"] = pc
            p["texto"] = list(pc.texto)
            pares = cl.pares_acoplados(pc)
            ws, ss, comps = [], [], []
            p["texto"].append("Pares de linhas vizinhas (estimativa por pares, meio homogêneo):")
            try:
                for i, (ze, zo) in enumerate(pares, 1):
                    w, s = la.sintetizar(la.TIPOS[0], ze, zo, sub, modelo=v["modelo"])
                    md = la.analisar(la.TIPOS[0], w, s, sub, modelo=v["modelo"])
                    comps.append(la.comprimento_quarto(md, pc.fr))
                    ws.append(w)
                    ss.append(s)
                    p["texto"].append(f"  {i}-{i + 1}: Z0e = {ze:.2f} Ω, Z0o = {zo:.2f} Ω → W = {w:.4f} mm, "
                                      f"S = {s:.4f} mm")
            except ValueError as e:
                p["texto"].append(f"⚠ {e} Sem layout.")
                ws = []
            if ws:
                larg = [ws[0]] + [0.5 * (a + b) for a, b in zip(ws[:-1], ws[1:])] + [ws[-1]]
                comp = float(np.mean(comps))
                p["texto"].append(f"Linhas de {comp:.3f} mm (λ/4 em fr = {pc.fr / 1e9:g} GHz, média dos modos); "
                                  "microstrip não é homogênea: ajustar W, S e C em EM.")
                info.update(N=v["n"], fr_GHz=pc.fr / 1e9)
                p["lay"] = circuitos.combline("filtro_combline", sub, larg, ss, comp, w0, pc.c, pc.cb, v["rvia"],
                                              info=info)
        else:
            pd = fc.ressonadores_derivacao(g, f0, fbw, z0)
            p["pd"] = pd
            p["texto"].append("Stubs em curto (Z0) acoplados por capacitores série (eqs. 8.136-8.141):")
            for i, (j, c) in enumerate(zip(pd.jz, pd.c)):
                p["texto"].append(f"  C{i}{i + 1}: Z0·J = {j:.4f}  C = {c * 1e12:.4f} pF")
            for i, (dc, l) in enumerate(zip(pd.dc, pd.l), 1):
                p["texto"].append(f"  stub {i}: ΔC = {dc * 1e12:.4f} pF  ℓ = {l:.4f}λ ({360 * l:.1f}°) = "
                                  f"{l * lg_de(w0):.3f} mm")
            p["texto"].append("Capacitores concentrados (SMD): sem layout de linhas.")
        return p

    def _redes(self, p, f):
        """(S ideal, S microstrip ou None)."""
        tipo, f0, z0, g, sub = p["tipo"], p["f0"], p["z0"], p["g"], p["sub"]
        perdas = self.perdas.get()
        s_ms = None
        if tipo == TIPOS[0]:
            s_id = abcd.abcd_para_s(fc.abcd_degraus_ideal(f, f0, p["tr"]), z0)
            s_ms = abcd.abcd_para_s(fc.abcd_degraus_microstrip(f, p["segs"], sub, perdas, True, p["w0"], p["w0"]),
                                    z0)
            p["s_conc"] = abcd.abcd_para_s(fc.abcd_concentrado(f, f0, g, z0), z0)
        elif tipo == TIPOS[1]:
            s_id = abcd.abcd_para_s(fc.abcd_elementos(f, f0, p["els"], z0), z0)
            s_ms = self._stubs_microstrip(p, f)
            p["s_conc"] = abcd.abcd_para_s(fc.abcd_concentrado(f, f0, g, z0, comeca_com_c=False), z0)
        elif tipo == TIPOS[2]:
            pa = p["pa"]
            s_id = fc.s_linhas_acopladas(f, fc.secoes_ideais(f0, pa.ze, pa.zo), z0)
            s_ms = fc.s_linhas_acopladas(f, pa.secoes, z0, sub, perdas)
        elif tipo in TIPOS[3:5]:
            s_id = abcd.abcd_para_s(fc.abcd_stubs_quarto(f, f0, p["zs"], z0, rejeita=tipo == TIPOS[3]), z0)
            if "stubs" in p:
                s_ms = self._stubs_microstrip(p, f)
        elif tipo == TIPOS[5]:
            s_id = abcd.abcd_para_s(fc.abcd_gaps_ideal(f, f0, p["pg"], z0), z0)
            if "gaps" in p:
                geo = lp.Geometria("Microstrip", p["w0"])
                w = 2 * np.pi * f
                mats = []
                for i, s in enumerate(p["gaps"]):
                    cg, cp = gap(p["w0"], s, sub)
                    mats += [abcd.shunt(1j * w * cp), abcd.series(1 / (1j * w * cg)), abcd.shunt(1j * w * cp)]
                    if i < len(p["comps"]):
                        mats.append(lp.abcd_trecho(f, geo, sub, p["comps"][i], perdas))
                s_ms = abcd.abcd_para_s(abcd.cascata(*mats), z0)
        elif tipo == TIPOS[7]:
            s_id = cl.s_combline(f, p["pc"])
        else:
            s_id = abcd.abcd_para_s(fc.abcd_derivacao_ideal(f, f0, p["pd"], z0), z0)
        return s_id, s_ms

    def _stubs_microstrip(self, p, f):
        """Linha principal (segmentos) com os stubs nas posições do layout, em microstrip real."""
        sub, perdas = p["sub"], self.perdas.get()
        stubs = sorted(p["stubs"], key=lambda st: st["x"])
        mats, pos, x, k = [], 0.0, 0.0, 0
        for comp, w in p["segs"]:
            geo = lp.Geometria("Microstrip", w)
            fim = x + comp
            while k < len(stubs) and stubs[k]["x"] <= fim + 1e-9:
                st = stubs[k]
                if st["x"] > pos:
                    mats.append(lp.abcd_trecho(f, geo, sub, st["x"] - pos, perdas))
                    pos = st["x"]
                comp_st = st["comp"] + (0.0 if st["curto"] else extensao_ponta_aberta(st["w"], sub))
                mats.append(lp.abcd_stub(f, lp.Geometria("Microstrip", st["w"]), sub, comp_st, curto=st["curto"],
                                         perdas=perdas, raio_via=p["v"]["rvia"] if st["curto"] else None))
                k += 1
            if fim > pos:
                mats.append(lp.abcd_trecho(f, geo, sub, fim - pos, perdas))
                pos = fim
            x = fim
        return abcd.abcd_para_s(abcd.cascata(*mats), p["z0"])

    def atualizar_previa(self):
        try:
            p = self._projeto()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        if p["lay"] is None:
            self.previa.mostrar_erro("Sem layout para este projeto (elementos concentrados ou impedâncias "
                                     "não realizáveis).")
        else:
            self.previa.mostrar(p["lay"])

    @protegido
    def calcular(self):
        p = self._projeto()
        v = p["v"]
        if p["lay"] is not None:
            self.previa.mostrar(p["lay"])
        f = np.linspace(v["fmax"] * 1e9 / 1400, v["fmax"] * 1e9, 1400)
        s_id, s_ms = self._redes(p, f)
        s_out = s_ms if s_ms is not None else s_id
        self._ultimo = (f, s_out, p["z0"], p["tipo"])
        self.estado.publicar_resultado("filtro_linhas", f, s_out, p["z0"])
        fg = f / 1e9
        ax = self.grafico.eixos(1, 2)
        if "s_conc" in p:
            ax[0].plot(fg, abcd.db(p["s_conc"][:, 1, 0]), color="#999", lw=0.8, label="|S21| concentrado")
        ax[0].plot(fg, abcd.db(s_id[:, 1, 0]), ls="--", label="|S21| ideal")
        if s_ms is not None:
            ax[0].plot(fg, abcd.db(s_ms[:, 1, 0]), label="|S21| microstrip")
            ax[0].plot(fg, abcd.db(s_ms[:, 0, 0]), lw=0.8, label="|S11| microstrip")
        else:
            ax[0].plot(fg, abcd.db(s_id[:, 0, 0]), lw=0.8, label="|S11| ideal")
        ax[0].axvline(v["f0"], color="gray", ls=":", lw=0.8)
        ax[0].set_ylim(-60, 2)
        ax[0].set_ylabel("|S| [dB]")
        ax[0].set_title(p["tipo"], fontsize=9)
        z = (abs(fg - v["f0"]) < 0.4 * v["f0"]) if p["tipo"] not in TIPOS[:2] else fg < 1.3 * v["f0"]
        ax[1].plot(fg[z], abcd.db(s_id[z, 1, 0]), ls="--", label="ideal")
        if s_ms is not None:
            ax[1].plot(fg[z], abcd.db(s_ms[z, 1, 0]), label="microstrip")
        ax[1].set_ylim(-5 if p["tipo"] != TIPOS[3] else -60, 0.5)
        ax[1].set_title("Detalhe da banda", fontsize=9)
        for a in ax:
            a.set_xlabel("f [GHz]")
            a.grid(True, alpha=0.3)
            a.legend(fontsize=7)
        self.grafico.desenhar()
        linhas = [p["tipo"], f"{v['aprox']} N = {v['n']}" + (f", ripple {v['ripple']:g} dB" if v["aprox"] ==
                                                            "Chebyshev" else ""), ""] + p["texto"]
        if p["tipo"] == TIPOS[2]:
            s_ref = s_ms if s_ms is not None else s_id
            i2 = int(np.argmin(np.abs(f - 2 * p["f0"])))
            janela = np.abs(f - 2 * p["f0"]) < 0.15 * p["f0"]
            linhas += ["", f"Espúria em torno de 2f0: |S21| máx = {abcd.db(s_ref[janela, 1, 0]).max():.1f} dB "
                           f"(ideal TEM: {abcd.db(s_id[i2, 1, 0]):.1f} dB em 2f0)"]
        linhas += ["", "Limites do substrato:"] + texto_limites(p["sub"], p["w0"], v["fmax"] * 1e9)
        linhas += ["", "Junções T, degraus e acoplamento entre stubs não modelados: ajustar em EM."]
        self.saida.escrever("\n".join(linhas))

    @protegido
    def publicar(self):
        p = self._projeto()
        if p["lay"] is None:
            raise ValueError("Este projeto não gera layout de linhas.")
        self.previa.mostrar(p["lay"])
        self.estado.publicar_layout(p["lay"])
        self.estado.msg(f"Layout '{p['lay'].nome}' enviado para Exportar (Ctrl+E).")

    def exportar_s2p(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s, z0, tipo = self._ultimo
            self.salvar_s2p("filtro_linhas", f, s, z0, tipo)
