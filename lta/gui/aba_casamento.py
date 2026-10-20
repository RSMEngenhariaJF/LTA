"""Aba Casamento de impedância (Pozar, cap. 5): rede L, stubs, transformadores multisseção e tapers;
redes Pi e T com Q especificado (Steer, seção 8.6)."""
import tkinter as tk
from tkinter import ttk

import numpy as np

from ..core import abcd
from ..core import casamento as cs
from ..core import filtros_convencionais as fc
from ..core import linhas_planares as lp
from ..core.descontinuidades import extensao_ponta_aberta, texto_limites
from ..core.microstrip import C0, analise, parametros_linha, sintese
from ..geometria import circuitos
from .comum import AbaBase, AreaTexto, Formulario, desenhar_smith, protegido
from .previa import PreviaEstrutura

METODOS = ("Stub simples em derivação", "Stub simples em série", "Stub duplo em derivação",
           "Rede L (elementos concentrados)", "Transformador binomial", "Transformador Chebyshev",
           "Taper exponencial", "Taper triangular", "Taper Klopfenstein",
           "Rede Pi com Q especificado (Steer 8.6)", "Rede T com Q especificado (Steer 8.6)")
PONTAS = ("Aberto", "Curto")


class AbaCasamento(AbaBase):
    TITULO = "Casamento de impedância (Pozar, cap. 5)"
    PREVIA = PreviaEstrutura

    def __init__(self, mestre, estado):
        super().__init__(mestre, estado)
        self.form = Formulario(self.painel, [
            ("metodo", "Método", METODOS, ""),
            ("rl", "Carga: R em f0", 60.0, "Ω"),
            ("xl", "Carga: X em f0 (R+L ou R+C)", -80.0, "Ω"),
            ("z0", "Z0 da linha / fonte", 50.0, "Ω"),
            ("f0", "Frequência de casamento f0", 2.0, "GHz"),
            ("ponta", "Stubs", PONTAS, ""),
            ("sol", "Solução (1 a 4)", 1, ""),
            ("q", "Q da rede Pi/T (> Q da rede L)", 5.0, ""),
            ("dd", "Stub duplo: espaçamento", 0.125, "λ"),
            ("d0", "Stub duplo: carga → 1º stub", 0.0, "λ"),
            ("n", "Nº de seções (transformador)", 3, ""),
            ("gm", "Γm (Chebyshev, Klopfenstein, banda)", 0.05, ""),
            ("lt", "Comprimento do taper", 1.0, "λg"),
            ("bw", "Banda p/ Bode-Fano (±)", 10.0, "%"),
            ("rvia", "Raio das vias (stubs em curto)", 0.2, "mm"),
            ("fmin", "f mín", 1.0, "GHz"),
            ("fmax", "f máx", 3.0, "GHz"),
        ], titulo="Carga e rede de casamento")
        self.form.pack(fill="x")
        self.corrigir = tk.BooleanVar(value=True)
        self.perdas = tk.BooleanVar(value=False)
        ttk.Checkbutton(self.painel, text="Corrigir ponta aberta, via e degraus (Garg-Bahl)",
                        variable=self.corrigir).pack(anchor="w", pady=(6, 0))
        ttk.Checkbutton(self.painel, text="Perdas na microstrip (condutor e dielétrico)",
                        variable=self.perdas).pack(anchor="w")
        self.botoes(("Calcular", self.calcular), ("Enviar p/ exportação", self.publicar),
                    ("Ex. 5.2: stub em curto", self.ex_stub), ("Ex. 5.4: stub duplo", self.ex_duplo),
                    ("Ex. 5.7: Chebyshev", self.ex_cheb), ("Ex. 5.8: Klopfenstein", self.ex_klop),
                    ("Steer Ex. 8.7: Pi, Q = 10", self.ex_pi),
                    ("Exportar S2P", self.exportar_s2p))
        self.saida = AreaTexto(self.painel, 22)
        self.saida.pack(fill="both", expand=True)
        self.saida.escrever(
            "Pozar, cap. 5. A carga é R em série com L (X > 0) ou C (X < 0),\n"
            "com ZL(f0) = R + jX, como nos exemplos do livro.\n\n"
            "Curvas: linhas TEM ideais (tracejado) e microstrip no substrato\n"
            "atual (contínuo), com ponta aberta real e, se marcado, perdas.\n"
            "Transformadores e tapers casam cargas reais (X é ignorado no\n"
            "projeto, mas entra na resposta).")
        self._ultimo = None
        self.observar_previa(self.form)
        estado.ouvir(lambda ev: ev == "substrato" and self.atualizar_previa())
        self.atualizar_previa()

    # ------------------------------------------------------------------ exemplos
    def ex_stub(self):
        self.form.set(dict(metodo=METODOS[0], rl=60, xl=-80, z0=50, f0=2, ponta=PONTAS[1], sol=1, fmin=1, fmax=3))
        self.calcular()

    def ex_duplo(self):
        self.form.set(dict(metodo=METODOS[2], rl=60, xl=-80, z0=50, f0=2, ponta=PONTAS[0], sol=2, dd=0.125,
                           d0=0, fmin=1, fmax=3))
        self.calcular()

    def ex_cheb(self):
        self.form.set(dict(metodo=METODOS[5], rl=100, xl=0, z0=50, f0=2, n=3, gm=0.05, fmin=0.2, fmax=3.8))
        self.calcular()

    def ex_pi(self):
        self.form.set(dict(metodo=METODOS[9], rl=500, xl=0, z0=50, f0=1, q=10, sol=1, fmin=0.5, fmax=1.5))
        self.calcular()

    def ex_klop(self):
        self.form.set(dict(metodo=METODOS[8], rl=50, xl=0, z0=100, f0=2, gm=0.02, lt=1.0, fmin=0.2, fmax=8))
        self.calcular()

    # ------------------------------------------------------------------ projeto
    def _projeto(self):
        v = self.form.get()
        if v["f0"] <= 0 or v["z0"] <= 0 or v["rl"] <= 0:
            raise ValueError("f0, Z0 e R devem ser positivos.")
        if not 0 < v["fmin"] < v["fmax"]:
            raise ValueError("Faixa de frequência inválida.")
        sub = self.estado.substrato
        f0, z0 = v["f0"] * 1e9, v["z0"]
        zl0 = complex(v["rl"], v["xl"])
        met = v["metodo"]
        curto = v["ponta"] == PONTAS[1]
        p = dict(v=v, sub=sub, f0=f0, z0=z0, zl0=zl0, met=met, curto=curto, texto=[], lay=None)
        w0 = float(sintese(z0, sub))
        _, ee0 = analise(w0, sub)
        lg = C0 / (f0 * np.sqrt(float(ee0))) * 1e3
        p.update(w0=w0, ee0=float(ee0), lg=lg)
        dl = extensao_ponta_aberta(w0, sub)
        if not self.corrigir.get():
            corr = 0.0
        elif curto:
            corr = lp.encurtamento_via(v["rvia"], lp.Geometria("Microstrip", w0), sub, f0) if v["rvia"] > 0 else 0.0
        else:
            corr = dl
        if met in METODOS[:3]:
            p["dl"] = dl
            if met == METODOS[2]:
                sols = cs.stub_duplo(zl0, z0, v["dd"], curto, v["d0"])
            else:
                sols = cs.stub_simples(zl0, z0, serie=met == METODOS[1], curto=curto)
            k = min(max(v["sol"], 1), len(sols)) - 1
            p["sols"], p["s"] = sols, sols[k]
            s = sols[k]
            if met == METODOS[2]:
                p["texto"] += [f"Stub duplo ({v['ponta'].lower()}s), espaçamento {v['dd']:g}λ, carga → 1º stub "
                               f"{v['d0']:g}λ (eq. 5.22):"]
                for i, so in enumerate(sols, 1):
                    p["texto"].append(f"  sol. {i}: b1 = {so.b1:+.3f}  b2 = {so.b2:+.3f}  ℓ1 = {so.l1:.3f}λ  "
                                      f"ℓ2 = {so.l2:.3f}λ")
                stubs = [dict(x=0.0, w=w0, comp=s.l2 * lg - corr, curto=curto),
                         dict(x=v["dd"] * lg, w=w0, comp=s.l1 * lg - corr, curto=curto)]
                d0 = v["d0"] * lg
                if d0 < w0:
                    p["texto"].append(f"  (layout: trecho de {w0 - d0:.3f} mm entre o 2º stub e o plano da carga "
                                      "para o stub caber na placa; incluído no modelo microstrip)")
                    d0 = w0
                segs = [(v["dd"] * lg, w0), (d0, w0)]
            else:
                tipo = "série" if met == METODOS[1] else "derivação"
                p["texto"] += [f"Stub simples em {tipo} ({v['ponta'].lower()}):"]
                for i, so in enumerate(sols, 1):
                    p["texto"].append(f"  sol. {i}: d = {so.d:.3f}λ  ℓ = {so.l:.3f}λ  "
                                      f"({'X' if tipo == 'série' else 'B'} do stub = {so.reativo:+.4g})")
                stubs = [dict(x=0.0, w=w0, comp=s.l * lg - corr, curto=curto)]
                d = s.d * lg
                if d < w0:
                    d += lg / 2
                    p["texto"].append("  (d < W: acrescentado λg/2 entre stub e carga, o que não muda o casamento "
                                      "em f0)")
                segs = [(d, w0)]
            p["texto"].append(f"Escolhida a sol. {k + 1}. λg em f0 = {lg:.3f} mm (W = {w0:.4f} mm).")
            for i, st in enumerate(stubs, 1):
                p["texto"].append(f"  stub {i}: {st['comp'] + corr:.3f} mm" +
                                  (f" − Δl {corr:.3f} = {st['comp']:.3f} mm" if corr else ""))
            if met != METODOS[1]:
                if min(st["comp"] for st in stubs) <= 0:
                    raise ValueError("Stub mais curto que a correção de ponta: escolha a outra solução.")
                p["stubs"], p["segs"] = stubs, segs
                info = {"tipo": f"casamento - {met}", "ZL_ohm": str(zl0), "f0_GHz": v["f0"], "W_mm": f"{w0:.4f}"}
                p["lay"] = circuitos.linha_com_stubs("casamento_stub", sub, segs, stubs, w0, w0, 5.0, 0.0,
                                                     raio_via=v["rvia"], info=info)
            else:
                p["texto"].append("Stub em série não é realizável em microstrip (Pozar o usa em coaxial/CPW):"
                                  " sem layout.")
        elif met == METODOS[3]:
            sols = cs.rede_l(zl0, z0, f0)
            k = min(max(v["sol"], 1), len(sols)) - 1
            p["sols"], p["s"] = sols, sols[k]
            p["texto"] += ["Rede L (eqs. 5.3 e 5.6):"]
            for i, so in enumerate(sols, 1):
                p["texto"].append(f"  sol. {i} ({so.topologia}): x = {so.x / z0:+.3f}, b = {so.b * z0:+.3f} → "
                                  f"série {cs.fmt_comp(so.serie)}, derivação {cs.fmt_comp(so.derivacao)}")
            p["texto"].append(f"Escolhida a sol. {k + 1} (componentes SMD: sem layout de linhas).")
        elif met in METODOS[9:]:
            pi = met == METODOS[9]
            sols = (cs.rede_pi if pi else cs.rede_t)(zl0, z0, f0, v["q"])
            k = min(max(v["sol"], 1), len(sols)) - 1
            p["sols"], p["s"] = sols, sols[k]
            rmaior, rmenor = max(zl0.real, z0), min(zl0.real, z0)
            p["texto"] += [f"Rede {'Pi' if pi else 'T'} = duas redes L através de RV = "
                           + (f"Rmaior/(Q² + 1) = {rmaior:g}/({v['q']:g}² + 1)" if pi else
                              f"Rmenor·(Q² + 1) = {rmenor:g}·({v['q']:g}² + 1)") + f" = {sols[0].rv:.4g} Ω",
                           f"  Q mínimo (rede L) = {cs.q_minimo(zl0, z0):.3f}; banda a 3 dB ≈ f0/Q = "
                           f"{v['f0'] / v['q'] * 1e3:.4g} MHz (o Q maior estreita a banda)"]
            nomes = ("B1", "X2", "B3") if pi else ("X1", "B2", "X3")
            for i, so in enumerate(sols, 1):
                val = "  ".join(f"{n} = {x:+.5g} {'S' if n[0] == 'B' else 'Ω'}" for n, x in zip(nomes, so.valores))
                p["texto"].append(f"  sol. {i}: {val}")
                p["texto"].append(f"          {so.descricao()}")
            p["texto"].append(f"Escolhida a sol. {k + 1} (fonte → carga; componentes SMD: sem layout de linhas). "
                              "Como no Steer, as quatro soluções diferem em DC, harmônicos e parasitas.")
        else:
            rl = zl0.real
            if v["xl"] != 0:
                p["texto"].append(f"⚠ X = {v['xl']:g} Ω ignorado no projeto (transformador casa carga real).")
            if met in METODOS[4:6]:
                if met == METODOS[4]:
                    zs = cs.binomial(z0, rl, v["n"])
                    bw = cs.banda_binomial(z0, rl, v["n"], v["gm"])
                    p["texto"].append(f"Binomial N = {v['n']} (eq. 5.53); banda |Γ| ≤ {v['gm']:g}: {100 * bw:.1f}%")
                else:
                    zs, tm = cs.chebyshev(z0, rl, v["n"], v["gm"])
                    bw = cs.banda_chebyshev(tm)
                    p["texto"].append(f"Chebyshev N = {v['n']}, Γm = {v['gm']:g} (eqs. 5.61-5.64): θm = "
                                      f"{np.degrees(tm):.1f}°, banda {100 * bw:.1f}%")
                trechos = []
                for i, z in enumerate(zs, 1):
                    w = float(sintese(z, sub))
                    _, ee = analise(w, sub)
                    comp = C0 / (4 * f0 * np.sqrt(float(ee))) * 1e3
                    trechos.append((comp, w))
                    p["texto"].append(f"  Z{i} = {z:.2f} Ω  W = {w:.4f} mm  λ/4 = {comp:.3f} mm")
                p["zs"], p["trechos"] = zs, trechos
                w_l = float(sintese(rl, sub))
                if self.corrigir.get():
                    trechos, dls = fc.correcao_degraus(trechos, sub, w0, w_l)
                    p["texto"].append("  degraus (Garg-Bahl 3.4.3): ℓ corrigidos = "
                                      + ", ".join(f"{c:.3f}" for c, _ in trechos) + " mm")
                p["trechos_fis"] = trechos
                x, ws = [], []
                pos = 0.0
                for comp, w in trechos:
                    x += [pos, pos + comp]
                    ws += [w, w]
                    pos += comp
            else:
                tipo = met.split()[-1].capitalize()
                n_pts = 201
                s_rel, zz = cs.perfil_taper(tipo, z0, rl, n_pts, v["gm"])
                ws = np.array([float(sintese(z, sub)) for z in zz])
                ee = np.array([float(analise(w, sub)[1]) for w in ws])
                l_mm = v["lt"] * C0 / (f0 * np.mean(np.sqrt(ee))) * 1e3
                x = s_rel * l_mm
                p["zs_taper"], p["ee_taper"] = zz, ee
                p["texto"].append(f"Taper {tipo.lower()} de {v['lt']:g} λg ({l_mm:.3f} mm), Z de {zz[0]:.2f} a "
                                  f"{zz[-1]:.2f} Ω, W de {ws[0]:.3f} a {ws[-1]:.3f} mm")
                if tipo == "Klopfenstein":
                    a = cs.parametro_a(z0, rl, v["gm"])
                    p["texto"].append(f"  A = {a:.3f}: banda passante para βL ≥ A, isto é, f ≥ "
                                      f"{a / (2 * np.pi * v['lt']) * v['f0']:.3f} GHz (eq. 5.78); degraus de "
                                      f"Γm nas pontas")
            w_out = float(sintese(rl, sub))
            info = {"tipo": f"casamento - {met}", "Z0_ohm": z0, "RL_ohm": rl, "f0_GHz": v["f0"]}
            p["x_perfil"], p["w_perfil"] = np.asarray(x), np.asarray(ws)
            p["lay"] = circuitos.linha_perfil("casamento_" + met.split()[-1].lower(), sub, x, ws, w0, w_out, 5.0,
                                              info=info)
        return p

    def _ideal(self, p, f):
        v, z0, f0 = p["v"], p["z0"], p["f0"]
        met = p["met"]
        if met in METODOS[:2]:
            return cs.abcd_stub_simples(f, f0, z0, p["s"], p["curto"], met == METODOS[1])
        if met == METODOS[2]:
            return cs.abcd_stub_duplo(f, f0, z0, p["s"], v["dd"], p["curto"], v["d0"])
        if met == METODOS[3] or met in METODOS[9:]:
            return p["s"].abcd(f)
        if met in METODOS[4:6]:
            return cs.abcd_multissecao(f, f0, p["zs"])
        zz = p["zs_taper"]
        th = np.pi * 2 * v["lt"] / (len(zz) - 1)                       # trechos de βL/(n-1) em f0
        zm = 0.5 * (zz[1:] + zz[:-1])
        return abcd.cascata(*[abcd.linha(z, th * np.asarray(f, dtype=float) / f0) for z in zm])

    def _microstrip(self, p, f):
        """Mesma rede com microstrip real (ponta aberta = comprimento físico + Δl)."""
        met, sub = p["met"], p["sub"]
        perdas = self.perdas.get()
        if met in (METODOS[1], METODOS[3]) or met in METODOS[9:]:
            return None
        if met in METODOS[:3]:
            geo = lp.Geometria("Microstrip", p["w0"])
            mats = []
            stubs = p["stubs"]
            for i, st in enumerate(stubs):
                comp = st["comp"] + (0.0 if st["curto"] else p["dl"])
                mats.append(lp.abcd_stub(f, geo, sub, comp, curto=st["curto"], perdas=perdas,
                                         raio_via=p["v"]["rvia"] if st["curto"] and p["v"]["rvia"] > 0 else None))
                if i < len(p["segs"]):
                    mats.append(lp.abcd_trecho(f, geo, sub, p["segs"][i][0], perdas))
            return abcd.cascata(*mats)
        x, w = p["x_perfil"], p["w_perfil"]
        if met in METODOS[4:6]:
            return fc.abcd_degraus_microstrip(f, p["trechos_fis"], sub, perdas, True, p["w0"],
                                              float(sintese(p["zl0"].real, sub)))
        else:
            comps = list(np.diff(x) * 1e-3)
            ws = list(0.5 * (w[1:] + w[:-1]))
        m = None
        for dl, wi in zip(comps, ws):
            z, ee, a = (q[:, 0] for q in parametros_linha(wi, sub, f, perdas, False))
            th = 2 * np.pi * np.asarray(f) * np.sqrt(ee) / C0 * dl - 1j * a * dl
            mi = abcd.linha(z, th)
            m = mi if m is None else m @ mi
        return m

    def atualizar_previa(self):
        try:
            p = self._projeto()
        except ValueError as e:
            self.previa.mostrar_erro(str(e))
            return
        if p["lay"] is None:
            self.previa.mostrar_erro("Sem layout para este método (componentes concentrados ou stub série).")
        else:
            self.previa.mostrar(p["lay"])

    @protegido
    def calcular(self):
        p = self._projeto()
        v, z0, f0 = p["v"], p["z0"], p["f0"]
        if p["lay"] is not None:
            self.previa.mostrar(p["lay"])
        f = np.linspace(v["fmin"] * 1e9, v["fmax"] * 1e9, 1201)
        zl = cs.carga_serie(p["zl0"], f0)(f)
        g_id = cs.gama(cs.zin_abcd(self._ideal(p, f), zl), z0)
        mm = self._microstrip(p, f)
        g_ms = None if mm is None else cs.gama(cs.zin_abcd(mm, zl), z0)
        gl = cs.gama(zl, z0)
        fg = f / 1e9
        ax = self.grafico.eixos(1, 2)
        ax[0].plot(fg, np.abs(gl), color="#999", lw=0.8, label="sem casamento")
        ax[0].plot(fg, np.abs(g_id), ls="--", label="TEM ideal")
        if g_ms is not None:
            ax[0].plot(fg, np.abs(g_ms), label="microstrip")
        ax[0].axhline(v["gm"], color="gray", ls=":", lw=0.8)
        ax[0].axvline(v["f0"], color="gray", ls=":", lw=0.8)
        ax[0].set_ylim(0, 1.02)
        ax[0].set_xlabel("f [GHz]")
        ax[0].set_ylabel("|Γ| na entrada")
        ax[0].set_title(p["met"], fontsize=9)
        ax[0].grid(True, alpha=0.3)
        ax[0].legend(fontsize=7)
        desenhar_smith(ax[1])
        ax[1].plot(gl.real, gl.imag, color="#999", lw=0.8)
        ax[1].plot(g_id.real, g_id.imag, ls="--", lw=1)
        if g_ms is not None:
            ax[1].plot(g_ms.real, g_ms.imag, lw=1.2)
        i0 = int(np.argmin(np.abs(f - f0)))
        ax[1].plot(gl[i0].real, gl[i0].imag, "ks", ms=5, label="ZL(f0)")
        ax[1].plot(g_id[i0].real, g_id[i0].imag, "o", color="C0", ms=5, label="Zin(f0)")
        ax[1].set_title("Carta de Smith (Z0 normalizada)", fontsize=9)
        ax[1].legend(fontsize=7, loc="lower left")
        self.grafico.desenhar()
        mrede = mm if mm is not None else self._ideal(p, f)
        self._ultimo = (f, abcd.abcd_para_s(mrede, z0), z0, p["met"])
        linhas = [f"Carga: {cs.descrever_carga(p['zl0'], f0)} (ZL(f0) = {p['zl0'].real:g} "
                  f"{p['zl0'].imag:+g}j Ω), Z0 = {z0:g} Ω", ""] + p["texto"] + [""]
        for nome, g in (("ideal", g_id), ("microstrip", g_ms)):
            if g is None:
                continue
            ok = np.abs(g) <= v["gm"]
            fx = _faixa(f, ok, f0)
            linhas.append(f"|Γ| ≤ {v['gm']:g} ({nome}): " + (f"{fx[0] / 1e9:.3f}-{fx[1] / 1e9:.3f} GHz "
                                                             f"({100 * (fx[1] - fx[0]) / f0:.1f}%)" if fx else
                                                             "não atende em f0"))
            linhas.append(f"  |Γ(f0)| = {abs(g[i0]):.4f}")
        try:
            bf = cs.bode_fano(p["zl0"], f0, f0 * (1 - v["bw"] / 100), f0 * (1 + v["bw"] / 100))
            linhas += [""] + bf.texto
        except ValueError:
            pass
        linhas += ["", "Limites do substrato:"] + texto_limites(p["sub"], p["w0"], v["fmax"] * 1e9)
        linhas += ["", "Junções T dos stubs e degraus não compensados: ajustar em EM."]
        self.saida.escrever("\n".join(linhas))

    @protegido
    def publicar(self):
        p = self._projeto()
        if p["lay"] is None:
            raise ValueError("Este método não gera layout de linhas.")
        self.previa.mostrar(p["lay"])
        self.estado.publicar_layout(p["lay"])
        self.estado.msg(f"Layout '{p['lay'].nome}' enviado para Exportar (Ctrl+E).")

    def exportar_s2p(self):
        if not self._ultimo:
            self.calcular()
        if self._ultimo:
            f, s, z0, met = self._ultimo
            self.salvar_s2p("casamento", f, s, z0, f"rede de casamento: {met} (porta 2 = plano da carga)")


def _faixa(f, ok, f0):
    i = int(np.argmin(np.abs(f - f0)))
    if not ok[i]:
        return None
    a = b = i
    while a > 0 and ok[a - 1]:
        a -= 1
    while b < len(f) - 1 and ok[b + 1]:
        b += 1
    return f[a], f[b]
