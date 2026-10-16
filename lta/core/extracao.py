"""Extração de parâmetros de linhas metamateriais ressonantes a partir de parâmetros S.

Referência: F. Martín (2015), Apêndice G (eqs. G.1-G.19) e seção 3.5.2.
Os dados (simulação EM sem perdas ou medição) são de UMA célula simétrica.
As condições usadas são as do livro:
  - interseção de S11 com o círculo de resistência unitária (r = 1): o ramo em
    derivação abre (T) ou o ramo em derivação do π abre;
  - interseção com o círculo de condutância unitária (g = 1): o ramo série se anula;
  - zero de transmissão (mín. |S21|) e zero de reflexão (mín. |S11|);
  - βl = ±π/2, obtido de cos βl = (A+D)/2 da matriz ABCD medida.
Todas as frequências podem ser fixadas manualmente (dicionário `manual`).
"""
from dataclasses import dataclass, field

import numpy as np

from . import abcd

# ---------------------------------------------------------------------------
# Modelos de circuito (unidades SI)
# ---------------------------------------------------------------------------
def _t(zs_braco, yp):
    return abcd.cascata(abcd.series(zs_braco), abcd.shunt(yp), abcd.series(zs_braco))


def _pi(zs, y_lado):
    return abcd.cascata(abcd.shunt(y_lado), abcd.series(zs), abcd.shunt(y_lado))


def _w(f):
    return 2 * np.pi * np.asarray(f, dtype=float)


def _tanque(w, l, c):
    """Impedância do tanque LC paralelo."""
    return 1j * w * l / (1 - w**2 * l * c)


def abcd_csrr(f, p):
    """Fig. 3.38a/b: braços L/2 (+ 2Cg); derivação C em série com (Lc ‖ Cc)."""
    w = _w(f)
    braco = 1j * w * p["L"] / 2 + (1 / (2j * w * p["Cg"]) if p.get("Cg") else 0)
    zp = 1 / (1j * w * p["C"]) + _tanque(w, p["Lc"], p["Cc"])
    return _t(braco * np.ones_like(w), 1 / zp)


def abcd_srr(f, p):
    """Fig. 3.35b: π com C/2 (‖ Lp') nos lados; série L' + (Ls' ‖ Cs')."""
    w = _w(f)
    y = 1j * w * p["C"] / 2 + (1 / (1j * w * p["Lp"]) if p.get("Lp") else 0)
    zs = 1j * w * p["L"] + _tanque(w, p["Ls"], p["Cs"])
    return _pi(zs, y * np.ones_like(w))


def abcd_osrr(f, p):
    """Fig. 3.53e: π com C nos lados; série Ls + Cs (OSRR)."""
    w = _w(f)
    zs = 1j * w * p["Ls"] + 1 / (1j * w * p["Cs"])
    return _pi(zs, 1j * w * p["C"])


def abcd_ocsrr(f, p):
    """Fig. 3.54e: T com L nos braços; derivação Lp ‖ Cp (OCSRR)."""
    w = _w(f)
    yp = 1j * w * p["Cp"] + 1 / (1j * w * p["Lp"])
    return _t(1j * w * p["L"], yp)


def abcd_crlh(f, p):
    """Fig. 3.20 (CL-loaded, Fig. 3.27): braços LR/2 + 2CL; derivação LL ‖ CR."""
    w = _w(f)
    braco = 1j * w * p["LR"] / 2 + 1 / (2j * w * p["CL"])
    yp = 1j * w * p["CR"] + 1 / (1j * w * p["LL"])
    return _t(braco, yp)


@dataclass
class Modelo:
    nome: str
    figura: str
    parametros: tuple        # chaves (L em H, C em F)
    abcd: callable
    descricao: str


MODELOS = {
    "CSRR (microstrip)": Modelo("CSRR (microstrip)", "Fig. 3.38a", ("L", "C", "Lc", "Cc"), abcd_csrr,
                                "Microstrip com CSRR no terra: braços L/2, derivação C + (Lc ‖ Cc)."),
    "CSRR + gap (microstrip)": Modelo("CSRR + gap (microstrip)", "Fig. 3.38b", ("L", "Cg", "C", "Lc", "Cc"),
                                      abcd_csrr, "CSRR e gap série: braços L/2 + 2Cg (CRLH)."),
    "SRR (CPW)": Modelo("SRR (CPW)", "Fig. 3.35b sem Lp", ("C", "L", "Ls", "Cs"), abcd_srr,
                        "CPW com SRRs: π com C/2; série L' + (Ls' ‖ Cs')."),
    "SRR + fita (CPW)": Modelo("SRR + fita (CPW)", "Fig. 3.35b", ("C", "Lp", "L", "Ls", "Cs"), abcd_srr,
                               "CPW com SRRs e fitas em derivação: π com C/2 ‖ Lp' (CRLH)."),
    "OSRR (CPW)": Modelo("OSRR (CPW)", "Fig. 3.53e", ("C", "Ls", "Cs"), abcd_osrr,
                         "CPW com OSRR em série: π com C; série Ls + Cs."),
    "OCSRR (CPW)": Modelo("OCSRR (CPW)", "Fig. 3.54e", ("L", "Lp", "Cp"), abcd_ocsrr,
                          "CPW com OCSRR em derivação: braços L; derivação Lp ‖ Cp."),
    "CRLH (CL-loaded)": Modelo("CRLH (CL-loaded)", "Fig. 3.20", ("LR", "CR", "LL", "CL"), abcd_crlh,
                               "Interdigital + stub com via (Fig. 3.27): braços LR/2 + 2CL, derivação LL ‖ CR."),
}


# ---------------------------------------------------------------------------
# Grandezas medidas
# ---------------------------------------------------------------------------
def _cruzamentos(f, y):
    """Frequências (interpoladas) em que y troca de sinal."""
    y = np.asarray(y, dtype=float)
    i = np.flatnonzero(np.sign(y[:-1]) * np.sign(y[1:]) < 0)
    return f[i] - y[i] * (f[i + 1] - f[i]) / (y[i + 1] - y[i])


def _zeros(f, y):
    """Zeros de y (troca de sinal com |y| decrescendo até o cruzamento), excluindo polos."""
    y = np.asarray(y, dtype=float)
    a = np.abs(y)
    out = []
    for i in np.flatnonzero(np.sign(y[:-1]) * np.sign(y[1:]) < 0):
        viz = a[max(i - 1, 0)], a[min(i + 2, len(y) - 1)]
        if max(a[i], a[i + 1]) <= min(viz) * 1.000001 or (i == 0 or i + 2 >= len(y)):
            out.append(f[i] - y[i] * (f[i + 1] - f[i]) / (y[i + 1] - y[i]))
    return np.array(out)


def _unir(*listas, tol=1e-4):
    """Une listas de frequências, descartando quase-repetidas (tolerância relativa)."""
    out = []
    for x in sorted(np.concatenate([np.atleast_1d(l) for l in listas])):
        if not out or abs(x / out[-1] - 1) > tol:
            out.append(x)
    return np.array(out)


def _minimo(f, y):
    """Mínimo de y com refinamento parabólico."""
    i = int(np.argmin(y))
    if 0 < i < len(f) - 1:
        a, b, c = y[i - 1:i + 2]
        den = a - 2 * b + c
        if den > 0:
            return f[i] + 0.5 * (a - c) / den * (f[1] - f[0])
    return f[i]


@dataclass
class Medida:
    f: np.ndarray
    s: np.ndarray
    z0: float = 50.0
    zin: np.ndarray = field(init=False)
    cos_bl: np.ndarray = field(init=False)

    def __post_init__(self):
        ok = np.all(np.isfinite(self.s.reshape(len(self.f), -1)), axis=1)
        self.f, self.s = self.f[ok], self.s[ok]          # descarta pontos singulares
        s11 = self.s[:, 0, 0]
        self.zin = self.z0 * (1 + s11) / (1 - s11)
        m = abcd.s_para_abcd(self.s, self.z0)
        a, b, c = 0.5 * (m[:, 0, 0] + m[:, 1, 1]), m[:, 0, 1], m[:, 1, 0]
        self.cos_bl = a
        with np.errstate(all="ignore"):
            self.zb = np.sqrt(b / c + 0j)                  # ZB da célula simétrica
            # equivalentes da célula simétrica: T (braço = (A-1)/C, Yp = C) e π (Zs = B, Y = (A-1)/B)
            self.t_braco, self.t_yp = (a - 1) / c, c
            self.pi_zs, self.pi_y = b, (a - 1) / b

    def em(self, f0, y):
        """Interpola y (complexo) em f0."""
        return np.interp(f0, self.f, np.real(y)) + 1j * np.interp(f0, self.f, np.imag(y))

    def derivacao_abre(self, topo):
        """Ramo em derivação abre: S11 sobre o círculo r = 1.

        Candidatos: cruzamentos de S11 com o círculo (critério do livro, adequado a
        dados distribuídos) e zeros do equivalente da célula (T: Yp = 0; π: Y = 0),
        que pegam os pontos em que S11 só tangencia o círculo.
        """
        return _unir(_cruzamentos(self.f, self.zin.real / self.z0 - 1),
                     _zeros(self.f, (self.t_yp if topo == "T" else self.pi_y).imag))

    def serie_anula(self, topo):
        """Ramo série se anula: S11 sobre o círculo g = 1 (cruzamentos e zeros do equivalente)."""
        return _unir(_cruzamentos(self.f, (self.z0 / self.zin).real - 1),
                     _zeros(self.f, (self.t_braco if topo == "T" else self.pi_zs).imag))

    def zeros_transmissao(self, max_n=3):
        """Mínimos locais internos de |S21| (notches), do mais profundo ao mais raso.

        Os extremos da faixa ficam de fora: com capacitor série, |S21| → 0 em DC sem
        que haja zero de transmissão.
        """
        y = 20 * np.log10(np.abs(self.s[:, 1, 0]) + 1e-15)
        idx = [i for i in range(1, len(y) - 1) if y[i] < y[i - 1] and y[i] <= y[i + 1]]
        idx = sorted(idx, key=lambda i: y[i])[:max_n]
        return [_minimo(self.f[i - 1:i + 2], y[i - 1:i + 2]) for i in idx]

    def bl_pi2(self):
        """Frequências com βl = ±π/2 (cos βl = 0) na banda passante."""
        return _cruzamentos(self.f, self.cos_bl.real)


# ---------------------------------------------------------------------------
# Extrações (Apêndice G). Cada função recebe as frequências já escolhidas (Hz).
# ---------------------------------------------------------------------------
def _csrr(m: Medida, fq, gap):
    fz, fo, f1 = fq["fz"], fq["fo"], fq["f1"]
    wo, wz, w1 = 2 * np.pi * fo, 2 * np.pi * fz, 2 * np.pi * f1
    chi = m.em(fo, m.zin).imag
    cond = {"fz [GHz] (mín. |S21|)": fz / 1e9, "fo [GHz] (r = 1)": fo / 1e9, "χ(fo) [Ω]": chi}
    if gap:
        ws = 2 * np.pi * fq["fs"]
        if abs(ws / wo - 1) > 0.05:
            L = chi / (wo - ws**2 / wo)                   # χ = ωo·L - 1/(ωo·Cg), ωs² = 1/(L·Cg)
        else:
            # quase balanceada: χ(fo) ≈ 0 não determina L. Em βl = ±π/2, A = 0 => Zp = -braço
            # e ZB² = braço·Zp = X², logo ZB(ω1) = |X do braço| = L·|ω1 - ωs²/ω1|/2.
            zb1 = abs(m.em(f1, m.zb).real)
            L = 2 * zb1 / abs(w1 - ws**2 / w1)
            cond["ZB(f1) [Ω] (balanceada)"] = zb1
        cg = 1 / (ws**2 * L)
        cond["fs [GHz] (g = 1)"] = fq["fs"] / 1e9
    else:
        L = chi / wo                                      # eq. G.2
        cg = None
    cond["f(βl = ±π/2) [GHz]"] = f1 / 1e9
    braco = w1 * L / 2 - (1 / (2 * w1 * cg) if cg else 0)     # Im do braço série em ω1
    g = w1 / (1 - w1**2 / wo**2) - 1 / (w1 * (1 / wz**2 - 1 / wo**2))
    lc = -braco / g                                           # braço = -Zp (eqs. G.5 / G.6)
    p = {"L": L, "C": (1 / wz**2 - 1 / wo**2) / lc, "Lc": lc, "Cc": 1 / (wo**2 * lc)}
    if gap:
        p["Cg"] = cg
    return p, cond


def _srr(m: Medida, fq, fita):
    fz, fs, f1 = fq["fz"], fq["fs"], fq["f1"]
    ws, wz, w1 = 2 * np.pi * fs, 2 * np.pi * fz, 2 * np.pi * f1
    b = (1 / m.em(fs, m.zin)).imag
    cond = {"fz [GHz] (mín. |S21|)": fz / 1e9, "fs [GHz] (g = 1)": fs / 1e9, "B(fs) [mS]": b * 1e3}
    if fita:
        wp = 2 * np.pi * fq["fp"]
        C = b / (ws - wp**2 / ws)                             # B = 2(ωC/2 - 1/(ωLp'))
        lp = 2 / (wp**2 * C)                                  # eq. G.11
        cond["fp [GHz] (r = 1)"] = fq["fp"] / 1e9
    else:
        C = b / ws                                            # eq. G.9
        lp = None
    cond["f(βl = ±π/2) [GHz]"] = f1 / 1e9
    b1 = w1 * C / 2 - (1 / (w1 * lp) if lp else 0)            # Im da admitância de cada lado em ω1
    # Zs(ωs) = 0 e Zs(ω1)·Y(ω1) = -1, com Zs = jω[L' + Ls'/(1 - ω²/ωz²)]
    ls = 1 / (w1 * b1) / (1 / (1 - w1**2 / wz**2) - 1 / (1 - ws**2 / wz**2))
    p = {"C": C, "L": -ls / (1 - ws**2 / wz**2), "Ls": ls, "Cs": 1 / (wz**2 * ls)}
    if fita:
        p["Lp"] = lp
    return p, cond


def _osrr(m: Medida, fq):
    fs, fr = fq["fs"], fq["frefl"]
    ws, w = 2 * np.pi * fs, 2 * np.pi * fr
    b = (1 / m.em(fs, m.zin)).imag
    C = b / (2 * ws)                                          # eq. G.12
    # S11(fr) = 0  <=>  ZB(fr) = Z0 no π: x(1 + Z0²ω²C²) = 2Z0²ωC, x = reatância série
    x = 2 * m.z0**2 * w * C / (1 + (m.z0 * w * C) ** 2)
    ls = x / (w - ws**2 / w)
    cond = {"fs [GHz] (g = 1)": fs / 1e9, "B(fs) [mS]": b * 1e3, "f zero de reflexão [GHz]": fr / 1e9}
    return {"C": C, "Ls": ls, "Cs": 1 / (ws**2 * ls)}, cond


def _ocsrr(m: Medida, fq):
    fp, fr = fq["fp"], fq["frefl"]
    wp, w = 2 * np.pi * fp, 2 * np.pi * fr
    chi = m.em(fp, m.zin).imag
    L = chi / (2 * wp)                                        # eq. G.16
    # S11(fr) = 0  <=>  ZB(fr) = Z0 no T: b(Z0² + ω²L²) = 2ωL, b = susceptância em derivação
    b = 2 * w * L / (m.z0**2 + (w * L) ** 2)
    cp = b / (w - wp**2 / w)
    cond = {"fp [GHz] (r = 1)": fp / 1e9, "χ(fp) [Ω]": chi, "f zero de reflexão [GHz]": fr / 1e9}
    return {"L": L, "Lp": 1 / (wp**2 * cp), "Cp": cp}, cond


def _crlh(m: Medida, fq):
    """fp (derivação abre) e χ(fp) -> LR, CL; fs (série nula) e B(fs) -> CR, LL.

    Quase balanceada (fs ≈ fp): χ e B se anulam; usa-se βl = ±π/2, onde A = 0:
    ZB(ω1) = |X do braço| e a susceptância em derivação vale 1/X.
    """
    fs, fp, f1 = fq["fs"], fq["fp"], fq["f1"]
    ws, wp, w1 = 2 * np.pi * fs, 2 * np.pi * fp, 2 * np.pi * f1
    cond = {"fs [GHz] (g = 1)": fs / 1e9, "fp [GHz] (r = 1)": fp / 1e9}
    if abs(ws / wp - 1) > 0.05:
        chi = m.em(fp, m.zin).imag
        b = (1 / m.em(fs, m.zin)).imag
        lr = chi / (wp - ws**2 / wp)                    # χ = 2·X_braço = LR(ω - ωs²/ω)
        cr = b / (ws - wp**2 / ws)                      # B = CR(ω - ωp²/ω)
        cond.update({"χ(fp) [Ω]": chi, "B(fs) [mS]": b * 1e3})
    else:
        zb1 = abs(m.em(f1, m.zb).real)
        x = zb1 * np.sign(w1 - ws**2 / w1)              # X do braço em ω1
        lr = 2 * x / (w1 - ws**2 / w1)
        cr = (1 / x) / (w1 - wp**2 / w1)
        cond.update({"f(βl = ±π/2) [GHz]": f1 / 1e9, "ZB(f1) [Ω] (balanceada)": zb1})
    return {"LR": lr, "CL": 1 / (ws**2 * lr), "CR": cr, "LL": 1 / (wp**2 * cr)}, cond


# modelo -> (função, frequências necessárias, topologia da célula)
_EXTRATORES = {
    "CSRR (microstrip)": (lambda m, fq: _csrr(m, fq, False), ("fz", "fo", "f1"), "T"),
    "CSRR + gap (microstrip)": (lambda m, fq: _csrr(m, fq, True), ("fz", "fo", "fs", "f1"), "T"),
    "SRR (CPW)": (lambda m, fq: _srr(m, fq, False), ("fz", "fs", "f1"), "pi"),
    "SRR + fita (CPW)": (lambda m, fq: _srr(m, fq, True), ("fz", "fs", "fp", "f1"), "pi"),
    "OSRR (CPW)": (_osrr, ("fs", "frefl"), "pi"),
    "OCSRR (CPW)": (_ocsrr, ("fp", "frefl"), "T"),
    "CRLH (CL-loaded)": (_crlh, ("fs", "fp", "f1"), "T"),
}
FREQUENCIAS = {"fz": "zero de transmissão (mín. |S21|)", "fo": "ramo em derivação abre (r = 1)",
               "fp": "ramo em derivação abre (r = 1)", "fs": "ramo série se anula (g = 1)",
               "f1": "βl = ±π/2", "frefl": "zero de reflexão (mín. |S11|)"}


def _candidatas(m: Medida, chave, topo, max_n=10):
    if chave == "fz":
        c = m.zeros_transmissao()
    elif chave in ("fo", "fp"):
        c = m.derivacao_abre(topo)
    elif chave == "fs":
        c = m.serie_anula(topo)
    elif chave == "f1":
        c = m.bl_pi2()
    else:
        y = 20 * np.log10(np.abs(m.s[:, 0, 0]) + 1e-15)
        idx = [i for i in range(1, len(y) - 1) if y[i] <= y[i - 1] and y[i] <= y[i + 1] and y[i] < -10]
        c = [_minimo(m.f[i - 1:i + 2], y[i - 1:i + 2]) for i in idx]
    return list(np.asarray(c, dtype=float)[:max_n])


def extrair(nome_modelo, f, s, z0=50.0, manual=None):
    """Extrai os parâmetros do modelo. Retorna (parâmetros, condições, erro rms de S).

    As frequências características têm em geral vários candidatos (a curva de S11
    cruza os círculos r = 1 e g = 1 mais de uma vez). Testam-se as combinações e
    fica a que dá todos os elementos positivos e o menor erro rms contra os dados.
    manual: frequências [Hz] fixadas pelo usuário (chaves de FREQUENCIAS).
    """
    if nome_modelo not in _EXTRATORES:
        raise ValueError(f"Modelo desconhecido: {nome_modelo}")
    func, chaves, topo = _EXTRATORES[nome_modelo]
    man = {k: v for k, v in (manual or {}).items() if v}
    m = Medida(np.asarray(f, dtype=float), np.asarray(s), z0)
    listas = []
    for k in chaves:
        c = [man[k]] if k in man else _candidatas(m, k, topo)
        if not c:
            raise ValueError(f"Não encontrei {FREQUENCIAS[k]} na faixa dos dados. "
                             "Amplie a faixa ou informe a frequência manualmente.")
        listas.append(c)
    melhor = None
    for comb in np.array(np.meshgrid(*listas, indexing="ij")).reshape(len(chaves), -1).T:
        fq = dict(zip(chaves, comb))
        try:
            with np.errstate(all="ignore"):
                p, cond = func(m, fq)
        except (ValueError, ZeroDivisionError):
            continue
        p = {k: float(np.real(v)) for k, v in p.items()}
        if not all(np.isfinite(v) and v > 0 for v in p.values()):
            continue
        err = erro_rms(nome_modelo, p, m.f, m.s, z0)
        if not np.isfinite(err):
            continue
        if melhor is None or err < melhor[2]:
            melhor = (p, cond, err)
    if melhor is None:
        raise ValueError("Nenhuma combinação de frequências deu elementos positivos: confira o modelo "
                         "escolhido, a faixa dos dados, ou informe as frequências manualmente.")
    return melhor


def s_modelo(nome_modelo, p, f, z0=50.0):
    return abcd.abcd_para_s(MODELOS[nome_modelo].abcd(f, p), z0)


def erro_rms(nome_modelo, p, f, s, z0=50.0):
    """Erro rms |S_modelo - S_dados| (S11 e S21), adimensional."""
    with np.errstate(all="ignore"):
        sm = s_modelo(nome_modelo, p, f, z0)
    d = np.concatenate([sm[:, 0, 0] - s[:, 0, 0], sm[:, 1, 0] - s[:, 1, 0]])
    d = d[np.isfinite(d)]
    return float(np.sqrt(np.mean(np.abs(d) ** 2))) if d.size else np.inf


# ---------------------------------------------------------------------------
# Transformação do modelo físico da CPW com SRR e fitas (eqs. 3.71-3.74)
# ---------------------------------------------------------------------------
def srr_cpw_transformado(L, M, Ls, Cs, Lp=np.inf):
    """Elementos do π equivalente (Fig. 3.35b) a partir do modelo físico (Fig. 3.35a).

    Retorna (L', Lp', Ls', Cs'). Validado contra as Tabelas 3.1 e 3.2 do livro.
    Na eq. 3.72, o fator (1 + M²/(2LpLs)) aparece ao quadrado: é o que reproduz a
    Tabela 3.1 e dá ωz = ωo/√(1 + M²/(2LpLs)) ≤ ωo, como afirma o texto.
    """
    wo2 = 1 / (Ls * Cs)
    a = 1 + L / (4 * Lp)
    b = 1 + M**2 / (2 * Lp * Ls)
    ls_ = 2 * M**2 * Cs * wo2 * a**2 / b
    cs_ = Ls / (2 * M**2 * wo2) * b**2 / a**2
    l_ = (2 + L / (2 * Lp)) * L / 2 - ls_
    lp_ = 2 * Lp + L / 2
    return l_, lp_, ls_, cs_
