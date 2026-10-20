"""Projeto robusto: tolerâncias, potência máxima e orçamento de erros da microstrip.

Referência: R. Garg, I. Bahl e M. Bianchi, Microstrip Lines and Slotlines, 4ª ed., 2024:
  - sensibilidades e pior caso (eqs. 2.97-2.105) e precisão de fabricação exigida (2.100-2.101);
    em linhas acopladas, eqs. 4.120-4.127; aqui as sensibilidades S = (B/A)·∂A/∂B são numéricas
    sobre os mesmos modelos usados no resto do programa;
  - potência média por aquecimento (eqs. 2.84-2.95) e de pico por ruptura (eq. 2.96);
  - orçamento de erros (Tabela 2.5): espessura, dispersão, tolerâncias, degrau e conectores.
"""
from dataclasses import dataclass, field, replace

import numpy as np

from .microstrip import ETA0, MU0, NP_DB, Substrato, alfa_condutor, alfa_dieletrico, analise, dispersao

# Tabela 2.2 (K em W/cm/°C convertido para W/m/°C; rigidez em kV/cm). Laminados de PCB: datasheets
# (conferir com o fornecedor).
MATERIAIS = {
    "Alumina (Tabela 2.2)": (30.0, 4000.0),
    "Safira (Tabela 2.2)": (40.0, 4000.0),
    "Quartzo fundido (Tabela 2.2)": (1.0, 10000.0),
    "GaAs (Tabela 2.2)": (30.0, 350.0),
    "Si 10³ Ω·cm (Tabela 2.2)": (90.0, 300.0),
    "BeO (Tabela 2.2)": (250.0, None),
    "SiC (Tabela 2.2)": (350.0, 2000.0),
    "FR-4 (Tabela 2.2)": (0.27, 200.0),
    "Poliestireno (Tabela 2.2)": (0.15, 280.0),
    "Rogers RO4003C (datasheet)": (0.71, 312.0),
    "Rogers RO3010 (datasheet)": (0.95, None),
    "Rogers RT5880 (datasheet)": (0.20, None),
}


# ---------------------------------------------------------------------------
# Sensibilidades e pior caso
# ---------------------------------------------------------------------------
def _ms(w, sub):
    z, e = analise(w, sub)
    return float(z), float(e)


def sensibilidades(w_mm, sub: Substrato, rel=1e-4):
    """{'W'|'h'|'εr': (S de Z0, S de εeff)} com S_B^A = (B/A)·∂A/∂B (eq. 2.99)."""
    z0, e0 = _ms(w_mm, sub)
    out = {}
    for nome, var in (("W", "w"), ("h", "h_mm"), ("εr", "er")):
        if var == "w":
            z1, e1 = _ms(w_mm * (1 + rel), sub)
        else:
            z1, e1 = _ms(w_mm, replace(sub, **{var: getattr(sub, var) * (1 + rel)}))
        out[nome] = ((z1 - z0) / z0 / rel, (e1 - e0) / e0 / rel)
    return out


@dataclass
class PiorCaso:
    dz_rel: float          # |ΔZ0/Z0| máximo
    de_rel: float          # |Δεeff/εeff| máximo
    vswr: float            # VSWR por ΔZ0 (linha de Z0 + ΔZ0 num sistema de Z0)
    sens: dict = field(default_factory=dict)


def pior_caso(w_mm, sub: Substrato, dw, dh, der):
    """Pior caso das eqs. 2.97-2.98 para tolerâncias ±ΔW [mm], ±Δh [mm], ±Δεr."""
    s = sensibilidades(w_mm, sub)
    rel = {"W": dw / w_mm, "h": dh / sub.h_mm, "εr": der / sub.er}
    dz = sum(abs(rel[k] * s[k][0]) for k in s)
    de = sum(abs(rel[k] * s[k][1]) for k in s)
    return PiorCaso(dz, de, 1 + dz, s)


def precisao_w(w_mm, sub: Substrato, dz_rel_max, de_rel_max, dh, der):
    """Tolerância de largura ±ΔW [mm] exigida para ΔZ0/Z0 e Δεeff/εeff máximos (eqs. 2.100-2.101);
    None se as tolerâncias de h e εr já consomem a especificação."""
    s = sensibilidades(w_mm, sub)
    rel = {"h": dh / sub.h_mm, "εr": der / sub.er}
    out = []
    for i, alvo in ((0, dz_rel_max), (1, de_rel_max)):
        resto = alvo - sum(abs(rel[k] * s[k][i]) for k in rel)
        if abs(s["W"][i]) < 1e-12:
            continue
        out.append(resto / abs(s["W"][i]) * w_mm if resto > 0 else None)
    if any(x is None for x in out):
        return None
    return min(out) if out else None


def monte_carlo(w_mm, sub: Substrato, dw, dh, der, n=2000, semente=1):
    """(Z0, εeff) para W, h e εr uniformes dentro das tolerâncias."""
    rng = np.random.default_rng(semente)
    w = w_mm + rng.uniform(-dw, dw, n)
    h = sub.h_mm + rng.uniform(-dh, dh, n)
    er = sub.er + rng.uniform(-der, der, n)
    z, e = np.empty(n), np.empty(n)
    for i in range(n):
        z[i], e[i] = _ms(w[i], replace(sub, h_mm=h[i], er=er[i]))
    return z, e


def pior_caso_acopladas(w_mm, s_mm, sub: Substrato, dw, ds, dh, der, modelo=None):
    """Variações máximas de Z0e, Z0o, Z0 = √(Z0e·Z0o) (eq. 4.125) e do acoplamento C (eq. 4.127)."""
    from . import linhas_acopladas as la
    modelo = modelo or la.MODELOS[0]

    def m(w=w_mm, s=s_mm, sb=sub):
        md = la.analisar(la.TIPOS[0], w, s, sb, modelo=modelo)
        return np.array([md.ze, md.zo, md.c])
    base = m()
    variacoes = [m(w=w_mm * (1 + 1e-4)), m(s=s_mm * (1 + 1e-4)),
                 m(sb=replace(sub, h_mm=sub.h_mm * (1 + 1e-4))), m(sb=replace(sub, er=sub.er * (1 + 1e-4)))]
    rel = [dw / w_mm, ds / s_mm, dh / sub.h_mm, der / sub.er]
    sens = [(v - base) / base / 1e-4 for v in variacoes]
    d = sum(abs(r * s) for r, s in zip(rel, sens))           # [ΔZe/Ze, ΔZo/Zo, ΔC/C]
    dz0 = 0.5 * (d[0] + d[1])
    c_db = -20 * np.log10(base[2])
    return dict(ze=base[0], zo=base[1], c=base[2], dze=d[0], dzo=d[1], dz0=dz0,
                dc_db=abs(20 * np.log10(1 + d[2])), c_db=c_db)


# ---------------------------------------------------------------------------
# Potência
# ---------------------------------------------------------------------------
@dataclass
class Potencia:
    dt_por_w: np.ndarray   # °C/W
    p_media: np.ndarray    # W
    p_pico: float          # W (ruptura do dielétrico)
    we: float              # mm (largura térmica equivalente)


def potencia(w_mm, sub: Substrato, f, k_termica, t_max=100.0, t_amb=25.0, e_ruptura_kv_cm=None):
    """Potência média (eqs. 2.84-2.95) e de pico (eq. 2.96) de uma microstrip.

    ΔT/W = 0,2303·h/K·(αc/We + αd/(2·Weff(f))), α em dB/m; We = η0·h/Za (Za em ar);
    Weff(f) = W + (Weff(0) - W)/(1 + (f/fp)²), fp = Z0/(2µ0h). Pico: V0²/(2Z0) com V0 = Erup·h.
    """
    f = np.atleast_1d(np.asarray(f, dtype=float))
    z0, ee = _ms(w_mm, sub)
    h = sub.h_mm * 1e-3
    we = ETA0 * h / (z0 * np.sqrt(ee))                       # Za = Z0·√εeff
    weff = w_mm * 1e-3 + (we - w_mm * 1e-3) / (1 + (f / (z0 / (2 * MU0 * h))) ** 2)
    ac = alfa_condutor(w_mm, sub, f)[:, 0] * NP_DB
    ad = alfa_dieletrico(w_mm, sub, f)[:, 0] * NP_DB
    dt = 0.2303 * h / k_termica * (ac / we + ad / (2 * weff))
    pp = (e_ruptura_kv_cm * 1e5 * h) ** 2 / (2 * z0) if e_ruptura_kv_cm else np.nan
    return Potencia(dt, (t_max - t_amb) / dt, float(pp), we * 1e3)


# ---------------------------------------------------------------------------
# Orçamento de erros (Tabela 2.5)
# ---------------------------------------------------------------------------
def orcamento(w_mm, sub: Substrato, f, dw, dh, der, vswr_conector=None, w2_mm=None):
    """Linhas (efeito, VSWR pela variação de Z0, Δεeff) como na Tabela 2.5."""
    from .descontinuidades import degrau
    z0, e0 = _ms(w_mm, sub)
    linhas = []
    if sub.t_um > 0:
        zt, et = _ms(w_mm, replace(sub, incluir_t=True))
        linhas.append((f"espessura t/h = {sub.t_um * 1e-3 / sub.h_mm:.4f}", max(z0, zt) / min(z0, zt), et - e0))
    zf, ef = (float(x[0, 0]) for x in dispersao(w_mm, sub, f))
    linhas.append((f"dispersão em {f / 1e9:g} GHz", max(z0, zf) / min(z0, zf), ef - e0))
    pc = pior_caso(w_mm, sub, dw, dh, der)
    linhas.append((f"tolerâncias ΔW ±{dw:g} mm, Δh ±{dh:g} mm, Δεr ±{der:g}", pc.vswr, pc.de_rel * e0))
    if w2_mm:
        dl = degrau(w_mm, w2_mm, sub)
        bl = 2 * np.pi * f * np.sqrt(e0) / 299792458.0 * dl * 1e-3
        g = abs(np.sin(bl)) / 2                                     # deslocamento de plano ≈ reatância
        linhas.append((f"degrau W2/W1 = {w2_mm / w_mm:.2f} (Δl = {dl:.3f} mm)", (1 + g) / (1 - g), None))
    if vswr_conector:
        linhas.append(("conector / transição", vswr_conector, None))
    vswr_total = np.prod([v for _, v, _ in linhas])
    return linhas, float(vswr_total)
