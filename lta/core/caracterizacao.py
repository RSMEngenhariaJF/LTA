"""Caracterização de linhas e substratos a partir de medidas (parâmetros S, ressonâncias, VSWR).

Referência: R. Garg, I. Bahl e M. Bianchi, Microstrip Lines and Slotlines, 4ª ed., 2024:
  - εr a partir da fase de uma linha de 50 Ω (eqs. 1.89-1.90): aqui εeff vem da fase medida e εr é
    obtido invertendo o modelo de Hammerstad-Jensen (mais preciso que a fórmula F do livro);
  - εeff pelo método de dois ressonadores λ/2 de comprimentos L1 e L2 (eqs. 1.93-1.95), que
    cancela os efeitos de ponta e de acoplamento;
  - α a partir do Q0 de um ressonador (eq. 1.96) e pelo método de duas linhas (comparação, que
    elimina transições e conectores);
  - Z0 por VSWR máximo e mínimo (eq. 1.91);
  - modelo RLGC da linha a partir de S (eqs. 2.143-2.148).
"""
from dataclasses import dataclass

import numpy as np

from . import abcd
from .microstrip import C0, NP_DB, Substrato, analise


@dataclass
class Linha:
    f: np.ndarray
    zc: np.ndarray        # impedância característica (complexa) [Ω]
    gama: np.ndarray      # γ = α + jβ [1/m]

    @property
    def eeff(self):
        return (self.gama.imag * C0 / (2 * np.pi * self.f)) ** 2

    @property
    def alfa_db_cm(self):
        return self.gama.real * NP_DB / 100

    def rlgc(self):
        """(R [Ω/m], L [H/m], G [S/m], C [F/m]) por unidade de comprimento (eqs. 2.143, 2.148)."""
        w = 2 * np.pi * self.f
        zs, yp = self.gama * self.zc, self.gama / self.zc
        return zs.real, zs.imag / w, yp.real, yp.imag / w


def _desenrolar_beta(exp_menos_gl, comp_m, f):
    """β a partir de e^(-γl): fase desenrolada e múltiplo de 2π escolhido para que a reta
    ajustada passe pela origem (β ∝ f numa linha quase-TEM)."""
    fase = -np.unwrap(np.angle(exp_menos_gl))
    if len(f) > 2:
        a = np.polyfit(f, fase, 1)[1]
        fase = fase - 2 * np.pi * np.round(a / (2 * np.pi))
    elif fase[0] < 0:
        fase = fase + 2 * np.pi * np.ceil(-fase[0] / (2 * np.pi))
    return fase / comp_m


def linha_de_s(f, s, comp_mm, z_ref=50.0):
    """Zc e γ de uma linha uniforme medida (simétrica e recíproca), eqs. 2.147-2.148:
    e^(-γl) = [(1 - S11² + S21²)/(2S21) ± K]⁻¹, K = √(((S11² - S21² + 1)² - (2S11)²)/(2S21)²),
    Zc = Zref·√(((1 + S11)² - S21²)/((1 - S11)² - S21²))."""
    f = np.asarray(f, dtype=float)
    s11, s21 = s[:, 0, 0], s[:, 1, 0]
    k = np.sqrt(((s11**2 - s21**2 + 1) ** 2 - (2 * s11) ** 2) / (2 * s21) ** 2)
    base = (1 - s11**2 + s21**2) / (2 * s21)
    e1, e2 = 1 / (base + k), 1 / (base - k)
    egl = np.where(np.abs(e1) <= np.abs(e2), e1, e2)
    zc = z_ref * np.sqrt(((1 + s11) ** 2 - s21**2) / ((1 - s11) ** 2 - s21**2))
    zc = np.where(zc.real < 0, -zc, zc)
    l = comp_mm * 1e-3
    alfa = -np.log(np.abs(egl)) / l
    return Linha(f, zc, alfa + 1j * _desenrolar_beta(egl, l, f))


def duas_linhas(f, s1, s2, comp1_mm, comp2_mm, z_ref=50.0):
    """γ pelo método de duas linhas: M1⁻¹·M2 = B⁻¹·L(ΔL)·B tem autovalores e^(±γΔL), sem depender
    das transições e conectores (iguais nas duas medidas)."""
    f = np.asarray(f, dtype=float)
    dl = (comp2_mm - comp1_mm) * 1e-3
    if dl <= 0:
        raise ValueError("A segunda linha deve ser mais longa que a primeira.")
    m = np.linalg.inv(abcd.s_para_abcd(s1, z_ref)) @ abcd.s_para_abcd(s2, z_ref)
    tr = m[:, 0, 0] + m[:, 1, 1]
    det = m[:, 0, 0] * m[:, 1, 1] - m[:, 0, 1] * m[:, 1, 0]
    raiz = np.sqrt(tr**2 / 4 - det + 0j)
    l1, l2 = tr / 2 + raiz, tr / 2 - raiz
    egl = np.where(np.abs(l1) <= np.abs(l2), l1, l2)        # e^(-γΔL)
    alfa = -np.log(np.abs(egl)) / dl
    return alfa + 1j * _desenrolar_beta(egl, dl, f)


def er_de_eeff(eeff, w_mm, h_mm, er_max=30.0):
    """εr do substrato que dá o εeff medido para a largura W (inverte Hammerstad-Jensen)."""
    lo, hi = 1.0001, er_max
    if not float(analise(w_mm, Substrato(lo, h_mm))[1]) < eeff < float(analise(w_mm, Substrato(hi, h_mm))[1]):
        raise ValueError(f"εeff = {eeff:.3f} fora do alcance do modelo para W/h = {w_mm / h_mm:.3f}.")
    for _ in range(80):
        m = 0.5 * (lo + hi)
        if float(analise(w_mm, Substrato(m, h_mm))[1]) < eeff:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def er_de_fase(f, fase_graus, comp_mm, w_mm, h_mm):
    """εr pela fase de transmissão φ (atraso, positivo e desenrolado) de uma linha de comprimento ℓ:
    εeff = (φ·c/(360·f·ℓ))² (eq. 1.90) e depois a inversão do modelo."""
    eeff = (fase_graus * C0 / (360 * f * comp_mm * 1e-3)) ** 2
    return er_de_eeff(eeff, w_mm, h_mm), eeff


def dois_ressonadores(f1, f2, comp1_mm, comp2_mm, n=1):
    """(εeff, ℓe [mm]) pelo método dos dois ressonadores λ/2 (eqs. 1.94-1.95), modo n."""
    if np.isclose(f1, f2):
        raise ValueError("As ressonâncias precisam ser diferentes.")
    eeff = (n * C0 * (f1 - f2) / (2 * f1 * f2 * (comp2_mm - comp1_mm) * 1e-3)) ** 2
    le = (f2 * comp2_mm - f1 * comp1_mm) / (f1 - f2)
    return float(eeff), float(le)


def alfa_de_q0(q0, lambda_g_mm):
    """α [dB/cm] = 27,3/(Q0·λg[cm]) (eq. 1.96; β/(2Q0) em dB)."""
    return 27.3 / (q0 * lambda_g_mm / 10)


def z0_de_vswr(z_ref, vswr_max, vswr_min, alta=True):
    """Z0 da linha medida em sistema de Zref: Zref·(VSWRmax·VSWRmin)^(±1/2) (eq. 1.91)."""
    return z_ref * (vswr_max * vswr_min) ** (0.5 if alta else -0.5)


def sintetico(f, z_linha, eeff, alfa_np_m, comp_mm, transicao=None, z_ref=50.0):
    """S de uma linha (com transições, se dadas) para testar os métodos sem medidas reais."""
    from .conectores import abcd_transicao
    f = np.asarray(f, dtype=float)
    th = 2 * np.pi * f * np.sqrt(eeff) / C0 * comp_mm * 1e-3 - 1j * alfa_np_m * comp_mm * 1e-3
    m = abcd.linha(z_linha, th)
    if transicao is not None:
        m = abcd.cascata(abcd_transicao(f, transicao, True), m, abcd_transicao(f, transicao, False))
    return abcd.abcd_para_s(m, z_ref)


# ---------------------------------------------------------------------------
# Calibração TRL e renormalização (Steer, Microwave and RF Design, seções 6.10-6.11)
# ---------------------------------------------------------------------------
def s_para_t(s):
    """Matriz de cascata T com [a1; b1] = T·[b2; a2]: T = (1/S21)·[[1, -S22], [S11, -Δ]]."""
    s11, s12, s21, s22 = s[..., 0, 0], s[..., 0, 1], s[..., 1, 0], s[..., 1, 1]
    t = np.empty(s.shape, dtype=complex)
    t[..., 0, 0] = 1 / s21
    t[..., 0, 1] = -s22 / s21
    t[..., 1, 0] = s11 / s21
    t[..., 1, 1] = -(s11 * s22 - s12 * s21) / s21
    return t


def t_para_s(t):
    t11, t12, t21, t22 = t[..., 0, 0], t[..., 0, 1], t[..., 1, 0], t[..., 1, 1]
    s = np.empty(t.shape, dtype=complex)
    s[..., 1, 0] = 1 / t11
    s[..., 1, 1] = -t12 / t11
    s[..., 0, 0] = t21 / t11
    s[..., 0, 1] = t22 - t21 * t12 / t11
    return s


@dataclass
class ResultadoTRL:
    s_dut: np.ndarray       # S corrigido do dispositivo (planos no fim das caixas de erro)
    gama: np.ndarray        # γ da linha [1/m]
    gama_r: np.ndarray      # coeficiente de reflexão do padrão Reflect
    caixa_a: np.ndarray     # S da caixa de erro da porta 1
    caixa_b: np.ndarray     # S da caixa de erro da porta 2


def trl(f, s_thru, s_line, gamma_m1, gamma_m2, s_dut_medido, comp_linha_mm, eeff_estimado, reflect_aberto=True):
    """Calibração TRL (Engen e Hoer, 1979; Steer, seção 6.11.5) com thru de comprimento zero.

    Cascatas T: Thru = X·Ȳ, Line = X·L·Ȳ com L = diag(e^γl, e^-γl); M = Line·Thru⁻¹ = X·L·X⁻¹
    dá as colunas de X (autovetores) a menos de escalas; a escolha do autovalor usa a fase
    estimada βl. O Reflect (desconhecido, igual nas duas portas, medido como Γm1 e Γm2) fixa a
    razão das escalas: u = ΓR·k2/k1 vem da porta 1, v = ΓR·k1/k2 da porta 2, ΓR = ±√(uv) com o
    sinal pelo tipo do padrão (aberto ≈ +1, curto ≈ -1). Entradas: S (Nf, 2, 2) e Γm (Nf,).
    """
    f = np.asarray(f, dtype=float)
    tt, tl, tm = s_para_t(s_thru), s_para_t(s_line), s_para_t(s_dut_medido)
    m = tl @ np.linalg.inv(tt)
    lam, vec = np.linalg.eig(m)
    alvo = np.exp(1j * 2 * np.pi * f * np.sqrt(eeff_estimado) / C0 * comp_linha_mm * 1e-3)
    d0 = np.abs(np.angle(lam[:, 0] / alvo))
    d1 = np.abs(np.angle(lam[:, 1] / alvo))
    troca = d1 < d0
    e = vec.copy()
    e[troca] = e[troca][:, :, ::-1]
    lam1 = np.where(troca, lam[:, 1], lam[:, 0])
    fmat = np.linalg.inv(e) @ tt
    g1, g2 = np.asarray(gamma_m1), np.asarray(gamma_m2)
    u = (e[:, 1, 0] - g1 * e[:, 0, 0]) / (g1 * e[:, 0, 1] - e[:, 1, 1])
    v = (g2 * fmat[:, 0, 0] + fmat[:, 0, 1]) / (fmat[:, 1, 1] + g2 * fmat[:, 1, 0])
    gr = np.sqrt(u * v)
    gr = np.where((gr.real > 0) != reflect_aberto, -gr, gr)
    r = gr / u                                                        # k1/k2
    d = np.zeros_like(e)
    d[:, 0, 0], d[:, 1, 1] = r, 1.0
    x = e @ d
    y = np.linalg.inv(x) @ tt
    t_dut = np.linalg.inv(x) @ tm @ np.linalg.inv(y)
    l_m = comp_linha_mm * 1e-3
    gama = np.log(lam1) / l_m
    gama = gama.real + 1j * (np.unwrap(np.angle(lam1)) / l_m)
    # caixas de erro em S: A (porta 1 → DUT) = X; B (DUT → porta 2) = Ȳ
    return ResultadoTRL(t_para_s(t_dut), gama, gr, t_para_s(x), t_para_s(y))


def renormalizar(s, z_antiga, z_nova):
    """S (N portas, referências reais) de z_antiga para z_nova (escalares ou por porta), via Z
    (Steer, seção 6.10): Z = R^½(I + S)(I - S)⁻¹R^½ e S' = R'^-½(Z - R')(Z + R')⁻¹R'^½."""
    n = s.shape[-1]
    eye = np.eye(n)
    r0 = np.broadcast_to(np.asarray(z_antiga, dtype=float), (n,))
    r1 = np.broadcast_to(np.asarray(z_nova, dtype=float), (n,))
    a0 = np.sqrt(r0)
    z = a0[:, None] * ((eye + s) @ np.linalg.inv(eye - s)) * a0[None, :]
    rn = np.diag(r1)
    a1 = np.sqrt(r1)
    return ((z - rn) @ np.linalg.inv(z + rn)) / a1[:, None] * a1[None, :]
