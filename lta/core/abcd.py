"""Método da matriz de transferência (ABCD) para estruturas periódicas.

Referência: F. Martín, Artificial Transmission Lines for RF and Microwave
Applications, Wiley, 2015, seção 2.3 (eqs. 2.15-2.42).

Todas as matrizes têm forma (Nf, 2, 2): uma matriz 2x2 por frequência.
"""
import numpy as np


def series(z):
    """Impedância série z (array por frequência)."""
    z = np.asarray(z, dtype=complex)
    m = np.zeros(z.shape + (2, 2), dtype=complex)
    m[..., 0, 0] = 1
    m[..., 0, 1] = z
    m[..., 1, 1] = 1
    return m


def shunt(y):
    """Admitância em derivação y (array por frequência)."""
    y = np.asarray(y, dtype=complex)
    m = np.zeros(y.shape + (2, 2), dtype=complex)
    m[..., 0, 0] = 1
    m[..., 1, 0] = y
    m[..., 1, 1] = 1
    return m


def linha(z0, theta):
    """Trecho de linha sem perdas: impedância z0 e comprimento elétrico theta [rad]."""
    theta = np.asarray(theta, dtype=complex)
    z0 = np.broadcast_to(np.asarray(z0, dtype=complex), theta.shape)
    m = np.empty(theta.shape + (2, 2), dtype=complex)
    m[..., 0, 0] = np.cos(theta)
    m[..., 0, 1] = 1j * z0 * np.sin(theta)
    m[..., 1, 0] = 1j * np.sin(theta) / z0
    m[..., 1, 1] = np.cos(theta)
    return m


def cascata(*mats):
    out = mats[0]
    for m in mats[1:]:
        out = out @ m
    return out


def potencia(m, n):
    """Cascata de n células idênticas."""
    return np.linalg.matrix_power(m, int(n))


def dispersao(m):
    """Relação de dispersão cosh(gamma*l) = (A+D)/2 (eq. 2.21).

    Retorna (alpha*l [Np], beta*l [rad]) com alpha >= 0 e 0 <= beta <= pi.
    """
    x = 0.5 * (m[..., 0, 0] + m[..., 1, 1])
    g = np.arccosh(x.astype(complex))
    return np.abs(g.real), np.abs(g.imag)


def impedancia_bloch(m):
    """Impedâncias de Bloch ZB+ e ZB- (eq. 2.26).

    ZB+ é a solução com parte real >= 0 (onda progressiva). Em células
    simétricas e sem perdas, é real na banda passante e imaginária na banda
    rejeitada.
    """
    a, b, d = m[..., 0, 0], m[..., 0, 1], m[..., 1, 1]
    raiz = np.sqrt((a + d) ** 2 - 4 + 0j)
    z1 = -2 * b / (a - d - raiz)
    z2 = -2 * b / (a - d + raiz)
    zp = np.where(z1.real >= z2.real, z1, z2)
    zm = np.where(z1.real >= z2.real, z2, z1)
    return zp, zm


def impedancia_iterativa(m):
    """Impedância iterativa ZI (eq. 2.40); igual a ZB para célula recíproca."""
    a, c, d = m[..., 0, 0], m[..., 1, 0], m[..., 1, 1]
    raiz = np.sqrt((a + d) ** 2 - 4 + 0j)
    z1 = (a - d + raiz) / (2 * c)
    z2 = (a - d - raiz) / (2 * c)
    return np.where(z1.real >= z2.real, z1, z2)


def abcd_para_s(m, z0=50.0):
    """Conversão ABCD -> S com impedância de referência z0 nas duas portas."""
    a, b, c, d = m[..., 0, 0], m[..., 0, 1], m[..., 1, 0], m[..., 1, 1]
    den = a + b / z0 + c * z0 + d
    s = np.empty(m.shape, dtype=complex)
    s[..., 0, 0] = (a + b / z0 - c * z0 - d) / den
    s[..., 0, 1] = 2 * (a * d - b * c) / den
    s[..., 1, 0] = 2 / den
    s[..., 1, 1] = (-a + b / z0 - c * z0 + d) / den
    return s


def abcd_para_s2(m, z01, z02):
    """ABCD -> S com impedâncias de referência reais diferentes nas portas 1 (z01) e 2 (z02)."""
    a, b, c, d = m[..., 0, 0], m[..., 0, 1], m[..., 1, 0], m[..., 1, 1]
    den = a * z02 + b + c * z01 * z02 + d * z01
    r = 2 * np.sqrt(z01 * z02)
    s = np.empty(m.shape, dtype=complex)
    s[..., 0, 0] = (a * z02 + b - c * z01 * z02 - d * z01) / den
    s[..., 0, 1] = r * (a * d - b * c) / den
    s[..., 1, 0] = r / den
    s[..., 1, 1] = (-a * z02 + b - c * z01 * z02 + d * z01) / den
    return s


def s_para_abcd(s, z0=50.0):
    """Conversão S -> ABCD (inversa de abcd_para_s), referência z0 nas duas portas."""
    s11, s12, s21, s22 = s[..., 0, 0], s[..., 0, 1], s[..., 1, 0], s[..., 1, 1]
    m = np.empty(s.shape, dtype=complex)
    m[..., 0, 0] = ((1 + s11) * (1 - s22) + s12 * s21) / (2 * s21)
    m[..., 0, 1] = z0 * ((1 + s11) * (1 + s22) - s12 * s21) / (2 * s21)
    m[..., 1, 0] = ((1 - s11) * (1 - s22) - s12 * s21) / (2 * s21 * z0)
    m[..., 1, 1] = ((1 - s11) * (1 + s22) + s12 * s21) / (2 * s21)
    return m


def db(x):
    return 20 * np.log10(np.maximum(np.abs(x), 1e-15))


def atraso_grupo(f, s21):
    """Atraso de grupo τg = -dφ/dω [s] a partir da fase de S21 (φ desenrolada)."""
    fase = np.unwrap(np.angle(s21))
    return -np.gradient(fase, 2 * np.pi * np.asarray(f, dtype=float))


# ---------------------------------------------------------------------------
# Célula genérica de elementos concentrados
# ---------------------------------------------------------------------------
TOPOLOGIAS = ("T", "Pi", "L", "L invertida", "Treliça (X)", "T em ponte")


def ramos(f, ls=0.0, cs=0.0, lp=0.0, cp=0.0):
    """Impedância do ramo série e admitância do ramo paralelo.

    Ramo série: Ls em série com Cs. Ramo paralelo: Cp em paralelo com Lp.
    Valor 0 significa elemento ausente (Cs=0 -> curto; Lp=0 -> aberto).
    Unidades SI (H, F, Hz).
    """
    w = 2 * np.pi * np.asarray(f, dtype=float)
    zs = 1j * w * ls + (1 / (1j * w * cs) if cs > 0 else 0)
    yp = 1j * w * cp + (1 / (1j * w * lp) if lp > 0 else 0)
    return np.asarray(zs, dtype=complex) * np.ones_like(w), np.asarray(yp, dtype=complex) * np.ones_like(w)


def celula(f, topologia, ls=0.0, cs=0.0, lp=0.0, cp=0.0):
    zs, yp = ramos(f, ls, cs, lp, cp)
    if topologia == "T":
        return cascata(series(zs / 2), shunt(yp), series(zs / 2))
    if topologia == "Pi":
        return cascata(shunt(yp / 2), series(zs), shunt(yp / 2))
    if topologia == "L":
        return cascata(series(zs), shunt(yp))
    if topologia == "L invertida":
        return cascata(shunt(yp), series(zs))
    if topologia == "Treliça (X)":
        return trelica(zs, 1 / yp)
    if topologia == "T em ponte":
        return t_em_ponte(zs, 1 / yp)
    raise ValueError(f"Topologia desconhecida: {topologia}")


def trelica(za, zb):
    """Rede em treliça (X) simétrica: braços série za e diagonais zb.

    cos(beta*l) = (zb + za)/(zb - za) e ZB = sqrt(za*zb) (Martín 2015, eqs. 6.13-6.14).
    Se za*zb é real e constante (ex.: L em série e C nas diagonais), ZB não depende
    da frequência: rede passa-tudo de resistência constante.
    """
    za = np.asarray(za, dtype=complex)
    zb = np.asarray(zb, dtype=complex)
    d = zb - za
    m = np.empty(za.shape + (2, 2), dtype=complex)
    m[..., 0, 0] = m[..., 1, 1] = (zb + za) / d
    m[..., 0, 1] = 2 * za * zb / d
    m[..., 1, 0] = 2 / d
    return m


def t_em_ponte(z1, z2):
    """T em ponte: braços série z1, derivação z2 e ponte z3 = 4*z2 (condição 6.27).

    Com essa condição, ZB = sqrt(2*z1*z2) (Martín 2015, eqs. 6.25-6.29, pp. 437-439).
    Calculado somando a admitância da ponte à matriz Y da célula T.
    """
    z1 = np.asarray(z1, dtype=complex)
    z2 = np.asarray(z2, dtype=complex)
    a = 1 + z1 / z2
    b = 2 * z1 + z1**2 / z2
    y3 = 1 / (4 * z2)
    y11 = a / b + y3            # célula T simétrica: Y11 = Y22 = A/B, Y12 = Y21 = -1/B
    y21 = -1 / b - y3
    det = y11**2 - y21**2
    m = np.empty(z1.shape + (2, 2), dtype=complex)
    m[..., 0, 0] = m[..., 1, 1] = -y11 / y21
    m[..., 0, 1] = -1 / y21
    m[..., 1, 0] = -det / y21
    return m


def banda_util(f, s, s11_max_db=-15.0, s21_min_db=-1.0):
    """Faixas em que |S11| <= s11_max_db e |S21| >= s21_min_db (simultaneamente).

    Retorna lista de (f_inicial, f_final), da mais larga para a mais estreita.
    """
    ok = (db(s[:, 0, 0]) <= s11_max_db) & (db(s[:, 1, 0]) >= s21_min_db)
    faixas = []
    i = 0
    while i < len(f):
        if ok[i]:
            j = i
            while j + 1 < len(f) and ok[j + 1]:
                j += 1
            faixas.append((f[i], f[j]))
            i = j + 1
        else:
            i += 1
    return sorted(faixas, key=lambda x: x[1] - x[0], reverse=True)


def descrever_banda(f, s, s11_max_db=-15.0, s21_min_db=-1.0, max_faixas=3):
    """Texto com as faixas úteis (em ordem de frequência), destacando a mais larga."""
    faixas = banda_util(f, s, s11_max_db, s21_min_db)
    criterio = f"Banda útil (|S11| ≤ {s11_max_db:g} dB e |S21| ≥ {s21_min_db:g} dB)"
    if not faixas:
        return f"{criterio}: nenhuma"
    maior = faixas[0]
    partes = []
    for a, b in sorted(faixas[:max_faixas]):
        txt = f"{a / 1e9:.3f} a {b / 1e9:.3f} GHz"
        if b >= f[-1]:
            txt += " (até o fim da faixa calculada)"
        if (a, b) == maior and len(faixas) > 1:
            txt += " [maior]"
        partes.append(txt)
    a, b = maior
    razao = f"; maior faixa com razão {b / a:.1f}:1" if a > 0 else ""
    extra = f"; +{len(faixas) - max_faixas} faixas estreitas" if len(faixas) > max_faixas else ""
    return f"{criterio}: " + "; ".join(partes) + razao + extra
