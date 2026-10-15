"""Perfis EBG definidos pelo usuário por uma expressão matemática em z.

A expressão pode descrever:
  - "K(z) [1/cm]": coeficiente de acoplamento; Z0(z) = Zref*exp(-2*integral(K dz))
                   (Martín 2015, eq. 2.100). O valor médio de K é removido para
                   que a linha termine em Zref.
  - "Z0(z) [Ω]":   impedância característica local.
  - "W(z) [mm]":   largura da trilha.

Variáveis: z (mm, de 0 a L), L (comprimento total, mm), Zref (Ω).
Funções auxiliares para projeto em múltiplas frequências:
  bragg(f_GHz) -> período [mm] cuja 1a banda (eq. 2.84) cai em f_GHz, usando o
                  εeff da linha de referência (aproximação: a banda real usa o
                  εeff médio do perfil e pode deslocar ~1-2%).
  janela(z)    -> janela de Hamming em [0, L] (reduz ondulação na banda passante).
  hann(z), ret(z, a, b) (1 entre a e b), degrau(z, a).

A avaliação é restrita (árvore sintática verificada): só números, operadores,
comparações e as funções listadas em FUNCOES.
"""
import ast

import numpy as np

from .microstrip import C0, Substrato, analise, sintese
from .modos_acoplados import Perfil

TIPOS_FUNCAO = ("K(z) [1/cm]", "Z0(z) [Ω]", "W(z) [mm]")

_FUNCOES_NUMPY = {
    "sin": np.sin, "cos": np.cos, "tan": np.tan, "arcsin": np.arcsin, "arccos": np.arccos,
    "arctan": np.arctan, "sinh": np.sinh, "cosh": np.cosh, "tanh": np.tanh, "exp": np.exp,
    "log": np.log, "log10": np.log10, "sqrt": np.sqrt, "abs": np.abs, "sign": np.sign,
    "floor": np.floor, "ceil": np.ceil, "mod": np.mod, "minimum": np.minimum, "maximum": np.maximum,
    "clip": np.clip, "where": np.where,
}
_CONSTANTES = {"pi": np.pi, "e": np.e}
_NOS_PERMITIDOS = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Call, ast.Name, ast.Load, ast.Constant,
                   ast.Compare, ast.BoolOp, ast.IfExp,
                   ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.Mod, ast.FloorDiv, ast.USub, ast.UAdd,
                   ast.Lt, ast.LtE, ast.Gt, ast.GtE, ast.Eq, ast.NotEq, ast.And, ast.Or, ast.Not)


def nomes_disponiveis():
    return sorted(set(_FUNCOES_NUMPY) | set(_CONSTANTES) | {"z", "L", "Zref", "bragg", "janela", "hann",
                                                           "ret", "degrau"})


def _validar(expr):
    try:
        arvore = ast.parse(expr.strip(), mode="eval")
    except SyntaxError as e:
        raise ValueError(f"Erro de sintaxe na expressão: {e.msg} (posição {e.offset}).") from None
    permitidos = set(nomes_disponiveis())
    for no in ast.walk(arvore):
        if not isinstance(no, _NOS_PERMITIDOS):
            raise ValueError(f"Construção não permitida na expressão: {type(no).__name__}.")
        if isinstance(no, ast.Name) and no.id not in permitidos:
            raise ValueError(f"Nome desconhecido: '{no.id}'. Disponíveis: {', '.join(nomes_disponiveis())}.")
        if isinstance(no, ast.Constant) and not isinstance(no.value, (int, float)):
            raise ValueError("Só são permitidas constantes numéricas.")
        if isinstance(no, ast.Call) and not isinstance(no.func, ast.Name):
            raise ValueError("Chamadas de função só por nome (ex.: sin(z)).")
    return compile(arvore, "<perfil>", "eval")


def avaliar(expr, z_mm, comprimento_mm, sub: Substrato, z_ref):
    """Avalia a expressão nos pontos z_mm e devolve um array do mesmo tamanho."""
    codigo = _validar(expr)
    _, ee_ref = analise(sintese(z_ref, sub), sub)
    ee_ref = float(ee_ref)

    def bragg(f_ghz):
        return C0 / (2 * np.asarray(f_ghz, dtype=float) * 1e9 * np.sqrt(ee_ref)) * 1e3

    def janela(z):
        return 0.54 - 0.46 * np.cos(2 * np.pi * np.asarray(z) / comprimento_mm)

    def hann(z):
        return 0.5 - 0.5 * np.cos(2 * np.pi * np.asarray(z) / comprimento_mm)

    def ret(z, a, b):
        z = np.asarray(z)
        return ((z >= a) & (z <= b)).astype(float)

    def degrau(z, a):
        return (np.asarray(z) >= a).astype(float)

    ambiente = {"__builtins__": {}, **_FUNCOES_NUMPY, **_CONSTANTES, "z": z_mm, "L": comprimento_mm,
                "Zref": z_ref, "bragg": bragg, "janela": janela, "hann": hann, "ret": ret, "degrau": degrau}
    with np.errstate(all="ignore"):
        try:
            val = eval(codigo, ambiente)  # noqa: S307 - árvore validada em _validar
        except (TypeError, ValueError, ZeroDivisionError) as e:
            raise ValueError(f"Erro ao avaliar a expressão: {e}") from None
    val = np.broadcast_to(np.asarray(val, dtype=float), z_mm.shape).copy()
    if not np.all(np.isfinite(val)):
        raise ValueError("A expressão gera valores infinitos ou indefinidos em algum ponto de z.")
    return val


def perfil(tipo_funcao, expr, comprimento_mm, passo_mm, sub: Substrato, z_ref=50.0):
    """Monta o Perfil (z em m, Z0 em ohm) a partir da expressão. Retorna (Perfil, avisos)."""
    if comprimento_mm <= 0 or passo_mm <= 0:
        raise ValueError("Comprimento e passo devem ser positivos.")
    n = int(round(comprimento_mm / passo_mm))
    if n < 10:
        raise ValueError("Passo grande demais para o comprimento (menos de 10 pontos).")
    if n > 20000:
        raise ValueError("Passo pequeno demais (mais de 20 000 pontos). Aumente o passo.")
    z_mm = np.linspace(0, comprimento_mm, n + 1)
    val = avaliar(expr, z_mm, comprimento_mm, sub, z_ref)
    avisos = []
    if tipo_funcao == TIPOS_FUNCAO[0]:
        k = val * 100.0                     # 1/cm -> 1/m
        k = k - np.trapezoid(k, z_mm) / comprimento_mm  # média nula: Z0(L) = Zref
        integral = np.concatenate([[0.0], np.cumsum(0.5 * (k[1:] + k[:-1]) * np.diff(z_mm * 1e-3))])
        z0 = z_ref * np.exp(-2 * integral)
    elif tipo_funcao == TIPOS_FUNCAO[1]:
        z0 = val
    elif tipo_funcao == TIPOS_FUNCAO[2]:
        if np.any(val <= 0):
            raise ValueError("A largura W(z) deve ser positiva em toda a linha.")
        z0, _ = analise(val, sub)
        z0 = np.asarray(z0, dtype=float)
    else:
        raise ValueError(f"Tipo de função desconhecido: {tipo_funcao}")
    if np.any(z0 <= 0):
        raise ValueError("A impedância Z0(z) resultante deve ser positiva.")
    if z0.max() / z0.min() > 20:
        avisos.append(f"Faixa de impedância muito ampla ({z0.min():.1f} a {z0.max():.1f} Ω).")
    for nome, zz in (("início", z0[0]), ("fim", z0[-1])):
        if abs(zz / z_ref - 1) > 0.02:
            avisos.append(f"Z0 no {nome} da linha = {zz:.1f} Ω difere de Zref = {z_ref:g} Ω "
                          "(haverá degrau na porta).")
    return Perfil("Personalizado", z_mm * 1e-3, z0, comprimento_mm * 1e-3, 1, {}), avisos


def estimativa_fourier(f, z, z0, eeff_med):
    """|S11| aproximado por -integral(K e^{-j2βz} dz) (eq. 2.95), válido p/ baixa refletividade."""
    k = -0.5 * np.diff(np.log(z0)) / np.diff(z)
    zm = 0.5 * (z[1:] + z[:-1])
    dz = np.diff(z)
    beta = 2 * np.pi * np.asarray(f)[:, None] * np.sqrt(eeff_med) / C0
    s11 = -np.sum(k[None, :] * np.exp(-2j * beta * zm[None, :]) * dz[None, :], axis=1)
    return np.minimum(np.abs(s11), 1.0)


def bandas_detectadas(f, s21, limiar_db=-10.0):
    """Faixas contíguas com |S21| abaixo do limiar: (f_início, f_fim, f_mín, S21_mín_dB)."""
    s21_db = 20 * np.log10(np.maximum(np.abs(s21), 1e-15))
    abaixo = s21_db < limiar_db
    bandas = []
    i = 0
    while i < len(f):
        if abaixo[i]:
            j = i
            while j + 1 < len(f) and abaixo[j + 1]:
                j += 1
            k = i + int(np.argmin(s21_db[i:j + 1]))
            bandas.append((f[i], f[j], f[k], float(s21_db[k])))
            i = j + 1
        else:
            i += 1
    return bandas


# Modelos prontos: (tipo de função, expressão, comprimento [mm], passo [mm], descrição)
MODELOS = {
    "Duas bandas: 3 e 5 GHz": (
        TIPOS_FUNCAO[0], "0.4*cos(2*pi*z/bragg(3)) + 0.4*cos(2*pi*z/bragg(5))", 150.0, 0.2,
        "Soma de duas senoides em K(z): cada uma gera uma banda (eq. 2.95)."),
    "Três bandas com janela Hamming: 3, 4 e 5 GHz": (
        TIPOS_FUNCAO[0], "janela(z)*(0.35*cos(2*pi*z/bragg(3)) + 0.35*cos(2*pi*z/bragg(4)) "
                         "+ 0.35*cos(2*pi*z/bragg(5)))", 260.0, 0.2,
        "Como na Fig. 2.20 do livro: a janela reduz a ondulação na banda passante."),
    "Supressão de 2f0 e 3f0 (f0 = 2,45 GHz)": (
        TIPOS_FUNCAO[0], "0.4*cos(2*pi*z/bragg(4.9)) + 0.4*cos(2*pi*z/bragg(7.35))", 130.0, 0.2,
        "EBG dupla para filtrar harmônicos de um circuito ativo (seção 2.5.1.4)."),
    "Chirp linear: banda larga 4 a 7 GHz": (
        TIPOS_FUNCAO[0], "1.2*janela(z)*cos(2*pi*(z/bragg(4) + (1/bragg(7) - 1/bragg(4))*z**2/(2*L)))",
        300.0, 0.1, "Período varia ao longo da linha: alarga a banda rejeitada (seção 2.5.1.1)."),
    "Livro Fig. 2.11 escrita como Z0(z)": (
        TIPOS_FUNCAO[1], "Zref*exp(2*0.064*10/pi*sin(2*pi*z/10))", 100.0, 0.1,
        "Mesmo perfil senoidal do livro (|K1| = 0,064 /mm, l = 10 mm)."),
    "Largura quase quadrada (bandas ímpares)": (
        TIPOS_FUNCAO[2], "1.19*exp(0.5*tanh(20*sin(2*pi*z/10)))", 100.0, 0.05,
        "W(z) em degraus suavizados: bandas em f e nas harmônicas ímpares; termina em W de 50 Ω."),
}
