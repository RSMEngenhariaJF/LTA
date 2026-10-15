"""Figuras didáticas e números dos exemplos do manual, calculados com o núcleo do LTA."""
import math
import sys
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure

RAIZ = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RAIZ))
IMG = Path(__file__).resolve().parent / "img"

from lta.core import abcd  # noqa: E402
from lta.core import carga_capacitiva as cc  # noqa: E402
from lta.core import modos_acoplados as cm  # noqa: E402
from lta.core.microstrip import C0, Substrato, analise, sintese  # noqa: E402

RO3010 = Substrato(10.2, 1.27, 35, 0.0022, "Rogers RO3010")
RO4003 = Substrato(3.55, 0.813, 35, 0.0027, "Rogers RO4003C")
FR4 = Substrato(4.4, 1.6, 35, 0.02, "FR-4")


def _salvar(fig, nome):
    IMG.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMG / f"fig_{nome}.png", dpi=200)


def fig_dispersao_conceito():
    """Linha levemente perturbada: bandas proibidas em beta*l = n*pi (Bragg)."""
    p = cc.ProjetoCarga(L=5e-9, C=2e-12, Cls=0.25e-12)
    f = np.linspace(1e6, 3 * p.fs * 1.02, 4000)
    al, bl = abcd.dispersao(cc.celula_distribuida(f, p))
    kl = 2 * np.pi * f * np.sqrt(p.L * p.C)
    dobrada = np.abs(((kl + np.pi) % (2 * np.pi)) - np.pi)
    fig = Figure(figsize=(7.5, 3.6), layout="constrained")
    ax1, ax2 = fig.subplots(1, 2)
    ax1.plot(dobrada / np.pi, f / 1e9, ":", color="gray", label="Linha sem perturbação (rebatida)")
    ax1.plot(bl / np.pi, f / 1e9, color="C0", label="Linha periódica")
    for n in (1, 2, 3):
        ax1.axhline(n * p.fs / 1e9, color="C3", lw=0.6, ls="--")
    ax1.set_xlabel("βl/π  (zona de Brillouin reduzida)")
    ax1.set_ylabel("f [GHz]")
    ax1.set_title("(a) Diagrama de dispersão", fontsize=10)
    ax1.legend(fontsize=7, loc="lower right")
    ax2.plot(al, f / 1e9, color="C1")
    ax2.set_xlabel("αl [Np por célula]")
    ax2.set_title("(b) Atenuação: bandas proibidas", fontsize=10)
    for a in (ax1, ax2):
        a.grid(True, alpha=0.3)
    _salvar(fig, "dispersao_conceito")
    return dict(fs=p.fs)


def fig_impacto_cls():
    """Efeito da capacitância de carga sobre dispersão e impedância de Bloch."""
    L, C = 3.0e-9, 0.3e-12
    f = np.linspace(1e7, 25e9, 3000)
    fig = Figure(figsize=(7.5, 3.6), layout="constrained")
    ax1, ax2 = fig.subplots(1, 2)
    tabela = []
    for i, r in enumerate((0.0, 1.0, 3.0, 9.0)):
        p = cc.ProjetoCarga(L, C, max(r * C, 1e-18))
        al, bl = abcd.dispersao(cc.celula_distribuida(f, p))
        zb, _ = abcd.impedancia_bloch(cc.celula_distribuida(f, p))
        primeira = f < (p.fs if r == 0 else p.fc * 1.3)
        ax1.plot(bl[primeira] / np.pi, f[primeira] / 1e9, color=f"C{i}", label=f"Cls/C = {r:g}")
        passa = (al < 1e-6) & (f < p.fs)
        ax2.plot(f[passa] / 1e9, zb[passa].real, color=f"C{i}", label=f"Cls/C = {r:g}")
        tabela.append((r, p.zb, p.fc / 1e9 if r else float("nan"), p.fator_slow_wave))
    ax1.set_xlabel("βl/π")
    ax1.set_ylabel("f [GHz]")
    ax1.set_title("(a) 1ª banda: inclinação = velocidade", fontsize=10)
    ax2.set_xlabel("f [GHz]")
    ax2.set_ylabel("Re ZB [Ω]")
    ax2.set_title("(b) Impedância de Bloch", fontsize=10)
    for a in (ax1, ax2):
        a.grid(True, alpha=0.3)
        a.legend(fontsize=7)
    _salvar(fig, "impacto_cls")
    return dict(L=L, C=C, zo=math.sqrt(L / C), tabela=tabela)


def fig_impacto_ebg():
    """Efeito de |K| e do número de períodos sobre a banda rejeitada (CMT)."""
    f = np.linspace(3e9, 9e9, 1200)
    fig = Figure(figsize=(7.5, 3.6), layout="constrained")
    ax1, ax2 = fig.subplots(1, 2, sharey=True)
    res_k, res_n = [], []
    for i, k in enumerate((0.2, 0.4, 0.64, 1.0)):
        p = cm.perfil("Senoidal", k * 100, 0.01, 10)
        w, ee = cm.larguras(p, RO3010)
        _, s21 = cm.resolver_cmt(f, p.z, p.z0, ee)
        b = cm.bandas_analiticas(p, cm.eeff_medio(p.z, ee))[0]
        ax1.plot(f / 1e9, abcd.db(s21), color=f"C{i}", label=f"|K| = {k:g} cm⁻¹")
        res_k.append((k, 20 * math.log10(b.s21_min), b.banda / 1e9, p.z0.min(), p.z0.max()))
    for i, n in enumerate((4, 8, 16)):
        p = cm.perfil("Senoidal", 40.0, 0.01, n)
        w, ee = cm.larguras(p, RO3010)
        _, s21 = cm.resolver_cmt(f, p.z, p.z0, ee)
        b = cm.bandas_analiticas(p, cm.eeff_medio(p.z, ee))[0]
        ax2.plot(f / 1e9, abcd.db(s21), color=f"C{i}", label=f"N = {n} períodos")
        res_n.append((n, 20 * math.log10(b.s21_min), b.banda / 1e9))
    ax1.set_title("(a) Amplitude |K| (N = 10)", fontsize=10)
    ax2.set_title("(b) Número de períodos (|K| = 0,4 cm⁻¹)", fontsize=10)
    ax1.set_ylabel("|S21| [dB]")
    for a in (ax1, ax2):
        a.set_xlabel("f [GHz]")
        a.set_ylim(-70, 3)
        a.grid(True, alpha=0.3)
        a.legend(fontsize=7)
    _salvar(fig, "impacto_ebg")
    return dict(k=res_k, n=res_n)


def fig_projeto_periodo():
    """Ábaco: frequência de Bragg (1ª banda) em função do período, linha de 50 Ω."""
    l = np.linspace(5, 40, 300)
    fig = Figure(figsize=(6.0, 3.4), layout="constrained")
    ax = fig.subplots()
    dados = []
    for i, s in enumerate((FR4, RO4003, RO3010)):
        w = sintese(50, s)
        _, ee = analise(w, s)
        ax.plot(l, cm.frequencia_bragg(l * 1e-3, ee) / 1e9, color=f"C{i}",
                label=f"{s.nome} (εr={s.er:g}, h={s.h_mm:g} mm): εeff={float(ee):.2f}")
        dados.append((s.nome, float(ee)))
    ax.set_xlabel("Período l [mm]")
    ax.set_ylabel("f Bragg [GHz]")
    ax.set_ylim(0, 12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)
    _salvar(fig, "projeto_periodo")
    return dados


# ---------------------------------------------------------------------------
# Exemplos de aplicação
# ---------------------------------------------------------------------------
def _k_da_razao(razao, l_m):
    """|K| [1/cm] do perfil senoidal com Zmax/Zref = razao (eq. 2.102)."""
    return math.pi * math.log(razao) / (2 * l_m) / 100


def _periodo_para(f_alvo, razao, n, sub):
    """Itera o período até a 1a banda (eq. 2.84) cair em f_alvo, com Zmax/Zref fixo."""
    w = sintese(50, sub)
    _, ee = analise(w, sub)
    l = C0 / (2 * f_alvo * math.sqrt(ee))
    for _ in range(6):
        p = cm.perfil("Senoidal", _k_da_razao(razao, l) * 100, l, n)
        w, ee_z = cm.larguras(p, sub)
        l = C0 / (2 * f_alvo * math.sqrt(cm.eeff_medio(p.z, ee_z)))
    p = cm.perfil("Senoidal", _k_da_razao(razao, l) * 100, l, n)
    w, ee_z = cm.larguras(p, sub)
    return p, w, ee_z, cm.eeff_medio(p.z, ee_z)


def exemplo_refletor():
    """Refletor a 2,45 GHz em FR-4, |S21| <= -30 dB no centro, Zmax/Zref = 1,5."""
    f0, razao, alvo_db = 2.45e9, 1.5, -30.0
    kl_min = math.acosh(10 ** (-alvo_db / 20))
    kl_periodo = math.pi * math.log(razao) / 2          # |K|*l, independe da frequência
    n = math.ceil(kl_min / kl_periodo)
    w50 = sintese(50, FR4)
    _, ee50 = analise(w50, FR4)
    l0 = C0 / (2 * f0 * math.sqrt(ee50))
    p, w, ee, eem = _periodo_para(f0, razao, n, FR4)
    b = cm.bandas_analiticas(p, eem)[0]
    f = np.linspace(1e9, 4e9, 1500)
    _, s21 = cm.resolver_cmt(f, p.z, p.z0, ee)
    fig = Figure(figsize=(6.0, 3.0), layout="constrained")
    ax = fig.subplots()
    ax.plot(f / 1e9, abcd.db(s21))
    ax.axvline(f0 / 1e9, color="C3", ls="--", lw=0.8)
    ax.axhline(alvo_db, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("f [GHz]")
    ax.set_ylabel("|S21| [dB]")
    ax.grid(True, alpha=0.3)
    _salvar(fig, "ex_refletor")
    return dict(f0=f0, razao=razao, alvo=alvo_db, kl_min=kl_min, kl_periodo=kl_periodo, w50=w50,
                ee50=float(ee50), l0=l0 * 1e3, n=n, k=_k_da_razao(razao, p.periodo),
                l=p.periodo * 1e3, comp=p.comprimento * 1e3, eem=eem, fmax=b.f_max / 1e9,
                s21=20 * math.log10(b.s21_min), bw=b.banda / 1e9, zmin=p.z0.min(), zmax=p.z0.max(),
                wmin=w.min(), wmax=w.max(), s21_num=float(abcd.db(s21).min()))


def exemplo_harmonico():
    """Supressão do 2o harmônico (4,9 GHz) de um oscilador a 2,45 GHz em RO4003C."""
    f0, razao, n = 2.45e9, 1.5, 8
    p, w, ee, eem = _periodo_para(2 * f0, razao, n, RO4003)
    b = cm.bandas_analiticas(p, eem)[0]
    f = np.linspace(0.5e9, 8e9, 2000)
    s11, s21 = cm.resolver_cmt(f, p.z, p.z0, ee)
    i0 = int(np.argmin(abs(f - f0)))
    i2 = int(np.argmin(abs(f - 2 * f0)))
    fig = Figure(figsize=(6.0, 3.0), layout="constrained")
    ax = fig.subplots()
    ax.plot(f / 1e9, abcd.db(s21), label="|S21|")
    ax.plot(f / 1e9, abcd.db(s11), label="|S11|", alpha=0.7)
    for x, t in ((f0, "f0"), (2 * f0, "2f0"), (3 * f0, "3f0")):
        ax.axvline(x / 1e9, color="gray", ls="--", lw=0.8)
        ax.text(x / 1e9, 1, t, ha="center", fontsize=8)
    ax.set_xlabel("f [GHz]")
    ax.set_ylabel("[dB]")
    ax.set_ylim(-60, 4)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)
    _salvar(fig, "ex_harmonico")
    return dict(f0=f0, razao=razao, n=n, k=_k_da_razao(razao, p.periodo), l=p.periodo * 1e3,
                comp=p.comprimento * 1e3, eem=eem, fmax=b.f_max / 1e9, s21=20 * math.log10(b.s21_min),
                bw=b.banda / 1e9, s21_f0=float(abcd.db(s21[i0])), s11_f0=float(abcd.db(s11[i0])),
                s21_2f0=float(abcd.db(s21[i2])), zmin=p.z0.min(), zmax=p.z0.max(), wmin=w.min(), wmax=w.max())


def exemplo_miniaturizacao():
    """Linha de 90° (lambda/4) de 50 Ohm a 1 GHz em RO4003C: convencional x slow-wave."""
    f0, zb, n, fs = 1e9, 50.0, 4, 12e9
    fase = 90 / n
    w50 = sintese(zb, RO4003)
    _, ee50 = analise(w50, RO4003)
    comp_conv = C0 / (f0 * math.sqrt(ee50)) / 4 * 1e3
    p = cc.projeto_por_fase(fase, f0, zb, fs)
    for zp in (20, 15, 12, 10):
        try:
            geo = cc.geometria(p, RO4003, zp, zb)
            if geo.l_patch < 0.8 * geo.periodo:
                break
        except ValueError:
            continue
    f = np.array([f0])
    m = abcd.potencia(cc.celula_distribuida(f, p), n)
    s = abcd.abcd_para_s(m, zb)
    _, bl = abcd.dispersao(cc.celula_distribuida(f, p))
    comp_sw = n * geo.periodo
    return dict(f0=f0, zb=zb, n=n, fs=fs, fase=fase, w50=w50, ee50=float(ee50), comp_conv=comp_conv,
                L=p.L * 1e9, C=p.C * 1e12, Cls=p.Cls * 1e12, zo=p.zo, fc=p.fc / 1e9, sw=p.fator_slow_wave,
                zp=zp, periodo=geo.periodo, w_linha=geo.w_linha, w_patch=geo.w_patch, l_patch=geo.l_patch,
                comp_sw=comp_sw, reducao=100 * (1 - comp_sw / comp_conv),
                fase_total=float(np.rad2deg(bl[0]) * n), s11=float(abcd.db(s[0, 0, 0])),
                fase_s21=float(np.angle(s[0, 1, 0], deg=True)))


def exemplo_livro_ebg():
    p = cm.perfil("Senoidal", 64.0, 0.01, 10)
    w, ee = cm.larguras(p, RO3010)
    eem = cm.eeff_medio(p.z, ee)
    b = cm.bandas_analiticas(p, eem)[0]
    return dict(eem=eem, fmax=b.f_max / 1e9, s21=20 * math.log10(b.s21_min), bw=b.banda / 1e9,
                s11=b.s11_max, wmin=w.min(), wmax=w.max(), w50=sintese(50, RO3010))


def exemplo_livro_carga():
    p = cc.projeto_por_fase(35, 1.4e9, 47.5, 22e9)
    sub = Substrato(2.4, 0.675, 35, 0.0015, "Arlon")
    geo = cc.geometria(p, sub, 10, 47.5)
    return dict(L=p.L * 1e9, C=p.C * 1e12, Cls=p.Cls * 1e12, zo=p.zo, fc=p.fc / 1e9, sw=p.fator_slow_wave,
                periodo=geo.periodo, w_linha=geo.w_linha, w_patch=geo.w_patch, l_patch=geo.l_patch)


def exemplo_personalizados():
    """Bandas obtidas com cada modelo de perfil personalizado (RO3010, 50 ohm)."""
    from lta.core import perfil_personalizado as pp
    res = []
    for nome, (tipo, expr, comp, passo, _) in pp.MODELOS.items():
        p, _ = pp.perfil(tipo, expr, comp, passo, RO3010, 50.0)
        w, ee = cm.larguras(p, RO3010)
        f = np.linspace(0.5e9, 18e9, 1800)
        _, s21 = cm.resolver_cmt(f, p.z, p.z0, ee)
        bandas = pp.bandas_detectadas(f, s21)
        res.append(dict(nome=nome, tipo=tipo, expr=expr, comp=comp, zmin=p.z0.min(), zmax=p.z0.max(),
                        wmin=w.min(), wmax=w.max(),
                        bandas=[(c / 1e9, d) for _, _, c, d in bandas]))
    return res


def exemplo_janelamento():
    """Faixa útil abaixo da banda rejeitada: perfil senoidal com e sem janela (RO3010, l = 10 mm)."""
    from lta.core import perfil_personalizado as pp
    f = np.linspace(0.1e9, 12e9, 1500)
    res = []
    for nome, expr, comp in (("Sem janela (10 períodos)", "1.28*cos(2*pi*z/10)", 100),
                             ("Hamming (10 períodos)", "2.0*janela(z)*cos(2*pi*z/10)", 100),
                             ("Hann (14 períodos)", "2.0*hann(z)*cos(2*pi*z/10)", 140)):
        p, _ = pp.perfil(pp.TIPOS_FUNCAO[0], expr, comp, 0.1, RO3010)
        w, ee = cm.larguras(p, RO3010)
        s11, s21 = cm.resolver_cmt(f, p.z, p.z0, ee)
        s = np.empty((f.size, 2, 2), complex)
        s[:, 0, 0], s[:, 1, 0], s[:, 0, 1], s[:, 1, 1] = s11, s21, s21, s11
        baixa = [x for x in sorted(abcd.banda_util(f, s)) if x[0] < 1e9]
        res.append(dict(nome=nome, expr=expr, zmin=p.z0.min(), zmax=p.z0.max(), s21=float(abcd.db(s21).min()),
                        f_util=baixa[0][1] / 1e9 if baixa else 0.0))
    return res


def exemplo_passa_tudo():
    """Banda útil da escada LC x treliça x T em ponte (L = 2,5 nH, C = 1 pF; ponte C = 1 pF, L = 1,25 nH)."""
    f = np.linspace(1e7, 20e9, 2001)
    casos = (("Escada T (Fig. 2.8)", "T", dict(ls=2.5e-9, cp=1e-12)),
             ("Treliça (X)", "Treliça (X)", dict(ls=2.5e-9, cp=1e-12)),
             ("T em ponte", "T em ponte", dict(cs=1e-12, lp=1.25e-9)))
    res = []
    for nome, top, el in casos:
        m = abcd.celula(f, top, **el)
        zb, _ = abcd.impedancia_bloch(m)
        al, _ = abcd.dispersao(m)
        s = abcd.abcd_para_s(abcd.potencia(m, 20), 50)
        a, b = abcd.banda_util(f, s)[0]
        res.append(dict(nome=nome, zb_min=zb[al < 1e-6].real.min(), zb_max=zb[al < 1e-6].real.max(),
                        f_util=b / 1e9, ate_fim=b >= f[-1]))
    return res


def gerar_tudo():
    return dict(conceito=fig_dispersao_conceito(), cls=fig_impacto_cls(), ebg=fig_impacto_ebg(),
                periodo=fig_projeto_periodo(), refletor=exemplo_refletor(), harmonico=exemplo_harmonico(),
                mini=exemplo_miniaturizacao(), livro_ebg=exemplo_livro_ebg(), livro_carga=exemplo_livro_carga(),
                pers=exemplo_personalizados(), janela=exemplo_janelamento(), passa_tudo=exemplo_passa_tudo())


if __name__ == "__main__":
    import pprint
    pprint.pprint(gerar_tudo(), width=110)
