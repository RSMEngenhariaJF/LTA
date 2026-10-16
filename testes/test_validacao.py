"""Validação contra valores publicados em F. Martín (2015), cap. 2.

Executar a partir da raiz do projeto:  python -m unittest discover testes
"""
import struct
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lta.core import abcd, carga_capacitiva as cc, modos_acoplados as cm  # noqa: E402
from lta.core import microstrip as ms  # noqa: E402
from lta.core.microstrip import Substrato, analise, sintese  # noqa: E402
from lta.exportacao import FORMATOS, exportar, simuladores  # noqa: E402
from lta.geometria import construtores  # noqa: E402

RO3010 = Substrato(10.2, 1.27, 35, 0.0022, "RO3010")
ARLON = Substrato(2.4, 0.675, 35, 0.0015, "Arlon")


class Microstrip(unittest.TestCase):
    def test_50ohm_ro3010(self):
        # livro: W = 1,2 mm para 50 ohm em RO3010, h = 1,27 mm (Fig. 2.16)
        self.assertAlmostEqual(sintese(50, RO3010), 1.2, delta=0.05)

    def test_inversa(self):
        w = sintese(73.3, RO3010)
        self.assertAlmostEqual(float(analise(w, RO3010)[0]), 73.3, places=4)


class ABCD(unittest.TestCase):
    f = np.linspace(1e8, 8e9, 200)

    def test_dispersao_independe_dos_planos(self):
        # eq. 2.33: T e L têm a mesma dispersão
        t = abcd.celula(self.f, "T", ls=2.5e-9, cp=1e-12)
        l_ = abcd.celula(self.f, "L", ls=2.5e-9, cp=1e-12)
        np.testing.assert_allclose(abcd.dispersao(t)[1], abcd.dispersao(l_)[1], atol=1e-9)

    def test_bloch_assimetrica(self):
        # eq. 2.32: Re iguais, Im = Zs/2 na célula L
        t = abcd.impedancia_bloch(abcd.celula(self.f, "T", ls=2.5e-9, cp=1e-12))[0]
        l_ = abcd.impedancia_bloch(abcd.celula(self.f, "L", ls=2.5e-9, cp=1e-12))[0]
        passa = self.f < 1 / (np.pi * np.sqrt(2.5e-9 * 1e-12))
        np.testing.assert_allclose(t[passa].real, l_[passa].real, rtol=1e-9)
        np.testing.assert_allclose(l_[passa].imag, (2 * np.pi * self.f * 2.5e-9 / 2)[passa], rtol=1e-9)
        self.assertAlmostEqual(t[0].real, 50.0, delta=0.1)

    def test_iterativa_igual_bloch(self):
        m = abcd.celula(self.f, "L", ls=2.5e-9, cp=1e-12)
        np.testing.assert_allclose(abcd.impedancia_bloch(m)[0], abcd.impedancia_iterativa(m), rtol=1e-9)


class ModosAcoplados(unittest.TestCase):
    def _caso(self, tipo, k):
        p = cm.perfil(tipo, k, 0.01, 10)
        w, ee = cm.larguras(p, RO3010)
        return p, ee, cm.eeff_medio(p.z, ee)

    def test_senoidal_fig_2_11(self):
        p, ee, em = self._caso("Senoidal", 64.0)
        self.assertAlmostEqual(p.z0.max(), 75.0, delta=0.5)
        self.assertAlmostEqual(p.z0.min(), 33.3, delta=0.3)
        b = cm.bandas_analiticas(p, em)[0]
        self.assertAlmostEqual(20 * np.log10(b.s21_min), -49.6, delta=0.1)
        self.assertAlmostEqual(b.f_max / 1e9, 5.66, delta=0.15)   # eeff difere do EM do livro
        self.assertAlmostEqual(b.banda / 1e9, 2.55, delta=0.1)
        f = np.linspace(4e9, 8e9, 400)
        s11, s21 = cm.resolver_cmt(f, p.z, p.z0, ee)
        a11, a21 = cm.resolver_abcd(f, p.z, p.z0, ee)
        self.assertAlmostEqual(abcd.db(s21).min(), -49.6, delta=1.0)
        np.testing.assert_allclose(abcd.db(s21), abcd.db(a21), atol=0.5)
        # sem perdas: |S11|^2 + |S21|^2 = 1
        np.testing.assert_allclose(abs(s11) ** 2 + abs(s21) ** 2, 1, atol=1e-6)

    def test_quadrada_fig_2_14(self):
        p, _, em = self._caso("Quadrada", 100.0)
        esperado = {1: -49.3, 3: -12.5, 5: -5.7}
        for b in cm.bandas_analiticas(p, em):
            if b.n in esperado:
                self.assertAlmostEqual(20 * np.log10(b.s21_min), esperado[b.n], delta=0.15)


class PassaTudo(unittest.TestCase):
    f = np.linspace(1e7, 30e9, 600)

    def test_trelica_resistencia_constante(self):
        # Martín eqs. 6.15-6.16: L em série e C nas diagonais -> ZB = sqrt(L/C) em todas as frequências
        m = abcd.celula(self.f, "Treliça (X)", ls=2.5e-9, cp=1e-12)
        zb, _ = abcd.impedancia_bloch(m)
        np.testing.assert_allclose(zb, 50.0, rtol=1e-9)
        s = abcd.abcd_para_s(abcd.potencia(m, 20), 50)
        self.assertLess(abcd.db(s[:, 0, 0]).max(), -200)
        np.testing.assert_allclose(np.abs(s[:, 1, 0]), 1.0, atol=1e-9)

    def test_t_em_ponte_contra_analise_nodal(self):
        w = 2 * np.pi * self.f
        z1, z2 = 1 / (1j * w * 1e-12), 1j * w * 1.25e-9
        m = abcd.t_em_ponte(z1, z2)
        for i in (0, 250, 599):
            y1, y2, y3 = 1 / z1[i], 1 / z2[i], 1 / (4 * z2[i])
            ym = np.array([[y1 + y3, -y3, -y1], [-y3, y1 + y3, -y1], [-y1, -y1, 2 * y1 + y2]])
            yr = ym[:2, :2] - np.outer(ym[:2, 2], ym[2, :2]) / ym[2, 2]
            ref = np.array([[-yr[1, 1] / yr[1, 0], -1 / yr[1, 0]],
                            [-(yr[0, 0] * yr[1, 1] - yr[0, 1] * yr[1, 0]) / yr[1, 0], -yr[0, 0] / yr[1, 0]]])
            np.testing.assert_allclose(m[i], ref, rtol=1e-9)
        zb, _ = abcd.impedancia_bloch(m)
        np.testing.assert_allclose(zb, np.sqrt(2 * 1.25e-9 / 1e-12), rtol=1e-9)   # eq. 6.29

    def test_regra_banda_util_escada(self):
        # escada LC: |S11| <= -15 dB garantido até ~0,55 fc, para qualquer número de células
        m = abcd.celula(self.f, "T", ls=2.5e-9, cp=1e-12)
        fc = 1 / (np.pi * np.sqrt(2.5e-9 * 1e-12))
        for n in (5, 20, 50):
            a, b = abcd.banda_util(self.f, abcd.abcd_para_s(abcd.potencia(m, n), 50))[0]
            self.assertGreaterEqual(b, 0.54 * fc)


class PerfilPersonalizado(unittest.TestCase):
    def test_duas_bandas(self):
        from lta.core import perfil_personalizado as pp
        tipo, expr, comp, passo, _ = pp.MODELOS["Duas bandas: 3 e 5 GHz"]
        p, avisos = pp.perfil(tipo, expr, comp, passo, RO3010)
        self.assertEqual(avisos, [])
        self.assertAlmostEqual(p.z0[0], 50, places=6)
        self.assertAlmostEqual(p.z0[-1], 50, delta=0.5)   # média de K removida
        _, ee = cm.larguras(p, RO3010)
        f = np.linspace(1e9, 7e9, 1200)
        _, s21 = cm.resolver_cmt(f, p.z, p.z0, ee)
        centros = [b[2] / 1e9 for b in pp.bandas_detectadas(f, s21)]
        self.assertEqual(len(centros), 2)
        self.assertAlmostEqual(centros[0], 3.0, delta=0.1)
        self.assertAlmostEqual(centros[1], 5.0, delta=0.1)

    def test_equivale_ao_senoidal(self):
        # o perfil do livro escrito como Z0(z) deve reproduzir o perfil senoidal interno
        from lta.core import perfil_personalizado as pp
        tipo, expr, comp, passo, _ = pp.MODELOS["Livro Fig. 2.11 escrita como Z0(z)"]
        p, _ = pp.perfil(tipo, expr, comp, passo, RO3010)
        ref = cm.perfil("Senoidal", 64.0, 0.01, 10)
        self.assertAlmostEqual(p.z0.max(), ref.z0.max(), delta=0.1)
        self.assertAlmostEqual(p.z0.min(), ref.z0.min(), delta=0.1)

    def test_expressao_segura(self):
        from lta.core import perfil_personalizado as pp
        for ruim in ('__import__("os").system("x")', "z.__class__", "open(1)", "(lambda: 1)()", '"a"'):
            with self.assertRaises(ValueError):
                pp.perfil(pp.TIPOS_FUNCAO[0], ruim, 100, 0.5, RO3010)


class CargaCapacitiva(unittest.TestCase):
    def test_inversao(self):
        p = cc.projeto_por_fc(5.6e9, 47.5, 22e9)
        self.assertAlmostEqual(p.fc, 5.6e9, delta=1)
        self.assertAlmostEqual(p.zb, 47.5, places=6)
        self.assertAlmostEqual(p.fs, 22e9, delta=1)

    def test_exemplo_livro(self):
        # livro (Fig. 2.34): L = 3,1 nH, C = 0,17 pF, Cls = 1,2 pF (valores arredondados/otimizados)
        p = cc.projeto_por_fase(35, 1.4e9, 47.5, 22e9)
        self.assertAlmostEqual(p.L * 1e9, 3.1, delta=0.35)
        self.assertAlmostEqual(p.C * 1e12, 0.17, delta=0.03)
        self.assertAlmostEqual(p.Cls * 1e12, 1.2, delta=0.15)


class PerdasDispersao(unittest.TestCase):
    def test_profundidade_pelicular(self):
        # livro, nota 24 (cap. 1): cobre, 2,06 µm a 1 GHz e 0,66 µm a 10 GHz
        d = ms.profundidade_pelicular(np.array([1e9, 10e9])) * 1e6
        self.assertAlmostEqual(d[0], 2.06, delta=0.05)
        self.assertAlmostEqual(d[1], 0.66, delta=0.02)

    def test_alfa_condutor_tira_larga(self):
        # eq. 1.90 tende a eq. 1.89 (αc = Rs/(Z0 W)) para W/h >> 1
        w = 40 * RO3010.h_mm
        z0, _ = analise(w, RO3010)
        ref = ms.resistencia_superficial(1e10, RO3010) / (z0 * w * 1e-3)
        self.assertAlmostEqual(float(ms.alfa_condutor(w, RO3010, 1e10)[0, 0]) / ref, 1.0, delta=0.15)

    def test_alfa_condutor_continuo(self):
        # emendas entre as expressões de Pucel e a eq. 1.90
        for u in (1 / (2 * np.pi), 2.0):
            a = ms.alfa_condutor(np.array([0.999, 1.001]) * u * RO3010.h_mm, RO3010, 1e10)[0]
            self.assertLess(abs(a[1] / a[0] - 1), 0.05)

    def test_alfa_dieletrico_tem(self):
        # com εeff = εr, a eq. 1.68 se reduz à eq. 1.67: αd = (β/2)·tanδ
        f = 5e9
        ad = float(ms.alfa_dieletrico(1.0, RO3010, f, eeff=RO3010.er)[0, 0])
        beta = 2 * np.pi * f * np.sqrt(RO3010.er) / ms.C0
        self.assertAlmostEqual(ad, beta / 2 * RO3010.tand, places=9)

    def test_dispersao_limites(self):
        w = sintese(50, RO3010)
        ee0 = float(analise(w, RO3010)[1])
        ee = ms.eeff_dispersivo(w, RO3010, np.array([1e6, 5e9, 20e9, 100e9]))[:, 0]
        self.assertAlmostEqual(ee[0], ee0, places=4)
        self.assertTrue(np.all(np.diff(ee) > 0))
        self.assertLess(ee[-1], RO3010.er)

    def test_ebg_dispersivo_fig_2_11(self):
        # com dispersão, o centro da 1ª banda se aproxima do valor do livro (5,66 GHz, solver EM)
        p = cm.perfil("Senoidal", 64.0, 0.01, 10)
        f = np.linspace(5e9, 6.5e9, 301)
        centros = []
        for disp in (False, True):
            z0, ee, al = cm.parametros_linha(f, p, RO3010, dispersivo=disp)
            _, s21 = cm.resolver_cmt(f, p.z, z0, ee, al)
            centros.append(f[np.argmin(abs(s21))] / 1e9)
        self.assertAlmostEqual(centros[1], 5.66, delta=0.08)
        self.assertLess(abs(centros[1] - 5.66), abs(centros[0] - 5.66))

    def test_cmt_com_perdas(self):
        p = cm.perfil("Senoidal", 64.0, 0.01, 10)
        f = np.linspace(1e9, 8e9, 300)
        z0, ee, al = cm.parametros_linha(f, p, RO3010, perdas=True, dispersivo=True)
        s11, s21 = cm.resolver_cmt(f, p.z, z0, ee, al)
        a11, a21 = cm.resolver_abcd(f, p.z, z0, ee, 50, al)
        ok = abcd.db(s21) > -30
        np.testing.assert_allclose(abcd.db(s21)[ok], abcd.db(a21)[ok], atol=0.5)
        self.assertTrue(np.all(abs(s11) ** 2 + abs(s21) ** 2 < 1))
        # sem as opções, igual ao modelo original
        z0, ee, al = cm.parametros_linha(f, p, RO3010)
        _, ee1 = cm.larguras(p, RO3010)
        np.testing.assert_allclose(cm.resolver_cmt(f, p.z, z0, ee, al)[1], cm.resolver_cmt(f, p.z, p.z0, ee1)[1])

    def test_carga_microstrip_sem_opcoes(self):
        p = cc.projeto_por_fase(35, 1.4e9, 47.5, 22e9)
        geo = cc.geometria(p, ARLON, 10, 47.5)
        f = np.linspace(1e8, 12e9, 50)
        np.testing.assert_allclose(cc.celula_microstrip(f, p, geo, ARLON), cc.celula_distribuida(f, p))
        al, _ = abcd.dispersao(cc.celula_microstrip(f, p, geo, ARLON, perdas=True))
        self.assertTrue(np.all(al[f < 0.5 * p.fc] > 0))

    def test_atraso_grupo_linha(self):
        # linha uniforme casada: τg = l·√εeff / c
        w = sintese(50, RO3010)
        z0, ee = (float(x) for x in analise(w, RO3010))
        f = np.linspace(1e8, 10e9, 2000)
        m = abcd.linha(z0, 2 * np.pi * f * np.sqrt(ee) / ms.C0 * 0.05)
        tau = abcd.atraso_grupo(f, abcd.abcd_para_s(m, z0)[:, 1, 0])
        np.testing.assert_allclose(tau, 0.05 * np.sqrt(ee) / ms.C0, rtol=1e-6)


class Floquet(unittest.TestCase):
    def test_linha_uniforme_fig_2_1(self):
        # perturbação infinitesimal: β estendido é a reta k0·√εeff e os harmônicos, retas deslocadas
        from lta.core import floquet as fl
        w = sintese(50, RO3010)
        z0, ee = (float(x) for x in analise(w, RO3010))
        l = 0.01
        f = np.linspace(1e8, 25e9, 1500)
        al, bl = abcd.dispersao(abcd.linha(z0, 2 * np.pi * f * np.sqrt(ee) / ms.C0 * l))
        ext = fl.desdobrar(bl)
        # atol: arccosh perde precisão perto de βl = mπ (cos βl = ±1)
        np.testing.assert_allclose(ext, 2 * np.pi * f * np.sqrt(ee) / ms.C0 * l, atol=1e-2)
        h = {x.n: x for x in fl.harmonicos(f, ext, l, 1)}
        f_ini, f_fim, t_ini, _ = fl.faixas_radiacao(f, h[-1])[0]
        self.assertAlmostEqual(f_ini / fl.inicio_radiacao(l, ee), 1, delta=0.01)
        self.assertAlmostEqual(t_ini, 180, delta=10)            # começa em backfire
        i = np.argmin(abs(f - fl.frequencia_broadside(l, ee)))
        self.assertAlmostEqual(h[-1].theta_graus[i], 90, delta=1)   # eq. 2.11 com βn = 0
        self.assertFalse(h[0].rapido.any())                     # fundamental é onda lenta

    def test_ebg_desdobrado_com_perdas(self):
        from lta.core import floquet as fl
        p = cm.perfil("Senoidal", 64.0, 0.01, 1)
        f = np.linspace(1e8, 25e9, 1500)
        z0, ee, al = cm.parametros_linha(f, p, RO3010, perdas=True, dispersivo=True)
        _, bl = abcd.dispersao(cm.abcd_perfil(f, p.z, z0, ee, al))
        ext = fl.desdobrar(bl)
        self.assertTrue(np.all(np.diff(ext) >= -1e-3))
        i = np.argmin(abs(f - 5.66e9))                          # 1ª banda: βl = π
        self.assertAlmostEqual(ext[i] / np.pi, 1.0, delta=0.02)


class AplicacoesEBG(unittest.TestCase):
    def test_cdl_exemplo_livro(self):
        # livro, p. 101: f0 = 9 GHz, Δf = 12 GHz, ψ = -0,5 ns/GHz -> l0 = 6,4 mm, L = 50·l0, C = -2080 1/m²
        from lta.core import aplicacoes_ebg as ap
        p = ap.projeto_cdl(9e9, 12e9, -0.5e-18, RO3010, n_periodos=50)
        self.assertAlmostEqual(p.l0 * 1e3, 6.4, delta=0.05)
        self.assertAlmostEqual(p.chirp / -2080, 1, delta=0.03)
        self.assertAlmostEqual(p.f_local(0) / 1e9, 9.0, places=6)
        self.assertAlmostEqual((p.f_local(-p.comprimento / 2) - p.f_local(p.comprimento / 2)) / 1e9, 12, places=6)
        perf = ap.perfil_cdl(p)
        self.assertAlmostEqual(perf.z0[0], 50, places=6)
        self.assertAlmostEqual(perf.z0[-1], 50, places=6)
        f = np.linspace(2e9, 16e9, 1200)
        _, ee = cm.larguras(perf, RO3010)
        s11, _ = cm.resolver_cmt(f, perf.z, perf.z0, ee)
        psi, _ = ap.inclinacao_atraso(f, abcd.atraso_grupo(f, s11), 4.2e9, 13.8e9)
        self.assertAlmostEqual(psi / p.psi, 1, delta=0.1)

    def test_cavidade_exemplo_livro(self):
        # livro, p. 93: φref1 = φref2 = -95,84° em 4,2 GHz -> βLm = 84,15° com n = -1
        from lta.core import aplicacoes_ebg as ap
        lm, n = ap.comprimento_cavidade(np.radians(-95.84), np.radians(-95.84), 1.0)
        self.assertEqual(n, -1)
        self.assertAlmostEqual(np.degrees(lm), 84.15, delta=0.02)

    def test_ressonador_em_f0(self):
        from lta.core import aplicacoes_ebg as ap
        r = ap.projeto_ressonador(4.2e9, 64.0, 3, RO3010)
        f = np.linspace(4.1e9, 4.3e9, 2001)
        _, ee = cm.larguras(r.perfil, RO3010)
        _, s21 = cm.resolver_cmt(f, r.perfil.z, r.perfil.z0, ee)
        fp, il, ql, _ = ap.fator_q(f, s21)
        self.assertAlmostEqual(fp / 4.2e9, 1, delta=1e-3)
        self.assertLess(il, 0.01)
        self.assertGreater(ql, 100)


class CRLH(unittest.TestCase):
    f = np.linspace(5e7, 8e9, 4000)

    def test_balanceada_eq_3_66(self):
        from lta.core import crlh
        c = crlh.balanceada_por_banda(1e9, 6e9, 50)
        self.assertTrue(c.balanceada)
        b = c.bandas()
        self.assertAlmostEqual(b["f_cl"] / 1e9, 1.0, places=6)
        self.assertAlmostEqual(b["f_cr"] / 1e9, 6.0, places=6)
        f0 = crlh.frequencia_transicao(c)
        w0, wr = 2 * np.pi * f0, c.wR                     # eq. 3.66 na forma do livro
        self.assertAlmostEqual(wr * (np.sqrt(1 + w0**2 / wr**2) - 1) / (2 * np.pi * 1e9), 1.0, places=6)
        self.assertAlmostEqual(abs(c.zb(np.array([f0 * 0.999]))[0]), 50, delta=0.1)   # eq. 3.63

    def test_dispersao_analitica_e_sinal(self):
        from lta.core import crlh
        for tipo in crlh.TIPOS:
            c = crlh.CelulaCRLH(2.5e-9, 1e-12, 3e-9, 2e-12, tipo)
            m = c.abcd(self.f)
            w = 2 * np.pi * self.f
            ok = (abs(w / c.ws - 1) > 0.01) & (abs(w / c.wp - 1) > 0.01)
            np.testing.assert_allclose(c.cos_bl(self.f)[ok], (0.5 * (m[:, 0, 0] + m[:, 1, 1])).real[ok],
                                       atol=1e-8)                       # eqs. 3.51 e 3.64
        c = crlh.CelulaCRLH(2.5e-9, 1e-12, 3e-9, 2e-12)
        bl, al = c.beta_l(self.f)
        g1, g2 = c.gap[0] / (2 * np.pi), c.gap[1] / (2 * np.pi)
        self.assertTrue(np.all(bl[(al < 1e-6) & (self.f < g1)] < 0))    # LH abaixo do gap
        self.assertTrue(np.all(bl[(al < 1e-6) & (self.f > g2)] > 0))    # RH acima
        self.assertTrue(np.all(al[(self.f > g1 * 1.001) & (self.f < g2 * 0.999)] > 0))  # eqs. 3.57-3.58
        eps, mu = c.meio_efetivo(self.f, 6e-3)
        lh = (al < 1e-6) & (bl < 0)
        self.assertTrue(np.all(eps[lh] < 0) and np.all(mu[lh] < 0))     # eqs. 3.61-3.62

    def test_sintese_por_fases(self):
        from lta.core import crlh
        c = crlh.balanceada_por_fases(0.9e9, np.radians(-30), 1.8e9, np.radians(30), 50)
        bl, _ = c.beta_l(np.array([0.9e9, 1.8e9]))
        np.testing.assert_allclose(np.degrees(bl), [-30, 30], atol=1.5)   # aproximação de grande λ


class Extracao(unittest.TestCase):
    # valores extraídos publicados no livro (Tabelas 3.1 e 3.3) e valores típicos para OSRR/OCSRR
    CASOS = {
        "CSRR (microstrip)": dict(L=5.78e-9, C=0.95e-12, Lc=1.65e-9, Cc=2.81e-12),
        "CSRR + gap (microstrip)": dict(L=5.78e-9, Cg=0.31e-12, C=6.43e-12, Lc=1.49e-9, Cc=3.10e-12),
        "SRR (CPW)": dict(C=2.07e-12, L=1.91e-9, Ls=0.24e-9, Cs=23.01e-12),
        "SRR + fita (CPW)": dict(C=2.07e-12, Lp=1.38e-9, L=4.91e-9, Ls=1.68e-9, Cs=3.90e-12),
        "OSRR (CPW)": dict(C=0.5e-12, Ls=6.16e-9, Cs=0.515e-12),
        "OCSRR (CPW)": dict(L=0.5e-9, Lp=1.55e-9, Cp=1.91e-12),
    }

    def test_ida_e_volta(self):
        from lta.core import extracao as ex
        f = np.linspace(0.2e9, 6e9, 3001)
        for nome, p in self.CASOS.items():
            q, _, err = ex.extrair(nome, f, ex.s_modelo(nome, p, f))
            for k in p:
                self.assertAlmostEqual(q[k] / p[k], 1, delta=0.005, msg=f"{nome}: {k}")
            self.assertLess(err, 2e-3)

    def test_touchstone_ida_e_volta(self):
        from lta.core import extracao as ex
        f = np.linspace(0.5e9, 5e9, 901)
        s = ex.s_modelo("CSRR + gap (microstrip)", self.CASOS["CSRR + gap (microstrip)"], f)
        with tempfile.TemporaryDirectory() as d:
            simuladores.touchstone(Path(d) / "c.s2p", f, s, 50)
            f2, s2, z0 = simuladores.ler_touchstone(Path(d) / "c.s2p")
            # formato MA em MHz, com dados em duas linhas
            linhas = ["! teste", "# MHz S MA R 50"]
            for fi, si in zip(f, s):
                v = [si[0, 0], si[1, 0], si[0, 1], si[1, 1]]
                linhas.append(f"{fi / 1e6:.6f} " + " ".join(f"{abs(x):.9e} {np.degrees(np.angle(x)):.9f}" for x in v[:2]))
                linhas.append(" ".join(f"{abs(x):.9e} {np.degrees(np.angle(x)):.9f}" for x in v[2:]))
            (Path(d) / "m.s2p").write_text("\n".join(linhas))
            f3, s3, _ = simuladores.ler_touchstone(Path(d) / "m.s2p")
        np.testing.assert_allclose(f2, f)
        np.testing.assert_allclose(s2, s, atol=1e-8)
        np.testing.assert_allclose(f3, f, rtol=1e-9)
        np.testing.assert_allclose(s3, s, atol=1e-7)
        self.assertEqual(z0, 50)

    def test_robustez_perdas(self):
        # dados com perdas no ressonador (R ‖ tanque), fora do modelo sem perdas do Apêndice G
        from lta.core import extracao as ex
        f = np.linspace(0.5e9, 5e9, 2001)
        w = 2 * np.pi * f
        p = self.CASOS["CSRR + gap (microstrip)"]
        tanque = 1 / (1 / (1j * w * p["Lc"]) + 1j * w * p["Cc"] + 1 / 5e3)
        braco = 1j * w * p["L"] / 2 + 1 / (2j * w * p["Cg"])
        m = abcd.cascata(abcd.series(braco), abcd.shunt(1 / (1 / (1j * w * p["C"]) + tanque)), abcd.series(braco))
        q, _, err = ex.extrair("CSRR + gap (microstrip)", f, abcd.abcd_para_s(m, 50))
        for k in p:
            self.assertAlmostEqual(q[k] / p[k], 1, delta=0.02, msg=k)
        self.assertLess(err, 0.02)

    def test_transformacao_srr_tabelas_3_1_e_3_2(self):
        # eqs. 3.71-3.74: modelo físico (Tabela 3.2) -> π equivalente (Tabela 3.1)
        from lta.core import extracao as ex
        l_, lp_, ls_, cs_ = ex.srr_cpw_transformado(1.94e-9, 1.23e-9, 17.16e-9, 0.31e-12, 0.20e-9)
        for calc, livro in ((l_, 4.91e-9), (lp_, 1.38e-9), (ls_, 1.68e-9), (cs_, 3.90e-12)):
            self.assertAlmostEqual(calc / livro, 1, delta=0.03)
        l_, lp_, ls_, cs_ = ex.srr_cpw_transformado(2.15e-9, 1.45e-9, 17.57e-9, 0.31e-12)
        for calc, livro in ((l_, 1.91e-9), (ls_, 0.24e-9), (cs_, 23.01e-12)):
            self.assertAlmostEqual(calc / livro, 1, delta=0.03)
        self.assertTrue(np.isinf(lp_))


class LinhaCSRR(unittest.TestCase):
    def test_eqs_3_86_a_3_92(self):
        from lta.core import csrr as cs
        c = cs.CircuitoCSRR(5.78e-9, 0.31e-12, 6.43e-12, 1.49e-9, 3.10e-12)     # Tabela 3.3, com gap
        f = np.linspace(0.3e9, 6e9, 3000)
        zbn, _ = abcd.impedancia_bloch(c.abcd(f))
        zba = c.zb(f)
        ok = np.abs(zbn) < 1e4
        np.testing.assert_allclose(np.abs(zba[ok]), np.abs(zbn[ok]), rtol=1e-8)        # eq. 3.90
        s21 = abcd.abcd_para_s(c.abcd(f), 50)[:, 1, 0]
        self.assertAlmostEqual(f[np.argmin(abs(s21))] / c.fz, 1, delta=0.003)          # eq. 3.86
        cs_, cpar = c.modelo_melhorado()
        self.assertAlmostEqual((2 * cs_ + cpar) / c.Cg, 1, places=9)                   # eq. 3.87
        self.assertAlmostEqual(cpar * (2 * cs_ + cpar) / cs_ / c.C, 1, places=9)       # eq. 3.88

    def test_projeto_balanceado(self):
        from lta.core import csrr as cs
        c = cs.projeto_balanceado(2e9, 50, 1.4e9, 3.5e-9)
        self.assertTrue(c.balanceada)                                                  # eq. 3.91
        self.assertAlmostEqual(c.fz / 1e9, 1.4, places=6)
        self.assertAlmostEqual(abs(c.zb(np.array([2e9 * 1.0001]))[0]), 50, delta=0.05)   # eq. 3.92
        from lta.core import extracao as ex                                            # projeto -> extração
        f = np.linspace(0.5e9, 4e9, 1501)
        q, _, err = ex.extrair("CSRR + gap (microstrip)", f, ex.s_modelo("CSRR + gap (microstrip)", c.params(), f))
        for k, v in c.params().items():
            self.assertAlmostEqual(q[k] / v, 1, delta=0.01, msg=k)

    def test_estimativa_geometrica_fig_3_39(self):
        # CSRR da Fig. 3.39 (c = 0,365, d = 0,24, rext = 3,661 mm): fo extraído ≈ 2,34 GHz (Tabela 3.3)
        from lta.core import csrr as cs
        _, _, fo = cs.estimativa_csrr(cs.GeometriaCSRR(3.661, 0.365, 0.24), RO3010)
        self.assertAlmostEqual(fo / 2.34e9, 1, delta=0.20)
        r = cs.raio_para_fo(2.34e9, 0.365, 0.24, RO3010, calibracao=2.34e9 / fo)
        self.assertAlmostEqual(r, 3.661, delta=0.01)

    def test_interdigital_bahl(self):
        from lta.core import csrr as cs
        c1 = cs.capacitancia_interdigital(6, 1.0, 0.1, RO3010)
        self.assertAlmostEqual(cs.capacitancia_interdigital(6, 2.0, 0.1, RO3010) / c1, 2, places=9)  # ∝ ℓ
        self.assertGreater(cs.capacitancia_interdigital(8, 1.0, 0.1, RO3010), c1)
        self.assertTrue(0.05e-12 < c1 < 2e-12)


class CLLoaded(unittest.TestCase):
    SUB = Substrato(2.2, 0.787, 35, 0.0009, "RT5880")

    def test_projeto_reproduz_elementos(self):
        from lta.core import cl_loaded as cl, crlh
        from lta.core.csrr import capacitancia_interdigital
        c = crlh.CelulaCRLH(2.45e-9, 0.5e-12, 3.38e-9, 0.68e-12)     # f0 ≈ 3,9 GHz (Fig. 3.28)
        self.assertAlmostEqual(crlh.frequencia_transicao(c) / 1e9, 3.9, delta=0.05)
        g, avisos = cl.projeto(c, self.SUB, n_dedos=14)
        self.assertEqual(avisos, [])
        _, ee = analise(g.w_linha, self.SUB)
        self.assertAlmostEqual(np.sqrt(c.LR / c.CR) * np.sqrt(float(ee)) * g.periodo * 1e-3 / ms.C0 / c.LR, 1, places=6)
        self.assertAlmostEqual(capacitancia_interdigital(14, g.comp_dedo, 0.1, self.SUB) / c.CL, 1, places=6)
        self.assertAlmostEqual((g.l_stub + g.l_via) / c.LL, 1, places=9)
        self.assertTrue(0.05e-9 < g.l_via < 0.5e-9)
        with self.assertRaises(ValueError):                           # LL menor que a via
            cl.projeto(crlh.CelulaCRLH(2.45e-9, 0.5e-12, 0.05e-9, 0.68e-12), self.SUB)

    def test_extracao_crlh(self):
        from lta.core import crlh, extracao as ex
        f = np.linspace(0.3e9, 8e9, 3001)
        for c in (crlh.CelulaCRLH(2.45e-9, 0.5e-12, 3.38e-9, 0.68e-12), crlh.balanceada_por_banda(1e9, 6e9, 50),
                  crlh.CelulaCRLH(2.5e-9, 1e-12, 3e-9, 2e-12)):
            p = dict(LR=c.LR, CR=c.CR, LL=c.LL, CL=c.CL)
            q, _, err = ex.extrair("CRLH (CL-loaded)", f, ex.s_modelo("CRLH (CL-loaded)", p, f))
            for k in p:
                self.assertAlmostEqual(q[k] / p[k], 1, delta=1e-3, msg=k)


class Multibanda(unittest.TestCase):
    def _confere(self, m_func, alvos, z):
        for f, fase, zz in alvos:
            m = m_func(np.array([f]))
            s21 = abcd.abcd_para_s(m, z)[0, 1, 0]
            self.assertAlmostEqual(np.degrees(-np.angle(s21)), np.degrees(fase), delta=0.01)
            self.assertAlmostEqual(abcd.impedancia_bloch(m)[0][0].real, zz, delta=0.01)

    def test_csrr_dual_livro(self):
        from lta.core import multibanda as mb
        # livro, p. 258: f1 = 0,9, f2 = 1,8, fz = 0,5 GHz, Zo = 35,35 Ω
        c = mb.csrr_dual(0.9e9, -np.pi / 2, 35.35, 1.8e9, np.pi / 2, 35.35, 0.5e9)[0]
        for calc, livro in ((c.L, 12.5e-9), (c.C, 24.9e-12), (c.Cg, 1.25e-12), (c.Lc, 3.38e-9), (c.Cc, 5.10e-12)):
            self.assertAlmostEqual(calc / livro, 1, delta=0.01)
        # livro, p. 259: trechos de linha com φL(f1) = 15° -> φ = -120° e +30°
        q1, q2 = mb.compensar_linhas(-np.pi / 2, np.pi / 2, 0.9e9, 1.8e9, 15)
        self.assertAlmostEqual(np.degrees(q1), -120, places=9)
        self.assertAlmostEqual(np.degrees(q2), 30, places=9)
        c = mb.csrr_dual(0.9e9, q1, 35.35, 1.8e9, q2, 35.35, 0.5e9)[0]
        for calc, livro in ((c.L, 9.45e-9), (c.C, 17.9e-12), (c.Cg, 1.01e-12), (c.Lc, 4.85e-9), (c.Cc, 2.95e-12)):
            self.assertAlmostEqual(calc / livro, 1, delta=0.01)
        self._confere(c.abcd, [(0.9e9, q1, 35.35), (1.8e9, q2, 35.35)], 35.35)
        self.assertAlmostEqual(c.fz / 1e9, 0.5, places=6)

    def test_crlh_dual_eqs_4_29_a_4_32(self):
        from lta.core import multibanda as mb
        c = mb.crlh_dual(2.4e9, -np.pi / 2, 35.35, 3.75e9, np.pi / 2, 35.35)
        w1, w2 = 2 * np.pi * 2.4e9, 2 * np.pi * 3.75e9
        self.assertTrue(c.balanceada)
        self.assertAlmostEqual(c.ws / np.sqrt(w1 * w2), 1, places=9)                    # eq. 4.31
        self.assertAlmostEqual(c.LR, 2 * 35.35 / (w2 - w1), delta=1e-15)                 # eq. 4.32a
        self.assertAlmostEqual(c.CL, (w2 - w1) / (2 * w1 * w2 * 35.35), delta=1e-18)     # eq. 4.32b
        self.assertAlmostEqual(c.LL / 0.83e-9, 1, delta=0.03)          # livro (p. 255): Lp' = 0,83 nH
        self._confere(c.abcd, [(2.4e9, -np.pi / 2, 35.35), (3.75e9, np.pi / 2, 35.35)], 35.35)
        c = mb.crlh_dual(1e9, -np.pi / 2, 30, 2e9, np.pi / 2, 60)                      # Z1 ≠ Z2: eq. 4.30
        self.assertAlmostEqual(c.ws * c.wp / (4 * np.pi**2 * 2e18), 1, places=9)

    def test_rh_plh(self):
        from lta.core import multibanda as mb
        r = mb.rh_plh(1e9, -np.pi / 2, 3.1e9, np.pi / 2, 50, 4)
        np.testing.assert_allclose(r.fase(np.array([1e9, 3.1e9])), [-np.pi / 2, np.pi / 2])   # eq. 4.20
        with self.assertRaises(ValueError):                                             # eq. 4.24
            mb.rh_plh(1e9, np.pi / 2, 3.5e9, 3 * np.pi / 2, 50)

    def test_ecrlh_quad_livro(self):
        from lta.core import multibanda as mb
        f = np.array([0.9, 1.17, 1.57, 1.8]) * 1e9
        e = mb.ecrlh_quad(f, 35.35)[0]
        livro = dict(Lhs=22.458e-9, Chs=0.766e-12, Lhp=0.883e-9, Chp=14.06e-12, Lvs=8.785e-9, Cvs=1.414e-12,
                     Lvp=1.915e-9, Cvp=8.986e-12)                                        # livro, p. 267
        for k, v in e.params().items():
            self.assertAlmostEqual(v / livro[k], 1, delta=0.05, msg=k)
        self._confere(e.abcd, [(fx, (-1) ** k * np.pi / 2, 35.35) for k, fx in enumerate(f, 1)], 35.35)


class Dispositivos(unittest.TestCase):
    def test_ideais(self):
        from lta.core import dispositivos as dv
        f = np.array([1e9])
        z = 50.0
        s = dv.branch_line(f, dv.linha_convencional(z / np.sqrt(2), np.pi / 2, 1e9),
                           dv.linha_convencional(z, np.pi / 2, 1e9))[0]
        np.testing.assert_allclose(np.abs(s[:, 0]) ** 2, [0, 0.5, 0.5, 0], atol=1e-12)
        self.assertAlmostEqual(abs(dv.diferenca_fase(s[1, 0], s[2, 0])), 90, places=6)
        zr = z * np.sqrt(2)
        s = dv.rat_race(f, dv.linha_convencional(zr, np.pi / 2, 1e9), dv.linha_convencional(zr, 3 * np.pi / 2, 1e9))[0]
        self.assertAlmostEqual(abs(dv.diferenca_fase(s[1, 0], s[2, 0])), 180, places=6)   # Δ
        self.assertAlmostEqual(dv.diferenca_fase(s[1, 3], s[2, 3]), 0, places=6)         # Σ
        s = dv.wilkinson(f, dv.linha_convencional(zr, np.pi / 2, 1e9))[0]
        np.testing.assert_allclose(np.abs(s) ** 2, [[0, .5, .5], [.5, 0, 0], [.5, 0, 0]], atol=1e-12)

    def test_branch_line_dual_band(self):
        from lta.core import dispositivos as dv
        f = np.array([0.9e9, 1.8e9])
        s = dv.branch_line(f, dv.linha_crlh_dual(50 / np.sqrt(2), 0.9e9, 1.8e9)[0],
                           dv.linha_crlh_dual(50, 0.9e9, 1.8e9)[0])
        for k in range(2):
            np.testing.assert_allclose(np.abs(s[k, :, 0]) ** 2, [0, 0.5, 0.5, 0], atol=1e-9)
        # banda LH inverte o sinal do balanço de fase (livro, Fig. 4.16: -93° e +92°)
        self.assertAlmostEqual(dv.diferenca_fase(s[0, 1, 0], s[0, 2, 0]), -dv.diferenca_fase(s[1, 1, 0], s[1, 2, 0]), 6)

    def test_limites_de_banda(self):
        from lta.core import dispositivos as dv, crlh
        self.assertAlmostEqual(dv.fase_limite_uma_celula() / np.pi, 0.7, delta=0.01)       # livro, p. 239
        self.assertAlmostEqual(dv.inclinacao_dual(-np.pi / 2, 1.0, 10000), np.pi / 2, places=6)   # eq. 4.11
        self.assertGreater(dv.inclinacao_dual(-np.pi / 2, 1.0, 1), dv.inclinacao_dual(-np.pi / 2, 1.0, 2))  # eq. 4.10
        c = crlh.balanceada_por_banda(1e9, 6e9, 50)
        f = np.array([crlh.frequencia_transicao(c)]) * np.array([0.9999, 1.0001])
        bl, _ = c.beta_l(f)
        self.assertAlmostEqual((bl[1] - bl[0]) / (2 * np.pi * (f[1] - f[0])) / dv.atraso_crlh_zero(c), 1, delta=0.01)

    def test_rat_race_banda_larga(self):
        # -90° LH no lugar de 270° alarga o balanço de fase (livro, Figs. 4.1 e 4.5)
        from lta.core import dispositivos as dv
        f = np.linspace(0.8e9, 2.4e9, 801)
        zr = 50 * np.sqrt(2)
        l90 = dv.linha_convencional(zr, np.pi / 2, 1.6e9)
        conv = dv.rat_race(f, l90, dv.linha_convencional(zr, 3 * np.pi / 2, 1.6e9))
        art = dv.rat_race(f, l90, dv.linha_plh(zr, -np.pi / 2, 1.6e9, 3)[0])
        banda = []
        for s in (conv, art):
            ok = np.abs(np.abs(dv.diferenca_fase(s[:, 1, 0], s[:, 2, 0])) - 180) < 5
            a, b = dv.faixa_em_torno(f, ok, 1.6e9)
            banda.append(b - a)
        self.assertGreater(banda[1], 3 * banda[0])


class Filtros(unittest.TestCase):
    def test_prototipos_tabelados(self):
        from lta.core import filtros as ft
        np.testing.assert_allclose(ft.prototipo("Chebyshev", 3, 0.5)[1:4], [1.5963, 1.0967, 1.5963], atol=1e-4)
        np.testing.assert_allclose(ft.prototipo("Chebyshev", 4, 0.1)[1:6], [1.1088, 1.3061, 1.7703, 0.8180, 1.3554],
                                   atol=2e-4)
        np.testing.assert_allclose(ft.prototipo("Butterworth", 3)[1:4], [1, 2, 1], atol=1e-12)
        self.assertAlmostEqual(ft.g_periodico(3), 1.521, delta=0.001)                # livro, p. 283

    def test_passa_faixa_chebyshev(self):
        from lta.core import filtros as ft
        f = np.linspace(0.85e9, 1.15e9, 12001)
        for n, rip in ((3, 0.5), (4, 0.1)):
            g = ft.prototipo("Chebyshev", n, rip)
            ceq, leq, _ = ft.ressonadores(g, 1e9, 0.1)
            s = ft.resposta_inversores(f, g, [ft.tanque(c, l) for c, l in zip(ceq, leq)])
            s21 = abcd.db(s[:, 1, 0])
            banda = np.flatnonzero(s21 >= -rip - 1e-6)
            fa, fb = f[banda[0]], f[banda[-1]]
            self.assertAlmostEqual((fb - fa) / np.sqrt(fa * fb), 0.1, delta=0.003)    # FBW de ripple
            self.assertAlmostEqual(-s21[banda].min(), rip, delta=1e-3)
            self.assertAlmostEqual(abcd.db(s[banda, 0, 0]).max(), 10 * np.log10(1 - 10 ** (-rip / 10)), delta=0.05)

    def test_dual_band_eqs_4_49_a_4_50(self):
        from lta.core import filtros as ft
        f1, f2 = 0.9e9, 1.8e9
        np.testing.assert_allclose(ft.mapa_dual(np.array([f1, f2]), f1, f2), [-1, 1])
        g = ft.prototipo("Chebyshev", 3, 0.1)
        ceq, leq, _ = ft.ressonadores(g, np.sqrt(f1 * f2), 0.1)
        ys = [ft.admitancia_dual(*ft.ressonador_dual(c, l, f1, f2)) for c, l in zip(ceq, leq)]
        s = ft.resposta_inversores(np.array([f1, f2, 1.3e9]), g, ys)
        np.testing.assert_allclose(np.abs(s[:2, 1, 0]), 1, atol=1e-9)                # centro das duas bandas
        self.assertLess(abcd.db(s[2, 1, 0]), -30)                                      # rejeição entre bandas

    def test_celula_hibrida(self):
        from lta.core import filtros as ft
        p = ft.celula_hibrida(1e9, 0.1, 1.521, 2e9)
        self.assertAlmostEqual(p["Cg"] / 1.59e-12, 1, delta=0.005)                      # eq. 4.46, Tabela 4.2
        w = 2 * np.pi * 1e9
        zb = 1 / (1j * w * p["C"]) + 1j * w * p["Lc"] / (1 - w**2 * p["Lc"] * p["Cc"])
        self.assertAlmostEqual((1 / (1 / (1j * w * p["Ld"]) + 1 / zb)).imag, 50, delta=1e-6)   # Zp(ω0) = jZ0
        self.assertAlmostEqual(1 / (2 * np.pi * np.sqrt(p["Lc"] * (p["C"] + p["Cc"]))), 2e9, delta=1)  # eq. 3.86
        f = np.linspace(0.5e9, 3e9, 25001)
        s = abcd.abcd_para_s(abcd.potencia(ft.abcd_hibrida(f, p), 3), 50)
        a, b, fc = ft.banda_3db(f, s[:, 1, 0])
        self.assertAlmostEqual(fc / 1e9, 1.0, delta=0.01)
        self.assertAlmostEqual((b - a) / fc, 0.1, delta=0.01)


class LWA(unittest.TestCase):
    def test_eq_4_51(self):
        from lta.core import floquet as fl
        f, d = 2.5e9, 0.02
        n = fl.celulas_lwa(0.076, f, d, 0.95)
        alfa = 0.076 * 2 * np.pi * f / ms.C0
        self.assertLessEqual(np.exp(-2 * alfa * n * d), 0.05)
        self.assertGreater(np.exp(-2 * alfa * (n - 1) * d), 0.05)
        self.assertAlmostEqual(fl.largura_feixe(0.4, ms.C0 / 0.04, 90), np.degrees(0.1), places=6)   # λ0/L


class Sintonia(unittest.TestCase):
    def test_mems_pull_in(self):
        from lta.core import sintonia as st
        g0, k, area = 3e-6, 10.0, 1e-8
        g = np.linspace(0.01, 0.999, 20001) * g0
        v = st.mems_tensao(g, g0, k, area)
        self.assertAlmostEqual(g[np.argmax(v)] / g0, 2 / 3, delta=1e-3)                 # máximo em 2g0/3
        self.assertAlmostEqual(v.max() / st.mems_pull_in(g0, k, area), 1, places=6)      # eq. 5.5

    def test_varactor(self):
        from lta.core import sintonia as st
        c = st.varactor(np.array([0, 5, 30]), 9e-12)
        np.testing.assert_allclose(st.tensao_varactor(c, 9e-12), [0, 5, 30], atol=1e-9)
        self.assertTrue(np.all(np.diff(c) < 0))
        self.assertAlmostEqual(st.hbv(1.3, 0.12e-12, 0.5), st.hbv(-1.3, 0.12e-12, 0.5))   # HBV simétrico

    def test_vlsrr_fig_5_7(self):
        from lta.core import sintonia as st
        f = np.linspace(1e9, 5e9, 8001)
        f0s, rej = [], []
        for cv in (1e-12, 9e-12):
            ceq = st.vlsrr_ceq(1.4e-12, cv)
            ls, cs, r, w0 = st.vlsrr_tanque(1.2e-9, 3e-9, ceq, 4.8)
            self.assertAlmostEqual(ls, 4 * 1.2e-9**2 / 3e-9, delta=1e-18)                # Ls' = 4M²/Lr
            self.assertAlmostEqual(w0 * r * cs, w0 * 3e-9 / 4.8, places=6)                # Q = ω0·Lr/Rs
            s = abcd.abcd_para_s(abcd.potencia(st.vlsrr_abcd(f, 3.3e-9, 1.3e-12, 1.2e-9, 3e-9, ceq, 4.8), 2), 50)
            i = np.argmin(abs(s[:, 1, 0]))
            self.assertAlmostEqual(f[i] / (w0 / 2 / np.pi), 1, delta=0.005)               # notch em f0
            rej.append(abcd.db(s[i, 1, 0]))
        self.assertLess(rej[0], rej[1])          # rejeição mais funda no topo da faixa (p. 373)

    def test_crlh_sintonizada(self):
        from lta.core import sintonia as st, crlh
        c = crlh.balanceada_por_banda(1e9, 6e9, 50)
        f = np.linspace(0.5e9, 6e9, 50)
        np.testing.assert_allclose(st.crlh_abcd_q(c, f), c.abcd(f))
        for x in (0.6, 1.7):
            cx = st.crlh_sintonizada(c, x)
            self.assertTrue(cx.balanceada)
            f0 = crlh.frequencia_transicao(cx)
            self.assertAlmostEqual(f0 / (crlh.frequencia_transicao(c) / np.sqrt(x)), 1, places=9)
            self.assertAlmostEqual(abs(cx.zb(np.array([f0 * 1.0001]))[0]), 50, delta=0.05)  # eq. 5.11
        self.assertAlmostEqual(st.figura_merito([10, 100], [-0.5, -2.0]), 45.0)             # eq. 5.12


class NLTL(unittest.TestCase):
    def _linha(self, cm=0.12):
        from lta.core import nltl as nl, sintonia as st
        return nl.LinhaNL(108.7e-12, 43.5e-15, lambda v: st.hbv(v, cm * 1e-12, 0.5))   # livro, Fig. 5.40

    def test_atraso_x_amplitude(self):
        from lta.core import nltl as nl
        lin = self._linha()
        tmin, tmax = nl.faixa_atraso(lin)
        self.assertAlmostEqual(tmin * 1e12, 2.174, delta=0.002)
        self.assertAlmostEqual(tmax * 1e12, 4.216, delta=0.002)
        amps = [nl.amplitude_soliton(lin, t) for t in (3.5e-12, 3.75e-12, 4e-12)]       # Fig. 5.41
        self.assertTrue(amps[0] > amps[1] > amps[2] > 0)
        for a, t in zip(amps, (3.5e-12, 3.75e-12, 4e-12)):
            self.assertAlmostEqual(nl.atraso_soliton(lin, a) / t, 1, places=6)

    def test_velocidade_simulada(self):
        # sóliton que emerge de um pulso: atraso medido no tempo contra a eq. 5.24
        from lta.core import nltl as nl
        lin = self._linha()
        fonte = lambda t: 2 * 3 * np.exp(-((t - 40e-12) / 8e-12) ** 2)                  # noqa: E731
        t, _, reg = nl.simular(lin, [120], fonte, 700e-12, 0.02e-12, nos=[60, 110])
        i1, i2 = np.argmax(reg[:, 0]), np.argmax(reg[:, 1])
        medido = (t[i2] - t[i1]) / 50
        self.assertAlmostEqual(medido / nl.atraso_soliton(lin, reg[i2, 1]), 1, delta=0.01)

    def test_multiplicador_tendencia(self):
        # livro, Fig. 5.44: CM maior -> máximo do 3º harmônico com menos seções
        from lta.core import nltl as nl
        ns = list(range(9, 46, 3))
        otimo = []
        for cm in (0.12, 0.24):
            lin = self._linha(cm)
            t_fim = 2 * max(ns) * nl.faixa_atraso(lin)[1] + 12 / 20e9                      # como na aba
            t, sai, _ = nl.simular(lin, ns, lambda t: -20 * np.sin(2 * np.pi * 20e9 * t), t_fim, 0.03e-12)
            ef = nl.eficiencia_harmonica(t, sai, 20e9, 3, lin.z0_linha, 100 / (2 * lin.z0_linha))
            otimo.append(ns[int(np.argmax(ef))])
            self.assertTrue(0.1 < ef.max() < 0.6)
        self.assertGreater(otimo[0], otimo[1])
        self.assertTrue(21 <= otimo[0] <= 39)                                             # livro: 27


class Diferencial(unittest.TestCase):
    def test_tabela_6_1(self):
        from lta.core import diferencial as df
        n, p = 1e-9, 1e-12
        for nome, (f0, le, ce, lc, cc, cr, fbw) in df.TABELA_6_1.items():
            el = (le * n, ce * p, lc * n, cc * p, cr * p)
            self.assertAlmostEqual(df.zero_comum(*el[:4]) / (f0 * 1e9), 1, delta=0.06, msg=nome)   # eq. 6.11
            _, _, b = df.banda_rejeicao(*el)
            self.assertAlmostEqual(100 * b, fbw, delta=2.1, msg=nome)        # Tabela 6.1 (máx. 2,03 pontos)

    def test_filtro_converge_e_modo_diferencial(self):
        from lta.core import diferencial as df
        n, p = 1e-9, 1e-12
        el = (9.0 * n, 0.8 * p, 0.45 * n, 4.8 * p, 0.09 * p)
        fa, fb, fbw = df.banda_rejeicao(*el)
        f = np.linspace(2.4e9, 3.6e9, 2401)
        fz = df.zero_comum(*el[:4])
        bw = []
        for k in (3, 10):
            rej = abcd.db(df.s_comum(f, k, *el)[:, 1, 0]) < -20
            i = np.argmin(abs(f - fz))
            a, b = i, i
            while rej[a - 1]:
                a -= 1
            while rej[b + 1]:
                b += 1
            bw.append((f[b] - f[a]) / np.sqrt(f[a] * f[b]))
        self.assertLess(bw[0], bw[1])
        self.assertAlmostEqual(bw[1], fbw, delta=0.01)            # muitas células -> banda máxima
        lo, co = df.linha_impar(50, 6.8, 7.7)
        sdd = df.s_diferencial(f, 10, lo, co)
        self.assertGreater(abcd.db(sdd[(f > fa) & (f < fb), 1, 0]).min(), -0.5)   # diferencial passa


class MIW(unittest.TestCase):
    def test_dispersao_circuito_e_dualidade(self):
        from lta.core import miw
        f0, k, l_res = 2.45e9, -0.0256, 10e-9
        c_res = 1 / ((2 * np.pi * f0) ** 2 * l_res)
        f = np.linspace(2.3e9, 2.6e9, 3001)
        al, bl = abcd.dispersao(miw.abcd_miw(f, l_res, c_res, k * l_res))
        ok = al < 1e-9
        np.testing.assert_allclose(miw.frequencia(bl[ok], f0, k), f[ok], rtol=1e-9)      # eq. 6.1
        fa, fb, fbw = miw.bordas(f0, k)
        self.assertAlmostEqual(100 * fbw, 5.13, delta=0.01)                              # livro, p. 429
        kc = miw.fator_acoplamento(cm_sobre_cc=0.0256)                                    # EIW dual (eq. 6.4)
        self.assertEqual(miw.bordas(f0, kc)[:2], (fa, fb))
        self.assertGreater(miw.frequencia(0.1, f0, k), miw.frequencia(3.0, f0, k))       # regressiva

    def test_atraso_e_ecos(self):
        from lta.core import miw
        f0, k, l_res = 2.45e9, -0.0256, 10e-9
        c_res = 1 / ((2 * np.pi * f0) ** 2 * l_res)
        fz = np.linspace(2.40e9, 2.50e9, 20001)
        cel = miw.abcd_miw(fz, l_res, c_res, k * l_res)
        zb = abs(abcd.impedancia_bloch(miw.abcd_miw(np.array([f0]), l_res, c_res, k * l_res))[0][0])
        s = abcd.abcd_para_s(abcd.potencia(cel, 6), zb)
        tg = abcd.atraso_grupo(fz, s[:, 1, 0])[np.argmin(abs(fz - f0))]
        self.assertAlmostEqual(tg / (6 / (2 * np.pi * f0 * abs(k))), 1, delta=1e-3)     # eq. 6.9 (ida)
        f = np.linspace(1e6, 5e9, 80000)
        t, env = miw.resposta_pulso(f, miw.s11_etiqueta(f, l_res, c_res, k * l_res, 10, "11", 1.0, zb), f0, 15e-9)
        t_bit = 10 * 2 / (2 * np.pi * f0 * abs(k))
        for m in (1, 2):
            janela = (t > (m - 0.3) * t_bit) & (t < (m + 0.3) * t_bit)
            tp = t[janela][np.argmax(env[janela])]
            self.assertAlmostEqual(tp / (m * t_bit), 1, delta=0.08)                       # ecos dos bits


class SIW(unittest.TestCase):
    def test_projeto(self):
        from lta.core import siw
        p = siw.projeto(6e9, 3.55, 0.8e-3)
        self.assertAlmostEqual(p.a, ms.C0 / (2 * 6e9 * np.sqrt(3.55)), places=12)        # eq. 6.34
        a_dw = p.w - p.d**2 / (0.95 * p.s)                                                  # Deslandes e Wu
        self.assertAlmostEqual(a_dw / p.a, 1, delta=0.01)
        self.assertAlmostEqual(p.d_max, 2 * p.a / (5 * np.sqrt(3)), places=12)             # eq. 6.32
        self.assertEqual(siw.verificacoes(p), [])
        self.assertTrue(siw.verificacoes(siw.projeto(6e9, 3.55, 0.8e-3, s=2.0e-3)))      # s > 2d
        f = np.array([5e9, 9e9])
        s = abcd.abcd_para_s(siw.abcd_siw(f, p, 0.02), abs(p.z_te(np.array([9e9]))[0]))
        self.assertLess(abcd.db(s[0, 1, 0]), -10)                                           # abaixo do corte
        self.assertAlmostEqual(abcd.db(s[1, 1, 0]), 0, places=9)
        self.assertAlmostEqual(p.lambda_g(9e9), 2 * np.pi / np.sqrt((2 * np.pi * 9e9 * np.sqrt(3.55) / ms.C0) ** 2
                                                                    - (np.pi / p.a) ** 2), places=12)


class NaoFoster(unittest.TestCase):
    def test_linha_enz(self):
        from lta.core import naofoster as nf
        lin = nf.LinhaENZ(31.5, 106e-12, 6.1e-12, 0.1)
        self.assertAlmostEqual(lin.c_ef * 1e12, 45.0, places=6)                            # Fig. 6.32
        self.assertAlmostEqual(lin.z_ef, 48.3, delta=0.1)
        self.assertGreater(lin.v_sobre_c, 1)
        f = np.linspace(1e6, 200e6, 2001)
        s = abcd.abcd_para_s(lin.abcd(f, 10), lin.z_ef)
        tg = abcd.atraso_grupo(f, s[:, 1, 0])
        np.testing.assert_allclose(tg[[100, 1000]] * ms.C0 / 1.0, 1 / lin.v_sobre_c, rtol=0.02)   # superluminal
        self.assertFalse(nf.LinhaENZ(31.5, 106e-12, 12e-12, 0.1).estavel)                  # C' < 0
        self.assertAlmostEqual(nf.capacitor_para_eps(31.5, 106e-12, 0.1, lin.eps_ef), 6.1e-12, delta=1e-18)


class LinhasConvencionais(unittest.TestCase):
    def test_stripline_pozar(self):
        # Pozar, Microwave Engineering, ex. 3.5: b = 3,16 mm, εr = 2,2, 50 Ω -> W ≈ 2,66 mm (fórmula ~1%)
        from lta.core import linhas_planares as lp
        w = lp.sintetizar("Stripline", 50, Substrato(2.2, 1.58), b=3.16)
        self.assertAlmostEqual(w, 2.66, delta=0.06)

    def test_limites_coplanares(self):
        from lta.core import linhas_planares as lp
        self.assertAlmostEqual(lp.cpw(1, 0.5, 1e5, 1.0)[0], 30 * np.pi * 1.27920, delta=0.05)   # CPW no ar, k = 0,5
        self.assertAlmostEqual(lp.cpw(0.5, 0.2, 1000, 10.2)[1], 5.6, places=4)                  # εeff → (εr+1)/2
        np.testing.assert_allclose(lp.gcpw(0.5, 0.2, 100, 10.2), lp.cpw(0.5, 0.2, 100, 10.2), rtol=1e-4)
        z1, e = lp.cpw(1, 0.5, 1e5, 4.0)
        z2, _ = lp.cps(1, 0.5, 1e5, 4.0)
        self.assertAlmostEqual(z1 * z2 / (ms.ETA0**2 / (4 * e)), 1, delta=0.002)              # Booker (duais)

    def test_sintese_ida_e_volta(self):
        from lta.core import linhas_planares as lp
        sub = Substrato(3.55, 0.813)
        for tipo, z, g, b in (("Microstrip", 50, 0, 0), ("Stripline", 50, 0, 1.6), ("CPW", 50, 0.2, 0),
                              ("CPW com terra (GCPW)", 50, 0.2, 0), ("Coplanar strips (CPS)", 100, 0.2, 0)):
            w = lp.sintetizar(tipo, z, sub, g=g, b=b)
            self.assertAlmostEqual(lp.analisar(lp.Geometria(tipo, w, g, b), sub)[0], z, places=4, msg=tipo)

    def test_elementos(self):
        from lta.core import linhas_planares as lp
        sub = Substrato(3.55, 0.813, tand=0.0)
        geo = lp.Geometria("Microstrip", lp.sintetizar("Microstrip", 50, sub))
        f = np.array([2.45e9])
        comp = lp.comprimento_eletrico(geo, sub, 2.45e9, 90)
        s_ab = abcd.abcd_para_s(lp.abcd_stub(f, geo, sub, comp, perdas=False), 50)
        s_cc = abcd.abcd_para_s(lp.abcd_stub(f, geo, sub, comp, curto=True, perdas=False), 50)
        self.assertLess(abcd.db(s_ab[0, 1, 0]), -60)                  # stub λ/4 aberto: curto em f0
        self.assertAlmostEqual(abs(s_cc[0, 1, 0]), 1, places=9)       # stub λ/4 em curto: aberto em f0
        geo_t = lp.Geometria("Microstrip", lp.sintetizar("Microstrip", np.sqrt(50 * 100), sub))
        m = lp.abcd_trecho(f, geo_t, sub, lp.comprimento_eletrico(geo_t, sub, 2.45e9, 90), perdas=False)
        self.assertLess(abcd.db(abcd.abcd_para_s2(m, 50, 100)[0, 0, 0]), -60)   # transformador λ/4 casado
        mm = lp.abcd_trecho(np.linspace(1e9, 3e9, 5), geo, sub, 10)
        np.testing.assert_allclose(abcd.abcd_para_s2(mm, 50, 50), abcd.abcd_para_s(mm, 50), atol=1e-12)

    def test_layout_sem_terra(self):
        from lta.geometria import construtores as k
        sub = Substrato(3.55, 0.813)
        with tempfile.TemporaryDirectory() as d:
            for tipo, terra in (("CPW", False), ("CPW com terra (GCPW)", True), ("Coplanar strips (CPS)", False)):
                lay = k.linha_convencional(tipo, "Trecho de linha", dict(w=1.2, g=0.2, comp=20, largura_terra=3,
                                                                          raio_via=0.15, passo_vias=1.2), sub,
                                           nome=f"l{int(terra)}{tipo[0]}")
                arqs = exportar(lay, d, list(FORMATOS), dict(f_min_ghz=1, f_max_ghz=10))
                self.assertEqual(any(a.name.endswith("_terra.gbl") for a in arqs), terra, tipo)
                dxf = (Path(d) / f"{lay.nome}.dxf").read_text()
                self.assertEqual("POLYLINE\n8\nTERRA\n" in dxf, terra, tipo)


class Exportacao(unittest.TestCase):
    def test_todos_formatos(self):
        p = cm.perfil("Senoidal", 64.0, 0.01, 10)
        w, _ = cm.larguras(p, RO3010)
        lays = [construtores.linha_largura_modulada(p.z * 1e3, w, RO3010, 5.0),
                construtores.linha_furos_terra(18.9, 4.725, 7, 50, RO3010),
                construtores.linha_carga_capacitiva(
                    cc.geometria(cc.projeto_por_fase(35, 1.4e9, 47.5, 22e9), ARLON, 10, 50), 2, ARLON)]
        from lta.core.csrr import GeometriaCSRR
        for serie in ("Gap simples", "Interdigital"):
            lays.append(construtores.linha_csrr(2, 8.0, 1.19, GeometriaCSRR(3.6, 0.3, 0.3), RO3010, serie,
                                                nome=f"csrr_{serie[0]}"))
        from lta.core import cl_loaded, crlh
        geo_cl, _ = cl_loaded.projeto(crlh.CelulaCRLH(2.45e-9, 0.5e-12, 3.38e-9, 0.68e-12), ARLON, n_dedos=14)
        lays.append(construtores.linha_cl_loaded(geo_cl, 3, ARLON))
        from lta.core import siw
        lays.append(construtores.linha_siw(siw.projeto(10e9, 3.55, 0.6e-3), 20, RO3010, 1.2, 3.0, 4.0))
        with tempfile.TemporaryDirectory() as d:
            for lay in lays:
                arqs = exportar(lay, d, list(FORMATOS), dict(f_min_ghz=1, f_max_ghz=10))
                for a in arqs:
                    self.assertGreater(a.stat().st_size, 0, a)
                txt = (Path(d) / f"{lay.nome}.kicad_pcb").read_text()
                self.assertEqual(txt.count("("), txt.count(")"))
                compile((Path(d) / f"{lay.nome}_hfss.py").read_text(), "hfss", "exec")
                self._gds_valido((Path(d) / f"{lay.nome}.gds").read_bytes())
                if lay.vias:
                    drl = (Path(d) / f"{lay.nome}_furos.drl").read_text()
                    self.assertEqual(sum(1 for x in drl.splitlines() if x.startswith("X")), len(lay.vias))
                    self.assertEqual(txt.count("(via "), len(lay.vias))
            f = np.linspace(1e9, 2e9, 3)
            s = np.tile(np.eye(2), (3, 1, 1)).astype(complex)
            simuladores.touchstone(Path(d) / "t.s2p", f, s)
            self.assertEqual(len((Path(d) / "t.s2p").read_text().splitlines()), 5)

    def _gds_valido(self, dados):
        i, tipos = 0, []
        while i < len(dados):
            n, t = struct.unpack(">HH", dados[i:i + 4])
            self.assertGreaterEqual(n, 4)
            tipos.append(t)
            i += n
        self.assertEqual(i, len(dados))
        self.assertEqual(tipos[0], 0x0002)
        self.assertEqual(tipos[-1], 0x0400)


if __name__ == "__main__":
    unittest.main()
