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


# ---------------------------------------------------------------------------
# D. M. Pozar, Microwave Engineering, 4ª ed. (2012): exemplos resolvidos
# ---------------------------------------------------------------------------
class PozarMicrostrip(unittest.TestCase):
    def test_ex_3_7_formulas_fechadas(self):
        from lta.core import descontinuidades as dc
        sub = Substrato(9.9, 0.5)
        w = dc.pozar_sintese(50, sub)
        self.assertAlmostEqual(w, 0.483, delta=0.002)                     # W = 0,483 mm
        self.assertAlmostEqual(float(dc.pozar_analise(w, sub)[1]), 6.665, delta=0.005)
        self.assertAlmostEqual(float(analise(w, sub)[0]), 50, delta=1.0)   # Hammerstad-Jensen concorda

    def test_ex_3_4_ondas_de_superficie(self):
        from lta.core import descontinuidades as dc
        lim = dc.limites_substrato(Substrato(2.55, 1.0))
        self.assertAlmostEqual(lim["TE1 (onda de superfície)"] * 1e-3 / ms.C0, 1 / (4 * np.sqrt(1.55)), places=9)
        self.assertAlmostEqual(lim["TM1 (onda de superfície)"] / lim["TE1 (onda de superfície)"], 2, places=9)

    def test_ponta_aberta(self):
        # Kirschning-Jansen-Koster contra a fórmula mais simples de Hammerstad (diferem até ~10%)
        from lta.core import descontinuidades as dc
        for er, u in ((2.2, 3.0), (10.2, 1.0), (4.4, 0.5)):
            sub = Substrato(er, 1.0)
            ee = float(analise(u, sub)[1])
            ham = 0.412 * (ee + 0.3) * (u + 0.264) / ((ee - 0.258) * (u + 0.8))
            self.assertAlmostEqual(dc.extensao_ponta_aberta(u, sub) / ham, 1, delta=0.1, msg=(er, u))

    def test_gap(self):
        from lta.core import descontinuidades as dc
        sub = Substrato(9.6, 0.635)
        cg = [dc.gap(0.635, s, sub)[0] for s in (0.07, 0.15, 0.3, 0.6)]
        self.assertTrue(all(a > b > 0 for a, b in zip(cg, cg[1:])))       # cai com o gap
        s = dc.gap_para_capacitancia(0.635, cg[1], sub)
        self.assertAlmostEqual(s, 0.15, delta=1e-4)


class PozarCasamento(unittest.TestCase):
    def test_ex_5_1_rede_l(self):
        from lta.core import casamento as cs
        sol = cs.rede_l(200 - 100j, 100, 500e6)
        self.assertAlmostEqual(sol[0].b * 100, 0.29, delta=0.005)
        self.assertAlmostEqual(sol[0].x / 100, 1.22, delta=0.01)
        self.assertEqual(sol[0].derivacao[0], "C")
        self.assertAlmostEqual(sol[0].derivacao[1] * 1e12, 0.92, delta=0.01)
        self.assertAlmostEqual(sol[1].derivacao[1] * 1e9, 46.1, delta=0.2)    # 2ª solução: L = 46,1 nH
        for s in sol:
            g = cs.gama(cs.zin_abcd(s.abcd(np.array([500e6])), 200 - 100j), 100)
            self.assertLess(abs(g[0]), 1e-9)

    def test_ex_5_2_e_5_3_stub_simples(self):
        from lta.core import casamento as cs
        s = cs.stub_simples(60 - 80j, 50, serie=False, curto=True)
        np.testing.assert_allclose([(x.d, x.l) for x in s], [(0.110, 0.095), (0.260, 0.405)], atol=0.002)
        s = cs.stub_simples(100 + 80j, 50, serie=True, curto=False)
        np.testing.assert_allclose([(x.d, x.l) for x in s], [(0.120, 0.397), (0.463, 0.103)], atol=0.002)
        zl = cs.carga_serie(60 - 80j, 2e9)
        for x in cs.stub_simples(60 - 80j, 50, curto=True):
            m = cs.abcd_stub_simples(np.array([2e9]), 2e9, 50, x, True, False)
            self.assertLess(abs(cs.gama(cs.zin_abcd(m, zl(np.array([2e9]))), 50)[0]), 1e-9)

    def test_ex_5_4_stub_duplo(self):
        from lta.core import casamento as cs
        s = cs.stub_duplo(60 - 80j, 50, 0.125, curto=False)
        np.testing.assert_allclose([(x.b1, x.b2, x.l1, x.l2) for x in s],
                                   [(1.314, 3.38, 0.146, 0.204), (-0.114, -1.38, 0.482, 0.350)], atol=0.003)
        with self.assertRaises(ValueError):                              # região proibida (g > 2 para d = λ/8)
            cs.stub_duplo(20 + 0j, 50, 0.125)

    def test_ex_5_6_e_5_7_multissecao(self):
        from lta.core import casamento as cs
        np.testing.assert_allclose(cs.binomial(100, 50, 3), [91.7, 70.7, 54.5], atol=0.06)
        self.assertAlmostEqual(cs.banda_binomial(100, 50, 3, 0.05), 0.70, delta=0.005)
        z, tm = cs.chebyshev(50, 100, 3, 0.05)
        np.testing.assert_allclose(z, [57.5, 70.7, 87.0], atol=0.06)
        self.assertAlmostEqual(np.degrees(tm), 44.7, delta=0.05)
        # Tabela 5.2 (exata), ZL/Z0 = 2, Γm = 0,05
        np.testing.assert_allclose(np.array(cs.chebyshev(1, 2, 2, 0.05)[0]), [1.2193, 1.6402], rtol=0.003)
        np.testing.assert_allclose(np.array(cs.chebyshev(1, 2, 4, 0.05)[0]), [1.1201, 1.2979, 1.5409, 1.7855],
                                   rtol=0.003)
        f = np.linspace(0.65, 1.35, 141) * 1e9                         # banda de 70% (θm = 44,7°)
        g = cs.gama(cs.zin_abcd(cs.abcd_multissecao(f, 1e9, z), 100.0), 50)
        self.assertLess(np.abs(g).max(), 0.055)

    def test_ex_5_8_klopfenstein(self):
        from lta.core import casamento as cs
        self.assertAlmostEqual(cs.parametro_a(100, 50, 0.02), 3.543, delta=0.002)
        s, z = cs.perfil_taper("Klopfenstein", 100, 50, 401, 0.02)
        self.assertAlmostEqual(np.log(z[0] / 100), -0.02, delta=1e-4)  # degrau Γ0/cosh A = -Γm na ponta
        np.testing.assert_allclose(cs.phi_klopfenstein([1.0], 2.0), (np.cosh(2) - 1) / 4, rtol=1e-6)
        # resposta exata (cascata) na banda βL ≥ A: |Γ| ≈ Γm
        bl = np.linspace(1.1 * 3.543, 12, 60)
        zm = 0.5 * (z[1:] + z[:-1])
        g = []
        for b in bl:
            m = abcd.cascata(*[abcd.linha(zz, np.array([b / 400])) for zz in zm])
            g.append(abs(cs.gama(cs.zin_abcd(m, 50.0), 100)[0]))
        self.assertLess(max(g), 0.024)
        self.assertLess(cs.gama_taper_analitico("Klopfenstein", bl, 100, 50, 0.02).max(), 0.0201)

    def test_bode_fano(self):
        from lta.core import casamento as cs
        b = cs.bode_fano(60 - 80j, 2e9, 1e9, 3e9)
        r, c = 60, 1 / (2 * np.pi * 2e9 * 80)
        w1, w2 = 2 * np.pi * 1e9, 2 * np.pi * 3e9
        self.assertAlmostEqual(np.log(1 / b.gama_min) * (1 / w1 - 1 / w2), np.pi * r * c, places=12)


class PozarAcopladas(unittest.TestCase):
    def test_ex_7_7_stripline(self):
        from lta.core import linhas_acopladas as la
        ze, zo = la.impedancias_acoplador(50, 0.1)
        self.assertAlmostEqual(ze, 55.28, delta=0.01)
        self.assertAlmostEqual(zo, 45.23, delta=0.01)
        w, s = la.sintetizar("Stripline", ze, zo, Substrato(2.2, 1.6), b=3.2)
        self.assertAlmostEqual(w, 2.59, delta=0.02)                     # W = 0,259 cm
        self.assertAlmostEqual(s, 0.98, delta=0.01)                     # S = 0,098 cm

    def test_ex_7_8_binomial(self):
        from lta.core import linhas_acopladas as la
        c = la.coeficientes_binomial(0.1, 3)
        np.testing.assert_allclose(c, [0.0125, 0.125, 0.0125], atol=1e-12)
        np.testing.assert_allclose(la.impedancias_acoplador(50, c[0]), (50.63, 49.38), atol=0.01)
        np.testing.assert_allclose(la.impedancias_acoplador(50, c[1]), (56.69, 44.10), atol=0.01)
        # TEM ideal: 20 dB em f0, isolamento perfeito e casamento em toda a banda
        secoes = [la.Secao(la.Modos(*la.impedancias_acoplador(50, x), 1.0, 1.0), ms.C0 / 4 / 3e9 * 1e3) for x in c]
        f = np.linspace(1e9, 5e9, 41)
        s = la.acoplador(f, secoes)
        self.assertAlmostEqual(-abcd.db(s[20, 2, 0]), 20.0, delta=0.05)          # teoria de acoplamento fraco
        self.assertLess(np.abs(s[:, 3, 0]).max(), 1e-9)
        self.assertLess(np.abs(s[:, 0, 0]).max(), 1e-9)

    def test_lange(self):
        from lta.core import linhas_acopladas as la
        c = 10 ** (-3 / 20)
        ze, zo = la.impedancias_lange(50, c)
        e4, o4 = la.modos_lange(ze, zo)
        self.assertAlmostEqual(np.sqrt(e4 * o4), 50, places=6)          # eq. 7.98
        self.assertAlmostEqual((e4 - o4) / (e4 + o4), c, places=6)      # eq. 7.99

    def test_microstrip_garg_bahl(self):
        from lta.core import linhas_acopladas as la
        sub = Substrato(10, 1.0)
        m = la.microstrip_acoplada(1.0, 0.5, sub)
        self.assertGreater(m.ze, m.zo)
        self.assertGreater(m.eee, m.eeo)                                # modo par mais lento (Pozar, p. 375)
        self.assertAlmostEqual(m.ze, 60, delta=3)                       # Fig. 7.30 (εr = 10, W/h = 1, S/h = 0,5)
        self.assertAlmostEqual(m.zo, 38, delta=3)
        w, s = la.sintetizar("Microstrip", 70, 36, sub, modelo="Garg-Bahl")
        r = la.microstrip_acoplada(w, s, sub)
        self.assertAlmostEqual(r.ze, 70, delta=0.01)
        self.assertAlmostEqual(r.zo, 36, delta=0.01)
        w, s = la.sintetizar("Microstrip", 70, 36, sub)                 # padrão: Kirschning-Jansen
        r = la.microstrip_acoplada_kj(w, s, sub)
        self.assertAlmostEqual(r.ze, 70, delta=0.01)
        self.assertAlmostEqual(r.zo, 36, delta=0.01)

    def test_lta_diretividade(self):
        # compensar o modo ímpar com capacitores periódicos recupera a diretividade da microstrip
        from lta.core import linhas_acopladas as la
        sub = Substrato(3.55, 0.813, 35, 0.0027)
        f = np.array([2e9])
        d = []
        for lta in (False, True):
            p = la.projeto_acoplador("Microstrip", 50, 10, 2e9, sub, 1, lta=lta, n_cel=8)
            c, di, _, _ = la.metricas(la.acoplador(f, p.secoes, 50, sub))
            self.assertAlmostEqual(c[0], 10, delta=0.15)
            d.append(di[0])
        self.assertLess(d[0], 20)
        self.assertGreater(d[1], 40)


class PozarFiltros(unittest.TestCase):
    def setUp(self):
        from lta.core import filtros
        self.g05 = filtros.prototipo("Chebyshev", 3, 0.5)

    def test_ex_8_5_kuroda(self):
        from lta.core import filtros, filtros_convencionais as fc
        g = filtros.prototipo("Chebyshev", 3, 3.0)
        np.testing.assert_allclose(g[1:4], [3.3487, 0.7117, 3.3487], atol=4e-4)   # Tabela 8.4 (4 casas)
        els = fc.richards_kuroda(g)
        np.testing.assert_allclose([e.z for e in els], [1.299, 4.350, 1.405, 4.350, 1.299], atol=0.002)
        f = np.linspace(0.3e9, 7e9, 25)
        ref = abcd.cascata(fc.abcd_ue(f, 4e9, 50, np.pi / 4), fc.abcd_richards(f, 4e9, g, 50),
                           fc.abcd_ue(f, 4e9, 50, np.pi / 4))
        np.testing.assert_allclose(abcd.abcd_para_s(fc.abcd_elementos(f, 4e9, els, 50), 50),
                                   abcd.abcd_para_s(ref, 50), atol=1e-9)      # identidade exata
        s = abcd.abcd_para_s(fc.abcd_elementos(np.array([4e9]), 4e9, els, 50), 50)
        self.assertAlmostEqual(abcd.db(s[0, 1, 0]), -3.0, delta=0.01)        # ripple de 3 dB em fc
        with self.assertRaises(ValueError):
            fc.richards_kuroda(filtros.prototipo("Butterworth", 5))

    def test_ex_8_6_degraus(self):
        from lta.core import filtros, filtros_convencionais as fc, linhas_planares as lp
        g = filtros.prototipo("Butterworth", 6)
        tr = fc.degraus(g, 50, 20, 120, comeca_com_c=True)
        np.testing.assert_allclose([np.degrees(t.bl) for t in tr], [11.8, 33.8, 44.3, 46.1, 32.4, 12.3], atol=0.15)
        sub = Substrato(4.2, 1.58)
        self.assertAlmostEqual(lp.sintetizar("Microstrip", 20, sub), 11.3, delta=0.25)
        self.assertAlmostEqual(lp.sintetizar("Microstrip", 120, sub), 0.428, delta=0.02)
        f = np.array([4e9])
        s = abcd.abcd_para_s(fc.abcd_degraus_ideal(f, 2.5e9, tr), 50)
        self.assertLess(abcd.db(s[0, 1, 0]), -20)                         # > 20 dB em 4 GHz

    def test_ex_8_7_linhas_acopladas(self):
        from lta.core import filtros_convencionais as fc
        p = fc.projeto_linhas_acopladas(self.g05, 2e9, 0.1, 50)
        np.testing.assert_allclose(p.jz, [0.3137, 0.1187, 0.1187, 0.3137], atol=1e-4)
        np.testing.assert_allclose(p.ze, [70.61, 56.64, 56.64, 70.61], atol=0.01)
        np.testing.assert_allclose(p.zo, [39.24, 44.77, 44.77, 39.24], atol=0.01)
        f = np.array([1.8e9, 2.0e9])
        s = fc.s_linhas_acopladas(f, fc.secoes_ideais(2e9, p.ze, p.zo))
        self.assertAlmostEqual(abcd.db(s[0, 1, 0]), -20, delta=1.5)        # ~20 dB em 1,8 GHz
        self.assertGreater(abcd.db(s[1, 1, 0]), -0.6)

    def test_lta_suprime_espuria(self):
        from lta.core import filtros_convencionais as fc
        sub = Substrato(3.55, 0.813, 35, 0.0)
        f = np.linspace(3.6e9, 4.4e9, 161)
        pico = []
        for lta in (False, True):
            p = fc.projeto_linhas_acopladas(self.g05, 2e9, 0.1, 50, sub, lta=lta, n_cel=6)
            s = fc.s_linhas_acopladas(f, p.secoes, 50, sub)
            pico.append(abcd.db(s[:, 1, 0]).max())
            s0 = fc.s_linhas_acopladas(np.array([2e9]), p.secoes, 50, sub)
            self.assertGreater(abcd.db(s0[0, 1, 0]), -1.0)
        self.assertGreater(pico[0], -3)                                  # espúria em 2f0 na microstrip
        self.assertLess(pico[1], -25)                                    # suprimida com LTA

    def test_ex_8_8_a_8_10(self):
        from lta.core import filtros_convencionais as fc
        np.testing.assert_allclose(fc.stubs_quarto(self.g05, 0.15, 50), [265.9, 387.0, 265.9], atol=0.1)
        p = fc.ressonadores_serie_gap(self.g05, 2e9, 0.1, 50)
        np.testing.assert_allclose(p.c * 1e12, [0.554, 0.192, 0.192, 0.554], atol=0.001)
        np.testing.assert_allclose(np.degrees(p.theta), [155.8, 166.5, 155.8], atol=0.1)
        q = fc.ressonadores_derivacao(self.g05, 2.5e9, 0.1, 50)
        np.testing.assert_allclose(q.c * 1e12, [0.2896, 0.0756, 0.0756, 0.2896], atol=2e-4)
        np.testing.assert_allclose(q.l * 360, [73.6, 83.2, 73.6], atol=0.1)


class PozarRessonadores(unittest.TestCase):
    def test_ex_6_2(self):
        from lta.core import descontinuidades as dc, ressonadores as rs
        sub = Substrato(2.08, 1.59, 35, 0.0004, "Teflon")
        r = rs.projeto(rs.TIPOS[0], 5e9, 50, sub, corrigir_ponta=False)
        self.assertAlmostEqual(r.w, 5.08, delta=0.03)
        self.assertAlmostEqual(r.eeff, 1.80, delta=0.01)
        self.assertAlmostEqual(r.comp, 22.4, delta=0.1)
        self.assertAlmostEqual(r.alfa_d, 0.024, delta=0.001)
        ac = float(dc.pozar_alfa_condutor(r.w, sub, 5e9))
        self.assertAlmostEqual(ac, 0.0724, delta=0.001)
        # O livro usa β = 151,0 rad/m, mas 2π·5 GHz·√1,80/c = 140,5 rad/m: com αc e αd do livro,
        # Q0 = 140,5/(2·0,0964) ≈ 729 (o 783 do texto herda o β inconsistente).
        self.assertAlmostEqual(r.beta, 140.5, delta=0.3)
        self.assertAlmostEqual(r.beta / (2 * (ac + r.alfa_d)), 729, delta=5)

    def test_ex_6_6_gap(self):
        from lta.core import ressonadores as rs
        alfa = 0.01 / ms.NP_DB * 100                                     # 0,01 dB/cm
        r = rs.Ressonador(rs.TIPOS[0], 5e9, 50, 4.0, 1.9, 21.75, alfa, 0.0)
        self.assertAlmostEqual(r.q0, 628, delta=2)
        a = rs.acoplamento_gap(r, Substrato(2.2, 1.6))
        self.assertAlmostEqual(a.bc, 0.05, delta=0.0005)
        self.assertAlmostEqual(a.c * 1e12, 0.032, delta=0.0005)
        self.assertAlmostEqual(a.f_carregada / 1e9, 4.918, delta=0.005)
        f = np.linspace(4.85e9, 5.0e9, 15001)
        z = rs.zin_acoplado(f, r, a)
        g = np.abs((z - 50) / (z + 50))
        self.assertLess(g.min(), 0.03)                                   # acoplamento crítico


class PozarDivisores(unittest.TestCase):
    def test_ex_7_1_juncao_t(self):
        from lta.core import dispositivos as dv
        z2, z3 = dv.impedancias_juncao_t(50, 0.5)                       # P3/P2 = 1/2 → 75 e 150 Ω
        self.assertAlmostEqual(z2, 75)
        self.assertAlmostEqual(z3, 150)
        s = dv.juncao_t(np.array([1e9]), 50, 150, 75)[0]
        self.assertAlmostEqual(abs(s[0, 0]), 0, places=5)
        self.assertAlmostEqual(s[1, 1].real, -0.666, delta=0.001)
        self.assertAlmostEqual(s[2, 2].real, -0.333, delta=0.001)

    def test_wilkinson_desigual_eq_7_37(self):
        from lta.core import dispositivos as dv
        f = np.array([1e9])
        z2, z3, r, r2, r3 = dv.impedancias_wilkinson(50, 2.0)
        lin = lambda z: dv.linha_convencional(z, np.pi / 2, 1e9)   # noqa: E731
        for tr in ((None, None), (lin(np.sqrt(50 * r2)), lin(np.sqrt(50 * r3)))):
            s = dv.wilkinson_desigual(f, lin(z2), lin(z3), 2.0, 50, *tr)[0]
            np.testing.assert_allclose(np.abs(s) ** 2, [[0, 1 / 3, 2 / 3], [1 / 3, 0, 0], [2 / 3, 0, 0]], atol=1e-9)

    def test_rede_exata_em_f0(self):
        # nós internos só com linhas λ/4 ideais exatamente em f0 (antes singular na redução de Kron)
        from lta.core import dispositivos as dv
        f = np.array([1e9])
        lin = dv.linha_convencional(50, np.pi / 2, 1e9)
        r = dv.Rede(5, [0, 3, 4], f)
        r.bipolo(0, 1, dv.linha_convencional(50 * np.sqrt(2), np.pi / 2, 1e9)(f))
        r.bipolo(0, 2, dv.linha_convencional(50 * np.sqrt(2), np.pi / 2, 1e9)(f))
        r.impedancia(1, 2, 100)
        r.bipolo(1, 3, lin(f))
        r.bipolo(2, 4, lin(f))
        np.testing.assert_allclose(np.abs(r.s(50)[0]) ** 2, [[0, .5, .5], [.5, 0, 0], [.5, 0, 0]], atol=1e-9)

    def test_slow_wave(self):
        from lta.core import dispositivos as dv
        sub = Substrato(3.55, 0.813, 35, 0.0)
        f = np.array([2e9])
        ramos = []
        for z in (50 / np.sqrt(2), 50):
            for tipo in dv.CARGAS_SW:
                try:
                    p = dv.projeto_slow_wave(z, np.pi / 2, 2e9, 3, 100, sub, 60 if tipo == dv.CARGAS_SW[0] else 10,
                                             tipo)
                except ValueError:
                    continue
                m = dv.linha_slow_wave(p, sub)(f)
                self.assertAlmostEqual(float(abcd.dispersao(m)[1][0]), np.pi / 2, places=6)
                self.assertAlmostEqual(abcd.impedancia_bloch(m)[0][0].real, z, places=4)
                self.assertLess(p.reducao, 0.6)                           # miniaturização
                if tipo == dv.CARGAS_SW[0]:
                    ramos.append(dv.linha_slow_wave(p, sub))
        s = dv.branch_line(f, *ramos)[0]
        np.testing.assert_allclose(np.abs(s[:, 0]) ** 2, [0, 0.5, 0.5, 0], atol=1e-9)


class LayoutsPozar(unittest.TestCase):
    def test_exportacao_circuitos(self):
        from lta.core import dispositivos as dv
        from lta.geometria import circuitos as ci
        from lta.geometria.layout import largura_porta_onda
        sub = Substrato(3.55, 0.813)
        p = dv.projeto_slow_wave(50 / np.sqrt(2), np.pi / 2, 2e9, 3, 100, sub, 60)
        q = dv.projeto_slow_wave(50, np.pi / 2, 2e9, 3, 100, sub, 60)
        lays = [ci.linha_com_stubs("cs", sub, [(10, 1.8)], [dict(x=0, w=1.8, comp=8, curto=True)], 1.8, 1.8, 5, 0),
                ci.linha_perfil("lp", sub, [0, 5, 5, 10], [1.8, 1.8, 1.0, 1.0], 1.8, 1.0),
                ci.acoplador_linhas("ac", sub, [(20, 1.5, 0.15)], 1.8, capacitores=[4]),
                ci.filtro_linhas_acopladas("fa", sub, [(22, 1.2, 0.2), (22, 1.6, 0.6), (22, 1.2, 0.2)], 1.8),
                ci.filtro_gaps("fg", sub, 1.8, [0.1, 0.2, 0.1], [40, 40]),
                ci.ressonador_gap("rg", sub, 1.8, 1.8, 40, 0.1, curto=True),
                ci.branch_line("bl", sub, 1.8, 3.0, 21, 1.8, 22),
                ci.branch_line("bs", sub, 1.8, 0, 0, 0, 0, sw_s=(p.geo, p.n), sw_p=(q.geo, q.n)),
                ci.rat_race("rr", sub, 1.8, 1.0, 21),
                ci.wilkinson("wk", sub, 1.8, 1.0, 22, 1.4, 22, 1.2, 2.5, t2=(1.4, 22), t3=(2.6, 21)),
                ci.juncao_t("jt", sub, 1.8, 0.8, 0.3, t2=(1.2, 22), t3=(0.7, 23))]
        with tempfile.TemporaryDirectory() as d:
            for lay in lays:
                arqs = exportar(lay, d, list(FORMATOS), dict(f_min_ghz=1, f_max_ghz=5))
                self.assertTrue(all(a.stat().st_size > 0 for a in arqs), lay.nome)
                compile((Path(d) / f"{lay.nome}_hfss.py").read_text(encoding="utf-8"), "hfss", "exec")
                txt = (Path(d) / f"{lay.nome}.kicad_pcb").read_text()
                self.assertEqual(txt.count("("), txt.count(")"))
                x0, y0, x1, y1 = lay.placa
                for pt in lay.portas:
                    self.assertIn(pt.x, (x0, x1), lay.nome)                # portas nas faces x
                    self.assertTrue(y0 < pt.y < y1, lay.nome)
                    wp = largura_porta_onda(lay, pt)
                    for o in lay.portas:                                 # portas de onda sem sobreposição
                        if o is not pt and o.lado == pt.lado:
                            self.assertLess(wp / 2 + largura_porta_onda(lay, o) / 2, abs(o.y - pt.y) + 1e-9,
                                            lay.nome)


# ---------------------------------------------------------------------------
# R. Garg, I. Bahl e M. Bianchi, Microstrip Lines and Slotlines, 4ª ed. (2024)
# ---------------------------------------------------------------------------
class GargBahlMicrostrip(unittest.TestCase):
    def test_espessura(self):
        # eqs. 2.124-2.127: t > 0 baixa Z0 e εeff (Tabela 2.5: Δεeff ≈ -0,03 a -0,08 em alumina)
        al = Substrato(9.7, 0.635, 10)
        w = sintese(50, al)
        z0, e0 = analise(w, al)
        zt, et = analise(w, Substrato(9.7, 0.635, 10, incluir_t=True))
        self.assertLess(zt, z0)
        self.assertTrue(-0.08 < et - e0 < -0.02)
        self.assertAlmostEqual(float(analise(w, Substrato(9.7, 0.635, 0.0, incluir_t=True))[0]), float(z0), places=9)

    def test_fT_qr_e_via(self):
        from lta.core import descontinuidades as dc
        self.assertAlmostEqual(dc.frequencia_ft(Substrato(3.8, 0.5)) / 1e9, 106.0, delta=0.5)   # < 0,5 mm a 100 GHz
        self.assertAlmostEqual(float(dc.q_radiacao(50, Substrato(3.8, 1.0), 1e9)), 4252, delta=10)  # 0,4252e4/(fh)²
        lv, _ = dc.via(0.1, Substrato(9.6, 0.38), 1e9)
        self.assertAlmostEqual(lv * 1e12, 67, delta=1)                                      # h = 380, r = 100 µm

    def test_curva_e_chanfro(self):
        from lta.core import descontinuidades as dc
        al = Substrato(9.7, 0.635)
        w = sintese(50, al)
        self.assertAlmostEqual(dc.chanfro_otimo(w, al), 72, delta=4)        # ~72% em alumina (Easter)
        cb, _ = dc.curva(w, al)
        self.assertGreater(cb, 0)
        self.assertGreater(dc.encurtamento_chanfro(w, al, 1e9), 0)

    def test_degrau_corrige_filtro(self):
        from lta.core import filtros, filtros_convencionais as fc
        sub = Substrato(4.2, 1.58, 0, 0)
        g = filtros.prototipo("Butterworth", 6)
        tr = fc.degraus(g, 50, 20, 120)
        segs = [(t.bl / (2 * np.pi) * 299.792458 / 2.5 / np.sqrt(float(analise(w, sub)[1])), w)
                for t, w in zip(tr, [sintese(t.z, sub) for t in tr])]
        cor, dls = fc.correcao_degraus(segs, sub)
        self.assertEqual(len(dls), 5)
        # o modelo com degraus aplicado aos trechos corrigidos volta aos comprimentos elétricos de projeto
        f = np.array([1.0e9, 2.5e9])
        a = abcd.abcd_para_s(fc.abcd_degraus_microstrip(f, segs, sub, com_degraus=False), 50)
        b = abcd.abcd_para_s(fc.abcd_degraus_microstrip(f, cor, sub, com_degraus=True), 50)
        np.testing.assert_allclose(a, b, atol=1e-9)

    def test_coplanares_e_slotline(self):
        from lta.core import linhas_planares as lp
        a, b = Substrato(3.55, 0.813, 35), Substrato(3.55, 0.813, 35, incluir_t=True)
        for tp in ("CPW", "CPW com terra (GCPW)", "Coplanar strips (CPS)"):
            geo = lp.Geometria(tp, 0.6, 0.2)
            z0, e0 = lp.analisar(geo, a)
            zt, et = lp.analisar(geo, b)
            self.assertLess(zt, z0, tp)                                   # espessura baixa Z0 e εeff (seção 7.3.3)
            self.assertLess(et, e0, tp)
        # perda no condutor da CPW (eq. 7.126) na ordem das medidas da Fig. 7.23 (0,4-0,8 dB/cm)
        st = Substrato(12.9, 0.5, 1.61, 0, "", 3.602e7)
        ac = lp.alfa_condutor_coplanar(lp.Geometria("CPW", 0.0712, 0.049), st, np.array([10e9])) * ms.NP_DB / 100
        self.assertTrue(0.4 < ac[0] < 0.9)
        z, e = lp.slotline(0.3, 0.787, 2.22, 3e9)
        self.assertTrue(80 < z < 110 and 1.1 < e < 1.6)
        w = lp.sintetizar("Slotline", 90, Substrato(2.22, 0.787), f=3e9)
        self.assertAlmostEqual(lp.analisar(lp.Geometria("Slotline", w), Substrato(2.22, 0.787), 3e9)[0], 90, places=3)

    def test_stub_com_via(self):
        from lta.core import linhas_planares as lp
        sub = Substrato(3.55, 0.813, 35, 0.0)
        geo = lp.Geometria("Microstrip", sintese(50, sub))
        f0 = 2e9
        comp = lp.comprimento_eletrico(geo, sub, f0, 90)
        d = lp.encurtamento_via(0.2, geo, sub, f0)
        y = lp.abcd_stub(np.array([f0]), geo, sub, comp - d, curto=True, perdas=False, raio_via=0.2)[0, 1, 0]
        self.assertLess(abs(y) * 50, 0.02)                                # λ/4 em curto compensado: aberto em f0


class GargBahlAcopladas(unittest.TestCase):
    def test_kirschning_jansen_fig_4_12(self):
        from lta.core import linhas_acopladas as la
        sub = Substrato(9.7, 0.64)
        m1 = la.microstrip_acoplada_kj(0.64, 0.064, sub, dispersivo=True)      # S/h = 0,1
        m2 = la.microstrip_acoplada_kj(0.64, 0.64, sub, dispersivo=True)       # S/h = 1
        self.assertAlmostEqual(m1.eee, 7.0, delta=0.1)
        self.assertAlmostEqual(m1.eeo, 5.55, delta=0.1)
        self.assertAlmostEqual(m2.eeo, 5.8, delta=0.1)
        ee30, eo30 = m1.disp(np.array([30e9]))
        self.assertAlmostEqual(float(ee30[0]), 8.35, delta=0.15)
        self.assertAlmostEqual(float(eo30[0]), 6.05, delta=0.15)
        self.assertAlmostEqual(float(m2.disp(np.array([30e9]))[1][0]), 6.85, delta=0.15)
        gb = la.microstrip_acoplada(0.64, 0.32, sub)                           # Garg-Bahl a ~3%
        kj = la.microstrip_acoplada_kj(0.64, 0.32, sub)
        self.assertAlmostEqual(gb.ze / kj.ze, 1, delta=0.03)
        self.assertAlmostEqual(gb.zo / kj.zo, 1, delta=0.04)

    def test_espessura_perdas_ponta(self):
        from lta.core import linhas_acopladas as la
        sub = Substrato(9.7, 0.64, 10)
        m0 = la.microstrip_acoplada_kj(0.64, 0.128, sub)
        mt = la.microstrip_acoplada_kj(0.64, 0.128, Substrato(9.7, 0.64, 10, incluir_t=True))
        self.assertLess(mt.zo, m0.zo)
        self.assertGreater(m0.eeo - mt.eeo, m0.eee - mt.eee)              # ímpar cai mais (seção 4.5.1)
        s = Substrato(3.55, 0.813, 35, 0.0)
        f = np.array([5e9])
        ae, ao = la.alfa_condutor_modos(1.8, 0.2, s, f)
        self.assertGreater(ao[0], ae[0])                                  # αco > αce (Fig. 4.19)
        ae, ao = la.alfa_condutor_modos(1.8, 5.0, s, f)
        iso = ms.alfa_condutor(1.8, s, f)[0, 0]
        self.assertAlmostEqual(ae[0] / iso, 1, delta=0.05)               # linhas afastadas: linha isolada
        dle, dlo = la.ponta_aberta_acoplada(0.64, 6.4, sub)
        from lta.core.descontinuidades import extensao_ponta_aberta
        self.assertAlmostEqual(dle, extensao_ponta_aberta(0.64, sub), delta=1e-3)   # Δle → Δloc para S grande
        self.assertAlmostEqual(dlo, dle, delta=1e-3)

    def test_marchand(self):
        from lta.core import balun
        for zl in (50, 100, 200):
            c, r = balun.acoplamento_otimo(2e9, 50, zl)
            self.assertAlmostEqual(c, 1 / np.sqrt(1 + zl / 50), places=4)
            self.assertLess(r, 1e-4)
        p = balun.projeto(2e9, 50, 200)
        m = balun.metricas(balun.s_marchand(np.array([2e9]), *p.secoes, 50, 200))
        self.assertAlmostEqual(m[1][0], -3.01, delta=0.01)
        self.assertAlmostEqual(m[2][0], -3.01, delta=0.01)
        self.assertAlmostEqual(abs(m[4][0]), 180, delta=0.01)


class GargBahlPratico(unittest.TestCase):
    def test_conector_e_de_embedding(self):
        from lta.core import conectores as co
        t = co.PREDEFINIDOS["Típica, VSWR 1,10 em 10 GHz (Tabela 2.5)"]
        s = abcd.abcd_para_s(co.abcd_transicao(np.array([10e9]), t), 50)
        self.assertAlmostEqual(float(co.vswr(s[:, 0, 0])[0]), 1.10, places=6)
        sub = Substrato(3.55, 0.813, 35, 0.0)
        f = np.linspace(1e9, 10e9, 50)
        tr = co.compensacao(t, sub)
        s0 = abcd.abcd_para_s(co.dispositivo_de_teste(f, t, sub, 30), 50)
        s1 = abcd.abcd_para_s(co.dispositivo_de_teste(f, t, sub, 30, comp_trecho=tr), 50)
        self.assertLess(np.abs(s1[:, 0, 0]).max(), 0.5 * np.abs(s0[:, 0, 0]).max())
        dut = co.abcd_linha_ms(f, sintese(70, sub), 20, sub)
        w0 = sintese(50, sub)
        med = abcd.cascata(co.abcd_transicao(f, t, True), co.abcd_linha_ms(f, w0, 5, sub), dut,
                           co.abcd_linha_ms(f, w0, 5, sub), co.abcd_transicao(f, t, False))
        np.testing.assert_allclose(co.de_embutir(f, abcd.abcd_para_s(med, 50), t, sub, 5.0),
                                   abcd.abcd_para_s(dut, 50), atol=1e-9)

    def test_caracterizacao(self):
        from lta.core import caracterizacao as ca, conectores as co
        sub = Substrato(3.55, 0.813)
        w = sintese(50, sub)
        z, ee = (float(x) for x in analise(w, sub))
        f = np.linspace(0.1e9, 10e9, 200)
        t = co.PREDEFINIDOS["Típica, VSWR 1,10 em 10 GHz (Tabela 2.5)"]
        g = ca.duas_linhas(f, ca.sintetico(f, z, ee, 2.0, 30, t), ca.sintetico(f, z, ee, 2.0, 60, t), 30, 60)
        np.testing.assert_allclose((g.imag * ms.C0 / (2 * np.pi * f)) ** 2, ee, rtol=1e-9)   # sem os conectores
        np.testing.assert_allclose(g.real, 2.0, rtol=1e-9)
        ln = ca.linha_de_s(f, ca.sintetico(f, z, ee, 2.0, 30), 30)
        np.testing.assert_allclose(ln.zc.real, z, rtol=1e-9)
        self.assertAlmostEqual(ca.er_de_eeff(ee, w, 0.813), 3.55, places=6)
        ee2, le = ca.dois_ressonadores(*[ms.C0 / (2 * (l + 0.3) * 1e-3 * np.sqrt(ee)) for l in (20.0, 25.0)], 20, 25)
        self.assertAlmostEqual(ee2, ee, places=9)
        self.assertAlmostEqual(le, 0.3, places=9)
        self.assertAlmostEqual(ca.z0_de_vswr(50, 1.5, 1.0), 50 * np.sqrt(1.5))

    def test_robustez(self):
        from lta.core import robustez as rb
        al = Substrato(9.7, 0.635, 10, 2e-4, "Alumina", 4.1e7)
        w = sintese(50, al)
        p = rb.potencia(w, al, np.array([2e9, 10e9, 20e9]), 30.0, 100, 25, 4000)
        np.testing.assert_allclose(p.p_media / np.array([12120, 5170, 3400]), 1, atol=0.3)   # Tabela 2.3
        s = rb.sensibilidades(w, al)
        self.assertAlmostEqual(s["W"][0], -s["h"][0], delta=0.01)       # Z0 depende de W/h
        pc = rb.pior_caso(w, al, 0.0025, 0.025, 0.25)
        self.assertGreater(pc.vswr, 1)
        dw = rb.precisao_w(w, al, pc.dz_rel, 1.0, 0.025, 0.25)
        self.assertAlmostEqual(dw, 0.0025, delta=1e-4)                   # eqs. 2.100-2.101 inversas de 2.97

    def test_dgs(self):
        from lta.core import dgs
        L, C = dgs.lc_de_resposta(8e9, 3.5e9)
        s = abcd.abcd_para_s(dgs.CelulaDGS(L, C).abcd(np.array([3.5e9])), 50)
        self.assertAlmostEqual(abcd.db(s[0, 1, 0]), -3.0103, places=3)   # eq. 6.8: fc é o corte de 3 dB
        el = dgs.passa_baixas(1.3e9, 50, 3, 0.01)
        self.assertAlmostEqual(el[0][1] * 1e9, 3.678, delta=0.2)          # livro: Ld = 3,678 nH, Cs = 2,33 pF
        self.assertAlmostEqual(el[1][1] * 1e12, 2.33, delta=0.06)
        f = np.linspace(1e9, 12e9, 600)
        ref = dgs.ModeloPi(2000, L, C, 5000, 0.05e-12, 0)
        m = dgs.extrair_pi(f, abcd.abcd_para_s(dgs.abcd_pi(f, ref), 50))
        for a, b in ((m.R, 2000), (m.L, L), (m.C, C), (m.Rp, 5000), (m.Cp, 0.05e-12)):
            self.assertAlmostEqual(a / b, 1, places=4)

    def test_layouts_novos(self):
        from lta.geometria import circuitos as ci
        sub = Substrato(3.55, 0.813)
        lays = [ci.linha_dgs("dg", sub, 1.8, 3, 10, 5, 5, 0.5, 2.4, stubs=[(15, 6, 1.0)]),
                ci.slotline("sl", sub, 0.3, 30),
                ci.fixture_conector("fx", sub, 1.8, 30, (0.8, 0.3)),
                ci.balun_marchand("bm", sub, 1.8, 1.2, 0.05, 23, 3.0, chanfro=True),
                ci.wilkinson("wc", sub, 1.8, 1.0, 22, 1.0, 22, chanfro_f0=2e9),
                ci.acoplador_linhas("ac", sub, [(20, 1.5, 0.15)], 1.8, chanfro=True)]
        with tempfile.TemporaryDirectory() as d:
            for lay in lays:
                arqs = exportar(lay, d, list(FORMATOS), dict(f_min_ghz=1, f_max_ghz=5))
                self.assertTrue(all(a.stat().st_size > 0 for a in arqs), lay.nome)
                self._gds_ok((Path(d) / f"{lay.nome}.gds").read_bytes())
        canto = ci.caminho([(0, 0), (10, 0), (10, 10)], 2.0, 60)
        self.assertEqual(len(canto), 3)
        self.assertEqual(len(canto[2]), 5)                                # canto com o chanfro (pentágono)

    def _gds_ok(self, dados):
        i = 0
        while i < len(dados):
            n = struct.unpack(">H", dados[i:i + 2])[0]
            self.assertGreaterEqual(n, 4)
            i += n
        self.assertEqual(i, len(dados))


class SteerSistemas(unittest.TestCase):
    """M. Steer, Microwave and RF Design: A Systems Approach, SciTech, 2010."""

    def test_trl_sintetico(self):
        from lta.core import caracterizacao as ca
        f = np.linspace(1e9, 8e9, 41)
        rng = np.random.default_rng(1)

        def caixa():
            m = abcd.cascata(abcd.series(1j * 2 * np.pi * f * 0.3e-9), abcd.shunt(1j * 2 * np.pi * f * 0.15e-12),
                             abcd.linha(50 + 5 * rng.random(), 2 * np.pi * f / 3e8 * 1.7 * 0.01 - 0.02j))
            return m
        xa, xb = caixa(), caixa()[:, ::-1, ::-1]                           # porta 2 espelhada
        dl, ee = 10.0, 2.7
        lin = abcd.linha(50, 2 * np.pi * f * np.sqrt(ee) / 3e8 * dl * 1e-3)
        dut = abcd.cascata(abcd.linha(30, 2 * np.pi * f * np.sqrt(ee) / 3e8 * 8e-3), abcd.series(5.0 + 0 * f))
        s = lambda m: abcd.abcd_para_s(m, 50)                              # noqa: E731
        sa, sb = s(xa), s(xb)
        for gr, aberto in ((1.0, True), (-1.0, False)):
            g1 = sa[:, 0, 0] + sa[:, 0, 1] * sa[:, 1, 0] * gr / (1 - sa[:, 1, 1] * gr)
            g2 = sb[:, 1, 1] + sb[:, 1, 0] * sb[:, 0, 1] * gr / (1 - sb[:, 0, 0] * gr)
            r = ca.trl(f, s(xa @ xb), s(xa @ lin @ xb), g1, g2, s(xa @ dut @ xb), dl, ee, aberto)
            np.testing.assert_allclose(r.s_dut, s(dut), atol=1e-9)
            np.testing.assert_allclose(r.gama_r, gr, atol=1e-9)

    def test_renormalizacao(self):
        from lta.core import caracterizacao as ca
        f = np.array([1e9, 2e9])
        m = abcd.cascata(abcd.series(10 + 20j + 0 * f), abcd.shunt(0.01j + 0 * f))
        s50 = abcd.abcd_para_s(m, 50)
        np.testing.assert_allclose(ca.renormalizar(s50, 50, 75), abcd.abcd_para_s(m, 75), atol=1e-12)
        np.testing.assert_allclose(ca.renormalizar(ca.renormalizar(s50, 50, [25, 100]), [25, 100], 50), s50,
                                   atol=1e-12)

    def test_estabilidade_phemt(self):
        from lta.core import amplificador as am
        f, s = am.s_de_tabela(am.PHEMT_STEER)
        e = am.estabilidade(s)
        self.assertAlmostEqual(e.k[0], 0.15178, places=4)                   # Tabela 11-6, 0,5 GHz
        self.assertAlmostEqual(e.delta[0], 0.62757, places=4)
        gtu = 10 * np.log10(am.ganhos(s)["GTU_max"])
        np.testing.assert_allclose(gtu[:6], [36.62, 31.07, 24.88, 21.26, 18.73, 17.00], atol=0.01)  # Tab. 11-7
        (cs_, rs), _ = am.circulos(s[:1])
        self.assertAlmostEqual(abs(cs_[0] - (0.92820 + 0.50682j)), 0, places=4)
        self.assertAlmostEqual(rs[0], 0.22434, places=4)
        s8 = s[f == 8e9]
        self.assertTrue(e.incondicional[f == 8e9][0])
        gs, gl = am.conjugado_simultaneo(s8)
        gin, gout = am.gamas_entrada_saida(s8, gs, gl)
        np.testing.assert_allclose(gin, np.conj(gs), atol=1e-12)            # entrada e saída conjugadas
        np.testing.assert_allclose(gout, np.conj(gl), atol=1e-12)
        self.assertAlmostEqual(float(am.gt(s8, gs, gl)[0]), float(am.ganhos(s8)["GMA"][0]), places=9)

    def test_combline_steer(self):
        from lta.core import combline as cl
        p = cl.projeto(1e9, 0.1, 3, ripple_eps=0.1, z_alvo=80)
        # Steer, seção 10.14: C1 = 1,6515 pF, C2 = 2,7171 pF, Z1 = 67,7535 Ω, acoplamento 442,38 Ω, Cb = 1,22306 pF
        for calc, livro in ((p.c[0] * 1e12, 1.6515), (p.c[1] * 1e12, 2.7171), (p.z_stub[0], 67.7535),
                            (p.z_acop[0], 442.38), (p.cb * 1e12, 1.22306)):
            self.assertAlmostEqual(calc / livro, 1, delta=0.01)
        s21 = 20 * np.log10(np.abs(cl.s_combline(np.array([0.5e9, 1e9, 1.5e9, 3e9]), p)[:, 1, 0]))
        self.assertGreater(s21[1], -0.01)
        self.assertTrue(np.all(s21[[0, 2, 3]] < -50))

    def test_rede_pi_ex87(self):
        from lta.core import casamento as cs
        sols = cs.rede_pi(500 + 0j, 50, 1e9, 10)
        self.assertEqual(len(sols), 4)
        self.assertAlmostEqual(sols[0].rv, 500 / 101, places=9)
        for so in sols:
            b1, x2, b3 = so.valores
            self.assertAlmostEqual(1 / abs(b1), 16.58, delta=0.01)             # o livro imprime 16,4
            self.assertAlmostEqual(1 / abs(b3), 50.0, places=6)
            # o livro arredonda RV = 4,95 Ω (exato 4,9505): |X2| = 64,433 ou 34,567 com erro ~0,006
            self.assertTrue(min(abs(abs(x2) - 64.433), abs(abs(x2) - 34.567)) < 0.01)
            zin = cs.zin_abcd(so.abcd(np.array([1e9])), 500)
            self.assertLess(abs(cs.gama(zin, 50))[0], 1e-12)
        for so in cs.rede_t(20 - 35j, 50, 1e9, 4):
            self.assertAlmostEqual(so.rv, 20 * 17, places=9)
            zin = cs.zin_abcd(so.abcd(np.array([1e9])), 20 - 35j)
            self.assertLess(abs(cs.gama(zin, 50))[0], 1e-12)
        with self.assertRaises(ValueError):
            cs.rede_pi(500 + 0j, 50, 1e9, 2.9)                             # Q < Q da rede L (3)

    def test_hibridas_concentradas(self):
        from lta.core import dispositivos as dv
        f = np.array([1e9])
        l35, (se, de) = dv.linha_concentrada(50 / np.sqrt(2), np.pi / 2, 1e9)
        self.assertAlmostEqual(se[1], 50 / np.sqrt(2) / (2 * np.pi * 1e9), places=20)   # L = Z/ω0
        s = dv.branch_line(f, l35, dv.linha_concentrada(50, np.pi / 2, 1e9)[0], 50)[0]
        self.assertAlmostEqual(abs(s[1, 0]) ** 2, 0.5, places=12)
        self.assertLess(abs(s[3, 0]), 1e-12)
        z = 50 * np.sqrt(2)
        l270, (se, de) = dv.linha_concentrada(z, 3 * np.pi / 2, 1e9)
        self.assertEqual((se[0], de[0]), ("C", "L"))                       # 270° → -90° passa-altas
        s = dv.rat_race(f, dv.linha_concentrada(z, np.pi / 2, 1e9)[0], l270, 50)[0]
        self.assertAlmostEqual(abs(s[1, 0]) ** 2 + abs(s[2, 0]) ** 2, 1, places=12)

    def test_acoplador_capacitores_pontas(self):
        from lta.core import linhas_acopladas as la
        f = np.array([2e9])
        conv = la.projeto_acoplador(la.TIPOS[0], 50, 10, 2e9, RO3010)
        pont = la.projeto_acoplador(la.TIPOS[0], 50, 10, 2e9, RO3010, pontas=True)
        d_conv = la.metricas(la.acoplador(f, conv.secoes, 50))[1][0]
        c, d, _, rl = (x[0] for x in la.metricas(la.acoplador(f, pont.secoes, 50)))
        self.assertLess(d_conv, 20)
        self.assertGreater(d, 60)
        self.assertGreater(rl, 60)
        self.assertAlmostEqual(c, 10, delta=0.05)
        self.assertGreater(pont.secoes[0].c_pontas, 0)

    def test_atenuadores(self):
        from lta.core import passivos as pv
        t = pv.projeto("T", 50, 50, 10)
        np.testing.assert_allclose(t.r, (25.975, 25.975, 35.136), atol=1e-3)
        pi = pv.projeto("Pi", 50, 50, 10)
        np.testing.assert_allclose(pi.r, (96.248, 96.248, 71.151), atol=1e-3)
        self.assertAlmostEqual(10 * np.log10(pv.k_minimo(50, 75)), 5.719, places=3)   # eq. 7.14
        for top in pv.TOPOLOGIAS:
            a = pv.projeto(top, 50, 75, 8)
            s = a.s(np.array([1e9]))[0]
            self.assertLess(abs(s[0, 0]), 1e-12)
            self.assertLess(abs(s[1, 1]), 1e-12)
            self.assertAlmostEqual(-20 * np.log10(abs(s[1, 0])), 8, places=9)
            self.assertAlmostEqual(sum(a.potencias(1.0)), 1, places=12)
        with self.assertRaises(ValueError):
            pv.projeto("T", 50, 75, 5)

    def test_marchand_miniaturizado(self):
        from lta.core import balun
        b = balun.projeto_miniaturizado(1e9, 50, 200, 45)
        m = balun.metricas(balun.s_marchand(np.array([1e9]), *b.secoes, 50, 200, caps=b.caps))
        self.assertLess(m[0][0], -60)
        self.assertLess(abs(m[3][0]), 0.01)
        self.assertAlmostEqual(abs(m[4][0]), 180, delta=0.1)
        conv = balun._secoes_ideais(0.5, 50, 1e9)[0].comp
        self.assertAlmostEqual(b.secoes[0].comp / conv, 0.5, places=6)


if __name__ == "__main__":
    unittest.main()
