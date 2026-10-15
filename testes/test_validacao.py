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


class Exportacao(unittest.TestCase):
    def test_todos_formatos(self):
        p = cm.perfil("Senoidal", 64.0, 0.01, 10)
        w, _ = cm.larguras(p, RO3010)
        lays = [construtores.linha_largura_modulada(p.z * 1e3, w, RO3010, 5.0),
                construtores.linha_furos_terra(18.9, 4.725, 7, 50, RO3010),
                construtores.linha_carga_capacitiva(
                    cc.geometria(cc.projeto_por_fase(35, 1.4e9, 47.5, 22e9), ARLON, 10, 50), 2, ARLON)]
        with tempfile.TemporaryDirectory() as d:
            for lay in lays:
                arqs = exportar(lay, d, list(FORMATOS), dict(f_min_ghz=1, f_max_ghz=10))
                for a in arqs:
                    self.assertGreater(a.stat().st_size, 0, a)
                txt = (Path(d) / f"{lay.nome}.kicad_pcb").read_text()
                self.assertEqual(txt.count("("), txt.count(")"))
                compile((Path(d) / f"{lay.nome}_hfss.py").read_text(), "hfss", "exec")
                self._gds_valido((Path(d) / f"{lay.nome}.gds").read_bytes())
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
