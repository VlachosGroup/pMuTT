# -*- coding: utf-8 -*-
"""
pMuTT.test_pMuTT_model_statmech
Tests for pMuTT module
"""
import unittest
from pMuTT import constants as c
from pMuTT.eos import IdealGasEOS, vanDerWaalsEOS


class TestIdealGas(unittest.TestCase):
    def setUp(self):
        self.ideal_gas = IdealGasEOS()

    def test_get_V(self):
        self.assertAlmostEqual(
            self.ideal_gas.get_V(P=c.P0('bar'), T=c.T0('K'), n=1.),
            c.V0('m3'))

    def test_get_T(self):
        self.assertAlmostEqual(
            self.ideal_gas.get_T(P=c.P0('bar'), V=c.V0('m3'), n=1.),
            c.T0('K'))

    def test_get_P(self):
        self.assertAlmostEqual(
            self.ideal_gas.get_P(T=c.T0('K'), V=c.V0('m3'), n=1.),
            c.P0('bar'))

    def test_get_n(self):
        self.assertAlmostEqual(
            self.ideal_gas.get_n(T=c.T0('K'), V=c.V0('m3'), P=c.P0('bar')),
            1.)


class TestvanDerWaalsEOS(unittest.TestCase):
    def setUp(self):
        self.van_der_waals = vanDerWaalsEOS(a=0.547, b=30.52e-6)
        self.van_der_waals_critical = vanDerWaalsEOS.from_critical(Tc=638.73248, 
                                                                   Pc=217.49762)
        self.T = 500.  # K
        self.V = 0.01  # m3
        self.n = 1.  # mol
        self.P = 4.1037990739174  # bar

    def test_get_V(self):
        self.assertAlmostEqual(
            self.van_der_waals.get_V(P=self.P, T=self.T, n=self.n),
            0.010028206752280285)
        self.assertAlmostEqual(
            self.van_der_waals_critical.get_V(P=self.P, T=self.T, n=self.n),
            0.01002819363297158)

    def test_get_T(self):
        self.assertAlmostEqual(
            self.van_der_waals.get_T(P=self.P, V=self.V, n=self.n),
            498.6261807103575)
        self.assertAlmostEqual(
            self.van_der_waals_critical.get_T(P=self.P, V=self.V, n=self.n),
            498.62682174170465)

    def test_get_P(self):
        self.assertAlmostEqual(
            self.van_der_waals.get_P(V=self.V, T=self.T, n=self.n),
            4.115256607566292)
        self.assertAlmostEqual(
            self.van_der_waals_critical.get_P(V=self.V, T=self.T, n=self.n),
            4.11525126335907)

    def test_get_n(self):
        self.assertAlmostEqual(
            self.van_der_waals.get_n(P=self.P, T=self.T, V=self.V),
            0.9971872586019559)
        self.assertAlmostEqual(
            self.van_der_waals_critical.get_n(P=self.P, T=self.T, V=self.V),
            0.997188563164668)


if __name__ == '__main__':
    unittest.main()
