# -*- coding: utf-8 -*-
"""
pmutt.test_pmutt
Tests for pmutt module
Created on Fri Jul 7 12:31:00 2018
"""
import unittest
import numpy as np
from ase.build import molecule
from pmutt.statmech import StatMech, presets
from pmutt.empirical import EmpiricalBase, GasPressureAdj


class TestEmpiricalBase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.EmpiricalBase = EmpiricalBase(
            name='H2O',
            elements={
                'H': 2,
                'O': 1
            },
            phase='g',
            potentialenergy=-14.2209,
            atoms=molecule('H2O'),
            symmetrynumber=2,
            spin=0,
            vib_wavenumbers=np.array([3825.434, 3710.2642, 1582.432]),
            **presets['idealgas'])


class TestGasPressureAdj(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.obj = GasPressureAdj()

    def test_get_SoR(self):
        P = np.logspace(-10, 10)
        exp_SoR = -np.log(P)
        np.testing.assert_array_almost_equal(exp_SoR, self.obj.get_SoR(P=P))

    def test_get_GoRT(self):
        P = np.logspace(-10, 10)
        exp_GoRT = np.log(P)
        np.testing.assert_array_almost_equal(exp_GoRT, self.obj.get_GoRT(P=P))


if __name__ == '__main__':
    unittest.main()
