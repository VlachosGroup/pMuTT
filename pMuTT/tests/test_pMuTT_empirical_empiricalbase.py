# -*- coding: utf-8 -*-
"""
pMuTT.test_pMuTT
Tests for pMuTT module
Created on Fri Jul 7 12:31:00 2018
"""
import unittest
import numpy as np
from ase.build import molecule
from pMuTT.statmech import StatMech, presets
from pMuTT.empirical import EmpiricalBase


class TestEmpiricalBase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.EmpiricalBase = EmpiricalBase(
            name='H2O',
            elements={'H': 2, 'O': 1},
            phase='g',
            potentialenergy=-14.2209,
            atoms=molecule('H2O'),
            symmetrynumber=2,
            spin=0,
            vib_wavenumbers=np.array([3825.434, 3710.2642, 1582.432]),
            **presets['idealgas'])


if __name__ == '__main__':
    unittest.main()
