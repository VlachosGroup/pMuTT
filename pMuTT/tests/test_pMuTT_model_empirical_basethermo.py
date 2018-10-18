# -*- coding: utf-8 -*-
"""
pMuTT.test_pMuTT
Tests for pMuTT module
Created on Fri Jul 7 12:31:00 2018
"""
import unittest
import numpy as np
from ase.build import molecule
from pMuTT.models.statmech import StatMech, presets
from pMuTT.models.empirical import BaseThermo

class TestBaseThermo(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.BaseThermo = BaseThermo(
            name = 'H2O',
            elements = {'H': 2, 'O': 1},
            phase = 'g',
            statmech_model = StatMech,
            potentialenergy = -14.2209,
            atoms = molecule('H2O'),
            symmetrynumber = 2,
            spin = 0,
            vib_wavenumbers = np.array([3825.434, 3710.2642, 1582.432]),
            **presets['idealgas'])

if __name__ == '__main__':
    unittest.main()