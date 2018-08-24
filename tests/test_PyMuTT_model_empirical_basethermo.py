# -*- coding: utf-8 -*-
"""
PyMuTT.test_PyMuTT
Tests for PyMuTT module
Created on Fri Jul 7 12:31:00 2018
"""
import unittest
import numpy as np
from ase.build import molecule
from ase.thermochemistry import IdealGasThermo
from PyMuTT.models.empirical import BaseThermo

class TestBaseThermo(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.BaseThermo = BaseThermo(
            name = 'H2O',
            elements = {'H': 2, 'O': 1},
            phase = 'g',
            thermo_model = IdealGasThermo,
            potentialenergy = -14.2209,
            geometry = 'nonlinear',
            atoms = molecule('H2O'),
            symmetrynumber = 2,
            spin = 0,
            vib_energies = np.array([0.47462, 0.46033, 0.19633])
            )

if __name__ == '__main__':
    unittest.main()