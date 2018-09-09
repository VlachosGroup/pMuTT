# -*- coding: utf-8 -*-
"""
PyMuTT.test_PyMuTT_model_statmech_vib
Tests for PyMuTT module
"""
import unittest
import numpy as np
from PyMuTT.models.statmech import vib

class TestHarmonicVib(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)

        self.vib_H2 = vib.HarmonicVib(vib_wavenumbers=[4306.1793])
        self.vib_H2O = vib.HarmonicVib(
            vib_wavenumbers=[3825.434, 3710.2642, 1582.432])
        self.T = 300 # K

    def test_get_q(self):
        self.assertAlmostEqual(self.vib_H2.get_q(T=self.T), 3.27680884e-05)
        self.assertAlmostEqual(self.vib_H2O.get_q(T=self.T), 3.1464834E-10)

    def test_get_CvoR(self):
        self.assertAlmostEqual(self.vib_H2.get_CvoR(T=self.T), 4.52088E-07)
        self.assertAlmostEqual(self.vib_H2O.get_CvoR(T=self.T), 0.02917545)

    def test_get_CpoR(self):
        self.assertAlmostEqual(self.vib_H2.get_CpoR(T=self.T), 4.52088E-07)
        self.assertAlmostEqual(self.vib_H2O.get_CpoR(T=self.T), 0.02917545)

    def test_get_ZPE(self):
        self.assertAlmostEqual(self.vib_H2.get_ZPE(), 0.26694909102484027)
        self.assertAlmostEqual(self.vib_H2O.get_ZPE(), 0.5652520248602153)        

    def test_get_UoRT(self):
        self.assertAlmostEqual(self.vib_H2.get_UoRT(T=self.T), 10.32605545)
        self.assertAlmostEqual(self.vib_H2O.get_UoRT(T=self.T), 21.8687737)

    def test_get_HoRT(self):
        self.assertAlmostEqual(self.vib_H2.get_HoRT(T=self.T), 10.32605545)
        self.assertAlmostEqual(self.vib_H2O.get_HoRT(T=self.T), 21.8687737)

    def test_get_SoR(self):
        self.assertAlmostEqual(self.vib_H2.get_SoR(T=self.T), 2.32489026469e-08)
        self.assertAlmostEqual(self.vib_H2O.get_SoR(T=self.T), 0.00434769)

    def test_get_AoRT(self):
        self.assertAlmostEqual(self.vib_H2.get_AoRT(T=self.T), 1.032605543E+01)
        self.assertAlmostEqual(self.vib_H2O.get_AoRT(T=self.T), 2.186442601E+01)
        
    def test_get_GoRT(self):
        self.assertAlmostEqual(self.vib_H2.get_GoRT(T=self.T), 1.032605543E+01)
        self.assertAlmostEqual(self.vib_H2O.get_GoRT(T=self.T), 2.186442601E+01)

if __name__ == '__main__':
    unittest.main()