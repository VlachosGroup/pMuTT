# -*- coding: utf-8 -*-
"""
pMuTT.test_pMuTT_model_statmech_vib
Tests for pMuTT module
"""
import unittest
import numpy as np
from pMuTT import constants as c
from pMuTT.models.statmech import vib

class TestHarmonicVib(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)

        self.vib_H2 = vib.HarmonicVib(vib_wavenumbers=[4306.1793])
        self.vib_H2O = vib.HarmonicVib(
                vib_wavenumbers=[3825.434, 3710.2642, 1582.432])
        self.vib_H2O_dict = {
            'class': "<class 'pMuTT.models.statmech.vib.HarmonicVib'>",
            'vib_wavenumbers': [3825.434, 3710.2642, 1582.432],
        }
        self.T = 300. #K

    def test_get_q(self):
        self.assertAlmostEqual(self.vib_H2.get_q(T=self.T), 
                               3.27680884e-05)
        self.assertAlmostEqual(self.vib_H2O.get_q(T=self.T), 
                               3.1464834E-10)

    def test_get_CvoR(self):
        self.assertAlmostEqual(self.vib_H2.get_CvoR(T=self.T), 4.52088E-07)
        self.assertAlmostEqual(self.vib_H2O.get_CvoR(T=self.T), 0.02917545)

    def test_get_CpoR(self):
        self.assertAlmostEqual(self.vib_H2.get_CpoR(T=self.T), 4.52088E-07)
        self.assertAlmostEqual(self.vib_H2O.get_CpoR(T=self.T), 0.02917545)

    def test_get_ZPE(self):
        self.assertAlmostEqual(self.vib_H2.get_ZPE(), 
                               0.26694909102484027)
        self.assertAlmostEqual(self.vib_H2O.get_ZPE(), 
                               0.5652520248602153)        

    def test_get_UoRT(self):
        self.assertAlmostEqual(self.vib_H2.get_UoRT(T=self.T), 10.32605545)
        self.assertAlmostEqual(self.vib_H2O.get_UoRT(T=self.T), 21.8687737)

    def test_get_HoRT(self):
        self.assertAlmostEqual(self.vib_H2.get_HoRT(T=self.T), 10.32605545)
        self.assertAlmostEqual(self.vib_H2O.get_HoRT(T=self.T), 21.8687737)

    def test_get_SoR(self):
        self.assertAlmostEqual(self.vib_H2.get_SoR(T=self.T), 
                               2.32489026469e-08)
        self.assertAlmostEqual(self.vib_H2O.get_SoR(T=self.T), 0.00434769)

    def test_get_AoRT(self):
        self.assertAlmostEqual(self.vib_H2.get_AoRT(T=self.T), 
                               1.032605543E+01)
        self.assertAlmostEqual(self.vib_H2O.get_AoRT(T=self.T), 
                               2.186442601E+01)
        
    def test_get_GoRT(self):
        self.assertAlmostEqual(self.vib_H2.get_GoRT(T=self.T), 
                               1.032605543E+01)
        self.assertAlmostEqual(self.vib_H2O.get_GoRT(T=self.T), 
                               2.186442601E+01)

    def test_to_dict(self):
        self.assertEqual(self.vib_H2O.to_dict(), self.vib_H2O_dict)

    def test_from_dict(self):
        self.assertEqual(vib.HarmonicVib.from_dict(self.vib_H2O_dict), 
                         self.vib_H2O)

class TestQRRHOVib(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)

        self.vib_H2 = vib.QRRHOVib(vib_wavenumbers=[4306.1793])
        self.vib_H2O = vib.QRRHOVib(
                vib_wavenumbers=[3825.434, 3710.2642, 1582.432])
        self.vib_H2O_dict = {
            'class': "<class 'pMuTT.models.statmech.vib.QRRHOVib'>",
            'vib_wavenumbers': [3825.434, 3710.2642, 1582.432],
            'alpha': 4,
            'Bav': 1.e-44,
            'v0': 100,
        }
        self.T = 300. # K

    def test_get_scaled_wavenumber(self):
        self.assertAlmostEqual(self.vib_H2O._get_scaled_wavenumber(1582.432),
                               0.999984052)
    def test_get_scaled_inertia(self):
        self.assertAlmostEqual(self.vib_H2O._get_scaled_inertia(1582.432),
                               1.76893896253E-49)

    def test_get_CvoR(self):
        self.assertAlmostEqual(self.vib_H2.get_CvoR(T=self.T), 
                               6.0337598762E-07)
        self.assertAlmostEqual(self.vib_H2O.get_CvoR(T=self.T), 
                               2.918349716E-02)

    def test_get_CpoR(self):
        self.assertAlmostEqual(self.vib_H2.get_CpoR(T=self.T), 
                               6.0337598762E-07)
        self.assertAlmostEqual(self.vib_H2O.get_CpoR(T=self.T), 
                               2.918349716E-02)

    def test_get_UoRT_RHHO(self):
        self.assertAlmostEqual(
            self.vib_H2O._get_UoRT_RRHO(T=self.T, vib_temperature=2276.767335),
            3.798453353623928)

    def test_get_UoRT(self):
        self.assertAlmostEqual(self.vib_H2.get_UoRT(T=self.T), 
                               10.3260525951174)
        self.assertAlmostEqual(self.vib_H2O.get_UoRT(T=self.T), 
                               21.868712644411)

    def test_get_HoRT(self):
        self.assertAlmostEqual(self.vib_H2.get_UoRT(T=self.T), 
                               10.3260525951174)
        self.assertAlmostEqual(self.vib_H2O.get_HoRT(T=self.T), 
                               21.868712644411)

    def test_get_SoR_H(self):
        self.assertAlmostEqual(
            self.vib_H2O._get_SoR_H(T=self.T, vib_temperature=2276.767335),
            0.0043471298500)

    def test_get_SoR_RRHO(self):
        self.assertAlmostEqual(
                self.vib_H2O._get_SoR_RRHO(T=self.T, 
                                           vib_inertia=1.768938963E-49),
                5.899139738E-02)

    def test_get_SoR(self):
        self.assertAlmostEqual(self.vib_H2.get_SoR(T=self.T), 
                               1.6315868671e-06)
        self.assertAlmostEqual(self.vib_H2O.get_SoR(T=self.T), 
                               0.00444131527822)

    def test_get_AoRT(self):
        self.assertAlmostEqual(self.vib_H2.get_AoRT(T=self.T), 
                               10.3260509635305)
        self.assertAlmostEqual(self.vib_H2O.get_AoRT(T=self.T), 
                               21.864271329132)
        
    def test_get_GoRT(self):
        self.assertAlmostEqual(self.vib_H2.get_GoRT(T=self.T), 
                               10.3260509635305)
        self.assertAlmostEqual(self.vib_H2O.get_GoRT(T=self.T), 
                               21.864271329132)

    def test_to_dict(self):
        self.assertEqual(self.vib_H2O.to_dict(), self.vib_H2O_dict)

    def test_from_dict(self):
        self.assertEqual(vib.QRRHOVib.from_dict(self.vib_H2O_dict), 
                self.vib_H2O)

class TestEinsteinVib(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)

        self.vib_Ag = vib.EinsteinVib(einstein_temperature=173.28913505677, 
                                      interaction_energy=0.5)
        self.vib_Ag_dict = {
            'class': "<class 'pMuTT.models.statmech.vib.EinsteinVib'>",
            'einstein_temperature': 173.28913505677,
            'interaction_energy': 0.5
        }
        self.T = 300. #K

    def test_get_q(self):
        self.assertAlmostEqual(self.vib_Ag.get_q(T=self.T), 6.8029252699e-09)

    def test_get_CvoR(self):
        self.assertAlmostEqual(self.vib_Ag.get_CvoR(T=self.T), 2.9179591404057)

    def test_get_CpoR(self):
        self.assertAlmostEqual(self.vib_Ag.get_CpoR(T=self.T), 2.9179591404057)

    def test_get_ZPE(self):
        self.assertAlmostEqual(self.vib_Ag.get_ZPE(), 0.5223993457128)

    def test_get_UoRT(self):
        self.assertAlmostEqual(self.vib_Ag.get_UoRT(T=self.T), 22.4238242106867)

    def test_get_HoRT(self):
        self.assertAlmostEqual(self.vib_Ag.get_HoRT(T=self.T), 22.4238242106867)

    def test_get_SoR(self):
        self.assertAlmostEqual(self.vib_Ag.get_SoR(T=self.T), 4.6878251875495)

    def test_get_AoRT(self):
        self.assertAlmostEqual(self.vib_Ag.get_AoRT(T=self.T), 17.7359990231372)
        
    def test_get_GoRT(self):
        self.assertAlmostEqual(self.vib_Ag.get_GoRT(T=self.T), 17.7359990231372)

    def test_to_dict(self):
        self.assertEqual(self.vib_Ag.to_dict(), self.vib_Ag_dict)

    def test_from_dict(self):
        self.assertEqual(vib.EinsteinVib.from_dict(self.vib_Ag_dict), 
                self.vib_Ag)

class TestHelperFunc(unittest.TestCase):
    def test_debye_to_einstein(self):
        self.assertAlmostEqual(vib.debye_to_einstein(215.), 173.28913505677)

    def test_einstein_to_debye(self):
        self.assertAlmostEqual(vib.einstein_to_debye(173.28913505677), 215.)

if __name__ == '__main__':
    unittest.main()