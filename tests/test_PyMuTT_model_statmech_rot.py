# -*- coding: utf-8 -*-
"""
PyMuTT.test_PyMuTT_model_statmech_rot
Tests for PyMuTT module
"""
import unittest
import numpy as np
from ase import Atoms
from ase.build import molecule
from PyMuTT.models.statmech import rot

class TestRigidRotor(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.rot_He = rot.RigidRotor(symmetrynumber=1, geometry='monatomic',
                                   rot_temperatures=[0.])
        self.rot_CO2 = rot.RigidRotor(symmetrynumber=2, geometry='linear',
                                    rot_temperatures=[0.561])
        self.rot_H2O = rot.RigidRotor(symmetrynumber=2, geometry='nonlinear',
                                    rot_temperatures=[40.1, 20.9, 13.4])
        self.T = 300 # K

    def test_get_q(self):
        self.assertAlmostEqual(self.rot_He.get_q(T=self.T), 0.)
        self.assertAlmostEqual(self.rot_CO2.get_q(T=self.T), 267.3796791)
        self.assertAlmostEqual(self.rot_H2O.get_q(T=self.T), 43.45393338)

    def test_get_CvoR(self):
        self.assertAlmostEqual(self.rot_He.get_CvoR(), 0.)
        self.assertAlmostEqual(self.rot_CO2.get_CvoR(), 1.)
        self.assertAlmostEqual(self.rot_H2O.get_CvoR(), 1.5)

    def test_get_CpoR(self):
        self.assertAlmostEqual(self.rot_He.get_CpoR(), 0.)
        self.assertAlmostEqual(self.rot_CO2.get_CpoR(), 1.)
        self.assertAlmostEqual(self.rot_H2O.get_CpoR(), 1.5)

    def test_get_UoRT(self):
        self.assertAlmostEqual(self.rot_He.get_UoRT(), 0.)
        self.assertAlmostEqual(self.rot_CO2.get_UoRT(), 1.)
        self.assertAlmostEqual(self.rot_H2O.get_UoRT(), 1.5)

    def test_get_HoRT(self):
        self.assertAlmostEqual(self.rot_He.get_HoRT(), 0.)
        self.assertAlmostEqual(self.rot_CO2.get_HoRT(), 1.)
        self.assertAlmostEqual(self.rot_H2O.get_HoRT(), 1.5)

    def test_get_SoR(self):
        self.assertAlmostEqual(self.rot_He.get_SoR(T=self.T), 0.)
        self.assertAlmostEqual(self.rot_CO2.get_SoR(T=self.T), 6.588669668)
        self.assertAlmostEqual(self.rot_H2O.get_SoR(T=self.T), 5.271701374)
        
    def test_get_AoRT(self):
        self.assertAlmostEqual(self.rot_He.get_AoRT(T=self.T), 0.)
        self.assertAlmostEqual(self.rot_CO2.get_AoRT(T=self.T), -5.588669668)
        self.assertAlmostEqual(self.rot_H2O.get_AoRT(T=self.T), -3.771701374)
        
    def test_get_GoRT(self):
        self.assertAlmostEqual(self.rot_He.get_GoRT(T=self.T), 0.)
        self.assertAlmostEqual(self.rot_CO2.get_GoRT(T=self.T), -5.588669668)
        self.assertAlmostEqual(self.rot_H2O.get_GoRT(T=self.T), -3.771701374)

class TestRotFunc(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.He = Atoms('He')
        self.CO2 = molecule('CO2')
        self.H2O = molecule('H2O')
        self.CH4 = molecule('CH4')
        self.NH3 = molecule('NH3')
        self.CO = molecule('CO')

    def test_get_rot_temperatures_from_atoms(self):
        self.assertListEqual(
            rot.get_rot_temperatures_from_atoms(self.He, geometry='monatomic'), 
                                                [0.])

        rot_Ts_CO2 = rot.get_rot_temperatures_from_atoms(self.CO2, 
            geometry='linear')
        self.assertTrue(len(rot_Ts_CO2), 1)
        self.assertTrue(np.isclose(rot_Ts_CO2[0], 0.545566039279433))

        exp_rot_Ts_H2O = [38.0937696114643, 20.650669844110, 13.3912565453767]
        rot_Ts_H2O = rot.get_rot_temperatures_from_atoms(self.H2O,
                                                         geometry='nonlinear')
        self.assertTrue(len(rot_Ts_H2O), 3)
        for rot_T in rot_Ts_H2O:
            self.assertTrue(any(np.isclose(rot_T, exp_rot_T)
                            for exp_rot_T in exp_rot_Ts_H2O))
        
    def test_get_geometry_from_atoms(self):
        self.assertEqual(rot.get_geometry_from_atoms(self.He), 'monatomic')
        self.assertEqual(rot.get_geometry_from_atoms(self.CO2), 'linear')
        self.assertEqual(rot.get_geometry_from_atoms(self.H2O), 'nonlinear')
        self.assertEqual(rot.get_geometry_from_atoms(self.CH4), 'nonlinear')
        self.assertEqual(rot.get_geometry_from_atoms(self.NH3), 'nonlinear')
        self.assertEqual(rot.get_geometry_from_atoms(self.CO), 'linear')

if __name__ == '__main__':
    unittest.main()