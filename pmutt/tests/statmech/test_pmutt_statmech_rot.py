# -*- coding: utf-8 -*-
"""
pmutt.test_pmutt_model_statmech_rot
Tests for pmutt module
"""
import unittest
import numpy as np
from ase import Atoms
from ase.build import molecule
from pmutt.statmech import rot


class TestRigidRotor(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.rot_He = rot.RigidRotor(symmetrynumber=1,
                                     geometry='monatomic',
                                     rot_temperatures=[0.])
        self.rot_CO2 = rot.RigidRotor(symmetrynumber=2,
                                      geometry='linear',
                                      rot_temperatures=[0.561])
        self.rot_H2O = rot.RigidRotor(symmetrynumber=2,
                                      geometry='nonlinear',
                                      rot_temperatures=[40.1, 20.9, 13.4])
        self.rot_H2O_pymatgen = rot.RigidRotor(symmetrynumber='pymatgen',
                                               geometry='nonlinear',
                                               rot_temperatures=[40.1, 20.9,
                                                                 13.4])
        self.T = 300  # K

        self.rot_CO2_dict = {
            'class': "<class 'pmutt.statmech.rot.RigidRotor'>",
            'symmetrynumber': 2,
            'geometry': 'linear',
            'rot_temperatures': [0.561],
        }

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

    def test_get_FoRT(self):
        self.assertAlmostEqual(self.rot_He.get_FoRT(T=self.T), 0.)
        self.assertAlmostEqual(self.rot_CO2.get_FoRT(T=self.T), -5.588669668)
        self.assertAlmostEqual(self.rot_H2O.get_FoRT(T=self.T), -3.771701374)

    def test_get_GoRT(self):
        self.assertAlmostEqual(self.rot_He.get_GoRT(T=self.T), 0.)
        self.assertAlmostEqual(self.rot_CO2.get_GoRT(T=self.T), -5.588669668)
        self.assertAlmostEqual(self.rot_H2O.get_GoRT(T=self.T), -3.771701374)
        self.assertAlmostEqual(self.rot_H2O_pymatgen.get_GoRT(T=self.T),
                               self.rot_H2O.get_GoRT(T=self.T))

    def test_to_dict(self):
        self.assertEqual(self.rot_CO2.to_dict(), self.rot_CO2_dict)

    def test_from_dict(self):
        self.assertEqual(rot.RigidRotor.from_dict(self.rot_CO2_dict),
                         self.rot_CO2)


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
        np.testing.assert_almost_equal(rot_Ts_CO2[0], 0.5456089507450469)

        exp_rot_Ts_H2O = [38.096765874734196,
		                  20.652294121365664,
		                  13.392309833894341]
        rot_Ts_H2O = rot.get_rot_temperatures_from_atoms(self.H2O,
                                                         geometry='nonlinear')
        self.assertTrue(len(rot_Ts_H2O), 3)
        np.testing.assert_array_almost_equal(rot_Ts_H2O, exp_rot_Ts_H2O)

    def test_get_geometry_from_atoms(self):
        self.assertEqual(rot.get_geometry_from_atoms(self.He), 'monatomic')
        self.assertEqual(rot.get_geometry_from_atoms(self.CO2), 'linear')
        self.assertEqual(rot.get_geometry_from_atoms(self.H2O), 'nonlinear')
        self.assertEqual(rot.get_geometry_from_atoms(self.CH4), 'nonlinear')
        self.assertEqual(rot.get_geometry_from_atoms(self.NH3), 'nonlinear')
        self.assertEqual(rot.get_geometry_from_atoms(self.CO), 'linear')


if __name__ == '__main__':
    unittest.main()
