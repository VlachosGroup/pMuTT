# -*- coding: utf-8 -*-
"""
pMuTT.test_pMuTT_model_statmech
Tests for pMuTT module
"""
import unittest
import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.thermochemistry import IdealGasThermo

from pMuTT import get_molecular_weight, _pass_expected_arguments
from pMuTT import constants as c
from pMuTT.statmech import trans, rot, elec, vib, StatMech

class TestStatMech(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        # Testing Ideal Gas Model
        CO2 = molecule('CO2')
        CO2_pMuTT_parameters = {
            'name': 'CO2',

            'trans_model': trans.IdealTrans,
            'n_degrees': 3,
            'molecular_weight': get_molecular_weight('CO2'),

            'rot_model': rot.RigidRotor,
            'rot_temperatures': rot.get_rot_temperatures_from_atoms(CO2, 
                geometry='linear'),
            'geometry': 'linear',
            'symmetrynumber': 2,

            'elec_model': elec.IdealElec,
            'potentialenergy': -22.994202,
            'spin': 0.,

            'vib_model': vib.HarmonicVib,
            'vib_wavenumbers': [3360., 954., 954., 1890.],
        }
        CO2_ase_parameters = {
            'atoms': CO2,
            'potentialenergy': -22.994202,
            'vib_energies': [c.wavenumber_to_energy(x) \
                for x in CO2_pMuTT_parameters['vib_wavenumbers']],
            'geometry':'linear',
            'symmetrynumber': 2,
            'spin': 0.
        }
        self.CO2_pMuTT = StatMech(**CO2_pMuTT_parameters)
        self.CO2_ASE = IdealGasThermo(**CO2_ase_parameters)
        
        self.T0 = c.T0('K') # K
        self.P0 = c.P0('Pa')
        self.V0 = c.V0('m3')

    def test_get_q(self):
        self.assertAlmostEqual(
            self.CO2_pMuTT.get_q(T=self.T0, ignore_q_elec=True, V=self.V0), 
                6.083051624373337e+25)

    def test_get_CvoR(self):
        self.assertAlmostEqual(
            self.CO2_pMuTT.get_CvoR(T=self.T0, V=self.V0), 
                2.9422622359004853)

    def test_get_CpoR(self):
        self.assertAlmostEqual(
            self.CO2_pMuTT.get_CpoR(T=self.T0, V=self.V0), 
                3.9422622359004853)

    def test_get_UoRT(self):
        self.assertAlmostEqual(
            self.CO2_pMuTT.get_UoRT(T=self.T0, V=self.V0), 
                -875.1095022368354)

    def test_get_HoRT(self):
        expected_HoRT_CO2 = \
            self.CO2_ASE.get_enthalpy(temperature=self.T0, verbose=False) \
            /c.R('eV/K')/self.T0
        calc_HoRT_CO2 = self.CO2_pMuTT.get_HoRT(T=self.T0)
        self.assertTrue(np.isclose(expected_HoRT_CO2, calc_HoRT_CO2))

    def test_get_SoR(self):
        expected_SoR_CO2 = \
            self.CO2_ASE.get_entropy(temperature=self.T0, pressure=self.P0,
                verbose=False)/c.R('eV/K')
        calc_SoR_CO2 = self.CO2_pMuTT.get_SoR(T=self.T0, V=self.V0)
        self.assertTrue(np.isclose(expected_SoR_CO2, calc_SoR_CO2))
        
    def test_get_AoRT(self):
        self.assertAlmostEqual(
            self.CO2_pMuTT.get_AoRT(T=self.T0, V=self.V0), 
                -900.5899966269182)
        
    def test_get_GoRT(self):
        expected_GoRT_CO2 = \
            self.CO2_ASE.get_gibbs_energy(temperature=self.T0, pressure=self.P0,
               verbose=False)/c.R('eV/K')/self.T0
        calc_GoRT_CO2 = self.CO2_pMuTT.get_GoRT(T=self.T0, V=self.V0)
        self.assertTrue(np.isclose(expected_GoRT_CO2, calc_GoRT_CO2))

if __name__ == '__main__':
    unittest.main()