# -*- coding: utf-8 -*-
"""
pmutt.test_pmutt_model_statmech
Tests for pmutt module
"""
import unittest
import numpy as np
from ase.build import molecule
from ase.thermochemistry import IdealGasThermo
from pmutt import get_molecular_weight
from pmutt import constants as c
from pmutt.statmech import trans, rot, elec, vib, StatMech


class TestStatMech(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        # Testing Ideal Gas Model
        CO2 = molecule('CO2')
        CO2_pmutt_parameters = {
            'name':
            'CO2',
            'elements': {
                'C': 1,
                'O': 2
            },
            'trans_model':
            trans.FreeTrans,
            'n_degrees':
            3,
            'molecular_weight':
            get_molecular_weight('CO2'),
            'rot_model':
            rot.RigidRotor,
            'rot_temperatures':
            rot.get_rot_temperatures_from_atoms(CO2, geometry='linear'),
            'geometry':
            'linear',
            'symmetrynumber':
            2,
            'elec_model':
            elec.GroundStateElec,
            'potentialenergy':
            -22.994202,
            'spin':
            0.,
            'vib_model':
            vib.HarmonicVib,
            'vib_wavenumbers': [3360., 954., 954., 1890.],
        }
        CO2_ase_parameters = {
            'atoms':
            CO2,
            'potentialenergy':
            -22.994202,
            'vib_energies': [
                c.wavenumber_to_energy(x) *
                c.convert_unit(initial='J', final='eV')
                for x in CO2_pmutt_parameters['vib_wavenumbers']
            ],
            'geometry':
            'linear',
            'symmetrynumber':
            2,
            'spin':
            0.
        }
        self.CO2_pmutt = StatMech(**CO2_pmutt_parameters)
        self.CO2_ASE = IdealGasThermo(**CO2_ase_parameters)

        self.T0 = c.T0('K')  # K
        self.P0 = c.P0('Pa')
        self.V0 = c.V0('m3')
        self.mw = get_molecular_weight({'C': 1, 'O': 2})

    def test_get_q(self):
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_q(T=self.T0, ignore_q_elec=True, V=self.V0),
            102.34984869477603)

    def test_get_CvoR(self):
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_CvoR(T=self.T0, V=self.V0), 2.9422622359004853)

    def test_get_Cv(self):
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_Cv(T=self.T0, V=self.V0, units='J/mol/K'),
            2.9422622359004853 * c.R('J/mol/K'))
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_Cv(T=self.T0, V=self.V0, units='J/g/K'),
            2.9422622359004853 * c.R('J/mol/K') / self.mw)

    def test_get_CpoR(self):
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_CpoR(T=self.T0, V=self.V0), 3.9422622359004853)

    def test_get_Cp(self):
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_Cp(T=self.T0, V=self.V0, units='J/mol/K'),
            3.9422622359004853 * c.R('J/mol/K'))
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_Cp(T=self.T0, V=self.V0, units='J/g/K'),
            3.9422622359004853 * c.R('J/mol/K') / self.mw)

    def test_get_EoRT(self):
        np.testing.assert_almost_equal(self.CO2_pmutt.get_EoRT(T=self.T0),
                                       -894.97476277965)
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_EoRT(T=self.T0, include_ZPE=True),
            -877.703643641077)

    def test_get_E(self):
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_E(T=self.T0, units='J/mol'),
            -894.97476277965 * c.R('J/mol/K') * self.T0)
        np.testing.assert_almost_equal(self.CO2_pmutt.get_E(T=self.T0,
                                                            units='J/mol',
                                                            include_ZPE=True),
                                       -877.703643641077 * c.R('J/mol/K') *
                                       self.T0,
                                       decimal=2)
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_E(T=self.T0, units='J/g'),
            -894.97476277965 * c.R('J/mol/K') * self.T0 / self.mw)

    def test_get_UoRT(self):
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_UoRT(T=self.T0, V=self.V0), -875.1095022368354)

    def test_get_U(self):
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_U(T=self.T0, V=self.V0, units='J/mol'),
            -875.1095022368354 * c.R('J/mol/K') * self.T0)
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_U(T=self.T0, V=self.V0, units='J/g'),
            -875.1095022368354 * c.R('J/mol/K') * self.T0 / self.mw)

    def test_get_HoRT(self):
        expected_HoRT_CO2 = \
            self.CO2_ASE.get_enthalpy(temperature=self.T0, verbose=False) \
            / c.R('eV/K')/self.T0
        calc_HoRT_CO2 = self.CO2_pmutt.get_HoRT(T=self.T0)
        np.testing.assert_almost_equal(expected_HoRT_CO2, calc_HoRT_CO2, 3)

    def test_get_H(self):
        expected_H_CO2 = \
            self.CO2_ASE.get_enthalpy(temperature=self.T0, verbose=False)
        calc_H_CO2 = self.CO2_pmutt.get_H(T=self.T0, units='eV')
        np.testing.assert_almost_equal(expected_H_CO2, calc_H_CO2, 3)

    def test_get_SoR(self):
        expected_SoR_CO2 = \
            self.CO2_ASE.get_entropy(temperature=self.T0, pressure=self.P0,
                                     verbose=False)/c.R('eV/K')
        calc_SoR_CO2 = self.CO2_pmutt.get_SoR(T=self.T0, V=self.V0)
        np.testing.assert_almost_equal(expected_SoR_CO2, calc_SoR_CO2, 3)

    def test_get_S(self):
        expected_S_CO2 = \
            self.CO2_ASE.get_entropy(temperature=self.T0, pressure=self.P0,
                                     verbose=False)
        calc_S_CO2 = self.CO2_pmutt.get_S(T=self.T0, V=self.V0, units='eV/K')
        np.testing.assert_almost_equal(expected_S_CO2, calc_S_CO2, 1)

    def test_get_FoRT(self):
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_FoRT(T=self.T0, V=self.V0), -900.6031596134445)

    def test_get_F(self):
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_F(T=self.T0, V=self.V0, units='J/mol'),
            -900.6031596134445 * c.R('J/mol/K') * self.T0)
        np.testing.assert_almost_equal(
            self.CO2_pmutt.get_F(T=self.T0, V=self.V0, units='J/g'),
            -900.6031596134445 * c.R('J/mol/K') * self.T0 / self.mw)

    def test_get_GoRT(self):
        expected_GoRT_CO2 = \
            self.CO2_ASE.get_gibbs_energy(temperature=self.T0,
                                          pressure=self.P0,
                                          verbose=False)/c.R('eV/K')/self.T0
        calc_GoRT_CO2 = self.CO2_pmutt.get_GoRT(T=self.T0, V=self.V0)
        np.testing.assert_almost_equal(expected_GoRT_CO2, calc_GoRT_CO2, 3)

    def test_get_G(self):
        expected_G_CO2 = \
            self.CO2_ASE.get_gibbs_energy(temperature=self.T0,
                                          pressure=self.P0,
                                          verbose=False)
        calc_G_CO2 = self.CO2_pmutt.get_G(T=self.T0, V=self.V0, units='eV')
        np.testing.assert_almost_equal(expected_G_CO2, calc_G_CO2, 3)


if __name__ == '__main__':
    unittest.main()
