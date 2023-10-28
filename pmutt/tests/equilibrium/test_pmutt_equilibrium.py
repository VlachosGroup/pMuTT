# -*- coding: utf-8 -*-
"""
pmutt.test_pmutt
Tests for pmutt module
Created on Wed Mar 8 2023
"""
import unittest
import os
import numpy.testing as npt
import numpy as np
from pmutt import equilibrium

class TestExamples(unittest.TestCase):

    def test_equilibrium_comp(self):
        os.chdir(os.path.dirname(__file__))
        filepath = 'thermdat_equilibrium_unittest.txt'
        network = {'CH3CH2CH3': 1, 'H2O': 0.7, 'H2': 0, 'CH2CHCH3': 0,
                   'CH4': 0, 'CHCH': 0, 'CH2CH2': 0, 'CH3CH3': 0,
                   'CO2': 0, 'CO': 0}
        equil1 = equilibrium.Equilibrium.from_thermdat(filepath, network)
        sol = equil1.get_net_comp(T=500, P=1.0)
        species_list = ['CH3CH2CH3', 'H2O', 'H2', 'CH2CHCH3', 'CH4',
                        'CHCH', 'CH2CH2', 'CH3CH3', 'CO2', 'CO']
        expected_moles = np.array([2.63077747e-02, 5.00399270e-09,
                                   2.55201909e-06, 1.46783124e-02,
                                   2.20446579e+00, 8.30523843e-10,
                                   1.05881731e-03, 4.65607288e-02,
                                   1.22663143e-01, 4.54673709e-01])
        expected_mol_frac = np.array([9.16516003e-03, 1.74330191e-09,
                                      8.89077988e-07, 5.11366257e-03,
                                      7.67996610e-01, 2.89339711e-10,
                                      3.68873090e-04, 1.62209285e-02,
                                      4.27336537e-02, 1.58400221e-01])
        self.assertEqual(species_list, sol.species)
        npt.assert_array_almost_equal(expected_moles, sol.moles)
        npt.assert_array_almost_equal(expected_mol_frac, sol.mole_frac)
        self.assertEqual(sol.T, 500)
        self.assertEqual(sol.P, 1.0)


if __name__ == '__main__':
    unittest.main()
