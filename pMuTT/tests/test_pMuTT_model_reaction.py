# -*- coding: utf-8 -*-
"""
pMuTT.test_pMuTT_model_reaction
Tests for pMuTT module
"""
import unittest
import numpy as np
from ase.build import molecule
from pMuTT import constants as c
from pMuTT.models import reaction as rxn
from pMuTT.models.empirical.nasa import Nasa
from pMuTT.models.statmech import StatMech, presets

class TestReaction(unittest.TestCase):
    def setUp(self):
        '''Reactions using Nasa polynomial'''
        self.H2O_nasa = Nasa(name='H2O', T_low=200., T_mid=1000., T_high=3500.,
                             elements={'H': 2, 'O': 1}, 
                             a_low=[4.19864056E+00, -2.03643410E-03, 6.52040211E-06, -5.48797062E-09, 1.77197817E-12, -3.02937267E+04, -8.49032208E-01],
                             a_high=[3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00])
        self.H2_nasa = Nasa(name='H2', T_low=200., T_mid=1000., T_high=3500.,
                            elements={'H': 2},
                            a_low=[2.34433112E+00, 7.98052075E-03, -1.94781510E-05, 2.01572094E-08, -7.37611761E-12, -9.17935173E+02, 6.83010238E-01],
                            a_high=[3.33727920E+00, -4.94024731E-05, 4.99456778E-07, -1.79566394E-10, 2.00255376E-14, -9.50158922E+02, -3.20502331E+00])
        self.O2_nasa = Nasa(name='O2', T_low=200., T_mid=1000., T_high=3500.,
                            elements={'O': 2},
                            a_low=[3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00],
                            a_high=[3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00])
        self.rxn_nasa = rxn.Reaction(
                reactants=[self.H2_nasa, self.O2_nasa],
                reactants_stoich=[1., 0.5],
                products=[self.H2O_nasa],
                products_stoich=[1.])

        self.rxn_nasa_dict ={
               'class': "<class 'pMuTT.models.reaction.Reaction'>",
               'products': [
                    {'T_high': 3500.0,
                     'T_low': 200.0,
                     'T_mid': 1000.0,
                     'a_high': [
                        3.03399249,
                        0.00217691804,
                        -1.64072518e-07,
                        -9.7041987e-11,
                        1.68200992e-14,
                        -30004.2971,
                        4.9667701],
                     'a_low': [
                        4.19864056,
                        -0.0020364341,
                        6.52040211e-06,
                        -5.48797062e-09,
                        1.77197817e-12,
                        -30293.7267,
                        -0.849032208],
                     'class': "<class 'pMuTT.models.empirical.nasa.Nasa'>",
                     'elements': {'H': 2, 'O': 1},
                     'name': 'H2O',
                     'notes': None,
                     'phase': None,
                     'references': None,
                     'statmech_model': None}],
                     'products_stoich': [1.0],
                'reactants': [
                    {'T_high': 3500.0,
                     'T_low': 200.0,
                     'T_mid': 1000.0,
                     'a_high': [
                        3.3372792,
                        -4.94024731e-05,
                        4.99456778e-07,
                        -1.79566394e-10,
                        2.00255376e-14,
                        -950.158922,
                        -3.20502331],
                     'a_low': [
                        2.34433112,
                        0.00798052075,
                        -1.9478151e-05,
                        2.01572094e-08,
                        -7.37611761e-12,
                        -917.935173,
                        0.683010238],
                     'class': "<class 'pMuTT.models.empirical.nasa.Nasa'>",
                     'elements': {'H': 2},
                     'name': 'H2',
                     'notes': None,
                     'phase': None,
                     'references': None,
                     'statmech_model': None},
                    {'T_high': 3500.0,
                     'T_low': 200.0,
                     'T_mid': 1000.0,
                     'a_high': [
                        3.28253784,
                        0.00148308754,
                        -7.57966669e-07,
                        2.09470555e-10,
                        -2.16717794e-14,
                        -1088.45772,
                        5.45323129],
                     'a_low': [
                        3.78245636,
                        -0.00299673416,
                        9.84730201e-06,
                        -9.68129509e-09,
                        3.24372837e-12,
                        -1063.94356,
                        3.65767573],
                     'class': "<class 'pMuTT.models.empirical.nasa.Nasa'>",
                     'elements': {'O': 2},
                     'name': 'O2',
                     'notes': None,
                     'phase': None,
                     'references': None,
                     'statmech_model': None}],
                'reactants_stoich': [1.0, 0.5],
                'transition_state': [None]}

        '''Reactions using StatMech'''
        ideal_gas_param = presets['idealgas']
        self.H2O_sm = StatMech(atoms=molecule('H2O'),
                               symmetrynumber=2,
                               vib_wavenumbers=[3825.434, 3710.2642, 1582.432],
                               potentialenergy=-6.7598,
                               spin=0.,
                               **ideal_gas_param)
        self.H2_sm = StatMech(atoms=molecule('H2'),
                              symmetrynumber=2,
                              vib_wavenumbers=[4306.1793],
                              potentialenergy=-14.2209,
                              spin=0.,
                              **ideal_gas_param)
        self.O2_sm = StatMech(atoms=molecule('O2'),
                              symmetrynumber=2,
                              vib_wavenumbers=[1556.],
                              potentialenergy=-9.862407,
                              spin=1.,
                              **ideal_gas_param)
        # This is an arbitrary transition state for testing
        self.H2O_TS_sm = StatMech(atoms=molecule('H2O'),
                                  symmetrynumber=1.,
                                  vib_wavenumbers=[4000., 3900., 1600.],
                                  potentialenergy=-5.7598,
                                  spin=0.,
                                  **ideal_gas_param)
        self.rxn_sm = rxn.Reaction(
                reactants=[self.H2_sm, self.O2_sm],
                reactants_stoich=[1., 0.5],
                products=[self.H2O_sm],
                products_stoich=[1.],
                transition_state=[self.H2O_TS_sm],
                transition_state_stoich=[1.])

        self.species_dict = {
            'H2O': self.H2O_sm,
            'H2': self.H2_sm,
            'O2': self.O2_sm,
            'H2O_TS': self.H2O_TS_sm
        }


    def test_compare_element_balance(self):
        self.assertIsNone(self.rxn_nasa.check_element_balance())

    def test_get_delta_CvoR(self):
        exp_sm_CvoR = self.H2O_sm.get_CvoR(T=c.T0('K')) \
                      - self.H2_sm.get_CvoR(T=c.T0('K')) \
                      - self.O2_sm.get_CvoR(T=c.T0('K'))*0.5
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_CvoR(T=c.T0('K')),
                exp_sm_CvoR)
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_CvoR(T=c.T0('K'), rev=True),
                -exp_sm_CvoR)


    def test_get_delta_CpoR(self):
        exp_nasa_CpoR = self.H2O_nasa.get_CpoR(T=c.T0('K')) \
                       - self.H2_nasa.get_CpoR(T=c.T0('K')) \
                       - self.O2_nasa.get_CpoR(T=c.T0('K'))*0.5
        exp_sm_CpoR = self.H2O_sm.get_CpoR(T=c.T0('K')) \
                      - self.H2_sm.get_CpoR(T=c.T0('K')) \
                      - self.O2_sm.get_CpoR(T=c.T0('K'))*0.5
        self.assertAlmostEqual(
                self.rxn_nasa.get_delta_CpoR(T=c.T0('K')),
                exp_nasa_CpoR)
        self.assertAlmostEqual(
                self.rxn_nasa.get_delta_CpoR(T=c.T0('K'), rev=True),
                -exp_nasa_CpoR)
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_CpoR(T=c.T0('K')),
                exp_sm_CpoR)
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_CpoR(T=c.T0('K'), rev=True),
                -exp_sm_CpoR)


    def test_get_delta_UoRT(self):
        exp_sm_UoRT = self.H2O_sm.get_UoRT(T=c.T0('K')) \
                      - self.H2_sm.get_UoRT(T=c.T0('K')) \
                      - self.O2_sm.get_UoRT(T=c.T0('K'))*0.5
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_UoRT(T=c.T0('K')),
                exp_sm_UoRT)
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_UoRT(T=c.T0('K'), rev=True),
                -exp_sm_UoRT)


    def test_get_delta_HoRT(self):
        exp_nasa_HoRT = self.H2O_nasa.get_HoRT(T=c.T0('K')) \
                       - self.H2_nasa.get_HoRT(T=c.T0('K')) \
                       - self.O2_nasa.get_HoRT(T=c.T0('K'))*0.5
        exp_sm_HoRT = self.H2O_sm.get_HoRT(T=c.T0('K')) \
                      - self.H2_sm.get_HoRT(T=c.T0('K')) \
                      - self.O2_sm.get_HoRT(T=c.T0('K'))*0.5
        self.assertAlmostEqual(
                self.rxn_nasa.get_delta_HoRT(T=c.T0('K')),
                exp_nasa_HoRT)
        self.assertAlmostEqual(
                self.rxn_nasa.get_delta_HoRT(T=c.T0('K'), rev=True),
                -exp_nasa_HoRT)
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_HoRT(T=c.T0('K')),
                exp_sm_HoRT)
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_HoRT(T=c.T0('K'), rev=True),
                -exp_sm_HoRT)


    def test_get_delta_SoR(self):
        exp_nasa_SoR = self.H2O_nasa.get_SoR(T=c.T0('K')) \
                      - self.H2_nasa.get_SoR(T=c.T0('K')) \
                      - self.O2_nasa.get_SoR(T=c.T0('K'))*0.5
        exp_sm_SoR = self.H2O_sm.get_SoR(T=c.T0('K')) \
                      - self.H2_sm.get_SoR(T=c.T0('K')) \
                      - self.O2_sm.get_SoR(T=c.T0('K'))*0.5
        self.assertAlmostEqual(
                self.rxn_nasa.get_delta_SoR(T=c.T0('K')),
                exp_nasa_SoR)
        self.assertAlmostEqual(
                self.rxn_nasa.get_delta_SoR(T=c.T0('K'), rev=True),
                -exp_nasa_SoR)
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_SoR(T=c.T0('K')),
                exp_sm_SoR)
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_SoR(T=c.T0('K'), rev=True),
                -exp_sm_SoR)


    def test_get_delta_AoRT(self):
        exp_sm_AoRT = self.H2O_sm.get_AoRT(T=c.T0('K')) \
                      - self.H2_sm.get_AoRT(T=c.T0('K')) \
                      - self.O2_sm.get_AoRT(T=c.T0('K'))*0.5
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_AoRT(T=c.T0('K')),
                exp_sm_AoRT)
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_AoRT(T=c.T0('K'), rev=True),
                -exp_sm_AoRT)


    def test_get_delta_GoRT(self):
        exp_nasa_GoRT = self.H2O_nasa.get_GoRT(T=c.T0('K')) \
                       - self.H2_nasa.get_GoRT(T=c.T0('K')) \
                       - self.O2_nasa.get_GoRT(T=c.T0('K'))*0.5
        exp_sm_GoRT = self.H2O_sm.get_GoRT(T=c.T0('K')) \
                      - self.H2_sm.get_GoRT(T=c.T0('K')) \
                      - self.O2_sm.get_GoRT(T=c.T0('K'))*0.5            
        self.assertAlmostEqual(
                self.rxn_nasa.get_delta_GoRT(T=c.T0('K')),
                exp_nasa_GoRT)
        self.assertAlmostEqual(
                self.rxn_nasa.get_delta_GoRT(T=c.T0('K'), rev=True),
                -exp_nasa_GoRT)
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_GoRT(T=c.T0('K')),
                exp_sm_GoRT)
        self.assertAlmostEqual(
                self.rxn_sm.get_delta_GoRT(T=c.T0('K'), rev=True),
                -exp_sm_GoRT)


    def test_get_CvoR_act(self):
        exp_sm_CvoR = self.H2O_TS_sm.get_CvoR(T=c.T0('K')) \
                      - self.H2_sm.get_CvoR(T=c.T0('K')) \
                      - self.O2_sm.get_CvoR(T=c.T0('K'))*0.5
        exp_sm_CvoR_rev = self.H2O_TS_sm.get_CvoR(T=c.T0('K')) \
                          - self.H2O_sm.get_CvoR(T=c.T0('K'))
        self.assertAlmostEqual(
                self.rxn_sm.get_CvoR_act(T=c.T0('K')),
                exp_sm_CvoR)
        self.assertAlmostEqual(
                self.rxn_sm.get_CvoR_act(T=c.T0('K'), rev=True),
                exp_sm_CvoR_rev)


    def test_get_CpoR_act(self):
        exp_sm_CpoR = self.H2O_TS_sm.get_CpoR(T=c.T0('K')) \
                      - self.H2_sm.get_CpoR(T=c.T0('K')) \
                      - self.O2_sm.get_CpoR(T=c.T0('K'))*0.5
        exp_sm_CpoR_rev = self.H2O_TS_sm.get_CpoR(T=c.T0('K')) \
                          - self.H2O_sm.get_CpoR(T=c.T0('K'))
        self.assertAlmostEqual(
                self.rxn_sm.get_CpoR_act(T=c.T0('K')),
                exp_sm_CpoR)
        self.assertAlmostEqual(
                self.rxn_sm.get_CpoR_act(T=c.T0('K'), rev=True),
                exp_sm_CpoR_rev)


    def test_get_UoRT_act(self):
        exp_sm_UoRT = self.H2O_TS_sm.get_UoRT(T=c.T0('K')) \
                      - self.H2_sm.get_UoRT(T=c.T0('K')) \
                      - self.O2_sm.get_UoRT(T=c.T0('K'))*0.5
        exp_sm_UoRT_rev = self.H2O_TS_sm.get_UoRT(T=c.T0('K')) \
                          - self.H2O_sm.get_UoRT(T=c.T0('K'))
        self.assertAlmostEqual(
                self.rxn_sm.get_UoRT_act(T=c.T0('K')),
                exp_sm_UoRT)
        self.assertAlmostEqual(
                self.rxn_sm.get_UoRT_act(T=c.T0('K'), rev=True),
                exp_sm_UoRT_rev)


    def test_get_HoRT_act(self):
        exp_sm_HoRT = self.H2O_TS_sm.get_HoRT(T=c.T0('K')) \
                      - self.H2_sm.get_HoRT(T=c.T0('K')) \
                      - self.O2_sm.get_HoRT(T=c.T0('K'))*0.5
        exp_sm_HoRT_rev = self.H2O_TS_sm.get_HoRT(T=c.T0('K')) \
                          - self.H2O_sm.get_HoRT(T=c.T0('K'))
        self.assertAlmostEqual(
                self.rxn_sm.get_HoRT_act(T=c.T0('K')),
                exp_sm_HoRT)
        self.assertAlmostEqual(
                self.rxn_sm.get_HoRT_act(T=c.T0('K'), rev=True),
                exp_sm_HoRT_rev)


    def test_get_SoR_act(self):
        exp_sm_SoR = self.H2O_TS_sm.get_SoR(T=c.T0('K')) \
                      - self.H2_sm.get_SoR(T=c.T0('K')) \
                      - self.O2_sm.get_SoR(T=c.T0('K'))*0.5
        exp_sm_SoR_rev = self.H2O_TS_sm.get_SoR(T=c.T0('K')) \
                          - self.H2O_sm.get_SoR(T=c.T0('K'))
        self.assertAlmostEqual(
                self.rxn_sm.get_SoR_act(T=c.T0('K')),
                exp_sm_SoR)
        self.assertAlmostEqual(
                self.rxn_sm.get_SoR_act(T=c.T0('K'), rev=True),
                exp_sm_SoR_rev)


    def test_get_AoRT_act(self):
        exp_sm_AoRT = self.H2O_TS_sm.get_AoRT(T=c.T0('K')) \
                      - self.H2_sm.get_AoRT(T=c.T0('K')) \
                      - self.O2_sm.get_AoRT(T=c.T0('K'))*0.5
        exp_sm_AoRT_rev = self.H2O_TS_sm.get_AoRT(T=c.T0('K')) \
                          - self.H2O_sm.get_AoRT(T=c.T0('K'))
        self.assertAlmostEqual(
                self.rxn_sm.get_AoRT_act(T=c.T0('K')),
                exp_sm_AoRT)
        self.assertAlmostEqual(
                self.rxn_sm.get_AoRT_act(T=c.T0('K'), rev=True),
                exp_sm_AoRT_rev)


    def test_get_GoRT_act(self):
        exp_sm_GoRT = self.H2O_TS_sm.get_GoRT(T=c.T0('K')) \
                      - self.H2_sm.get_GoRT(T=c.T0('K')) \
                      - self.O2_sm.get_GoRT(T=c.T0('K'))*0.5
        exp_sm_GoRT_rev = self.H2O_TS_sm.get_GoRT(T=c.T0('K')) \
                          - self.H2O_sm.get_GoRT(T=c.T0('K'))
        self.assertAlmostEqual(
                self.rxn_sm.get_GoRT_act(T=c.T0('K')),
                exp_sm_GoRT)
        self.assertAlmostEqual(
                self.rxn_sm.get_GoRT_act(T=c.T0('K'), rev=True),
                exp_sm_GoRT_rev)

    def test_get_A(self):
        exp_sm_SoR = self.H2O_TS_sm.get_SoR(T=c.T0('K')) \
                     - self.H2_sm.get_SoR(T=c.T0('K')) \
                     - self.O2_sm.get_SoR(T=c.T0('K'))*0.5
        exp_sm_A = np.exp(1.5)*c.kb('J/K')*c.T0('K')/c.h('J s') \
                   *np.exp(exp_sm_SoR)
        exp_sm_SoR_rev = self.H2O_TS_sm.get_SoR(T=c.T0('K')) \
                         - self.H2O_sm.get_SoR(T=c.T0('K'))
        exp_sm_A_rev = np.exp(1.)*c.kb('J/K')*c.T0('K')/c.h('J s') \
                       *np.exp(exp_sm_SoR_rev)
        self.assertAlmostEqual(
                self.rxn_sm.get_A(T=c.T0('K')), exp_sm_A)
        self.assertAlmostEqual(
                self.rxn_sm.get_A(T=c.T0('K'), rev=True), exp_sm_A_rev)

    def test_from_string(self):
        reaction_str = 'H2+0.5O2=H2O_TS=H2O'
        self.assertEqual(rxn.Reaction.from_string(reaction_str=reaction_str, 
                                                  species=self.species_dict),
                         self.rxn_sm)

    def test_to_dict(self):
        self.maxDiff = None
        self.assertEqual(self.rxn_nasa.to_dict(), self.rxn_nasa_dict)


    def test_from_dict(self):
        self.assertEqual(rxn.Reaction.from_dict(self.rxn_nasa_dict), 
                self.rxn_nasa)

class TestHelperReaction(unittest.TestCase):
    def test__parse_reaction(self):
        reaction_str = 'H2+0.5O2=H2O'
        expected_output = (['H2', 'O2'], [1., 0.5], ['H2O'], [1.], None, None)
        self.assertTupleEqual(rxn._parse_reaction(reaction_str=reaction_str),
                              expected_output)
        reaction_str = ' H2 + 0.5 O2 = H2O '
        expected_output = (['H2', 'O2'], [1., 0.5], ['H2O'], [1.], None, None)
        self.assertTupleEqual(rxn._parse_reaction(reaction_str=reaction_str),
                              expected_output)
        reaction_str = ' H2 + 0.5 O2 = H2O_TS = H2O '
        expected_output = (['H2', 'O2'], [1., 0.5], ['H2O'], [1.], 
                           ['H2O_TS'], [1.])
        self.assertTupleEqual(rxn._parse_reaction(reaction_str=reaction_str),
                              expected_output)
    
    def test__parse_reaction_side(self):
        reaction_str = 'H2+0.5O2'
        expected_output = (['H2', 'O2'], [1., 0.5])
        self.assertTupleEqual(
                rxn._parse_reaction_side(reaction_str=reaction_str),
                expected_output)

if __name__ == '__main__':
    unittest.main()