# -*- coding: utf-8 -*-
"""
pmutt.test_pmutt_model_empirical_lsr
Tests for pmutt module
"""
import unittest
import numpy as np
from pprint import pprint
from pmutt import constants as c
from pmutt.reaction import Reaction
from pmutt.empirical.shomate import Shomate
from pmutt.statmech.lsr import LSR
from pmutt.statmech import StatMech, presets


class TestLSR(unittest.TestCase):
    def setUp(self):
        slope = 0.5
        intercept = 10.  # kcal/mol
        del_E = -1.  # kcal/mol
        del_E_eV = del_E * c.convert_unit(initial='kcal/mol',
                                          final='eV/molecule')
        E_surf = -2.  # kcal/mol
        E_surf_eV = E_surf * c.convert_unit(initial='kcal/mol',
                                            final='eV/molecule')
        E_gas = -3.  # kcal/mol
        E_gas_eV = E_gas * c.convert_unit(initial='kcal/mol',
                                          final='eV/molecule')

        species = {
            'A(g)_shomate':
            Shomate(name='A(g)_shomate',
                    T_low=100.,
                    T_high=500.,
                    a=np.zeros(8)),
            'A(g)_statmech':
            StatMech(),
            '*':
            StatMech(),
            'A*':
            StatMech(U=del_E_eV, H=del_E_eV, **presets['constant']),
            'surf':
            StatMech(U=E_surf_eV, H=E_surf_eV, **presets['constant']),
            'gas':
            StatMech(U=E_gas_eV, H=E_gas_eV, **presets['constant'])
        }
        reaction_shomate = Reaction.from_string('A(g)_shomate + * = A*',
                                                species)
        reaction_statmech = Reaction.from_string('A(g)_statmech + * = A*',
                                                 species)

        self.lsr_const = LSR(slope=slope,
                             intercept=intercept,
                             reaction=del_E,
                             surf_species=E_surf,
                             gas_species=E_gas)
        self.lsr_shomate = LSR(slope=slope,
                               intercept=intercept,
                               reaction=reaction_shomate,
                               surf_species=species['surf'],
                               gas_species=species['gas'])
        self.lsr_statmech = LSR(slope=slope,
                                intercept=intercept,
                                reaction=reaction_statmech,
                                surf_species=species['surf'],
                                gas_species=species['gas'])

    def test_get_UoRT(self):
        self.assertAlmostEqual(self.lsr_const.get_UoRT(),
                               4.5 / c.R('kcal/mol/K') / c.T0('K'),
                               places=2)
        self.assertAlmostEqual(self.lsr_shomate.get_UoRT(),
                               4.5 / c.R('kcal/mol/K') / c.T0('K'),
                               places=2)
        self.assertAlmostEqual(self.lsr_statmech.get_UoRT(),
                               4.5 / c.R('kcal/mol/K') / c.T0('K'),
                               places=2)

    def test_get_HoRT(self):
        self.assertAlmostEqual(self.lsr_const.get_UoRT(),
                               4.5 / c.R('kcal/mol/K') / c.T0('K'),
                               places=2)
        self.assertAlmostEqual(self.lsr_shomate.get_UoRT(),
                               4.5 / c.R('kcal/mol/K') / c.T0('K'),
                               places=2)
        self.assertAlmostEqual(self.lsr_statmech.get_UoRT(),
                               4.5 / c.R('kcal/mol/K') / c.T0('K'),
                               places=2)

    def test_get_FoRT(self):
        self.assertAlmostEqual(self.lsr_const.get_UoRT(),
                               4.5 / c.R('kcal/mol/K') / c.T0('K'),
                               places=2)
        self.assertAlmostEqual(self.lsr_shomate.get_UoRT(),
                               4.5 / c.R('kcal/mol/K') / c.T0('K'),
                               places=2)
        self.assertAlmostEqual(self.lsr_statmech.get_UoRT(),
                               4.5 / c.R('kcal/mol/K') / c.T0('K'),
                               places=2)

    def test_get_GoRT(self):
        self.assertAlmostEqual(self.lsr_const.get_UoRT(),
                               4.5 / c.R('kcal/mol/K') / c.T0('K'),
                               places=2)
        self.assertAlmostEqual(self.lsr_shomate.get_UoRT(),
                               4.5 / c.R('kcal/mol/K') / c.T0('K'),
                               places=2)
        self.assertAlmostEqual(self.lsr_statmech.get_UoRT(),
                               4.5 / c.R('kcal/mol/K') / c.T0('K'),
                               places=2)


if __name__ == '__main__':
    unittest.main()
