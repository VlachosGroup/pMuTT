# -*- coding: utf-8 -*-
"""
pmutt.test_constants
Tests for pmutt.constants file
Created on Fri Jul 7 12:31:00 2018
"""
import os
import unittest

import numpy as np
import pandas as pd

from pmutt import constants as c


class TestConstants(unittest.TestCase):
    def setUp(self):
        # Read excel sheet with answers
        os.chdir(os.path.dirname(__file__))
        self.ans = pd.read_excel('test_constants.xlsx',
                                 sheet_name='ans',
                                 index_col=0,
                                 header=0)

    def test_R(self):
        # Test R for all units
        self.assertAlmostEqual(c.R('J/mol/K'), self.ans.at['test_R', 0])
        self.assertAlmostEqual(c.R('kJ/mol/K'), self.ans.at['test_R', 1])
        self.assertAlmostEqual(c.R('L kPa/mol/K'), self.ans.at['test_R', 2])
        self.assertAlmostEqual(c.R('cm3 kPa/mol/K'), self.ans.at['test_R', 3])
        self.assertAlmostEqual(c.R('m3 Pa/mol/K'), self.ans.at['test_R', 4])
        self.assertAlmostEqual(c.R('cm3 MPa/mol/K'), self.ans.at['test_R', 5])
        self.assertAlmostEqual(c.R('m3 bar/mol/K'), self.ans.at['test_R', 6])
        self.assertAlmostEqual(c.R('L bar/mol/K'), self.ans.at['test_R', 7])
        self.assertAlmostEqual(c.R('L torr/mol/K'), self.ans.at['test_R', 8])
        self.assertAlmostEqual(c.R('cal/mol/K'), self.ans.at['test_R', 9])
        self.assertAlmostEqual(c.R('kcal/mol/K'), self.ans.at['test_R', 10])
        self.assertAlmostEqual(c.R('L atm/mol/K'), self.ans.at['test_R', 11])
        self.assertAlmostEqual(c.R('cm3 atm/mol/K'), self.ans.at['test_R', 12])
        self.assertAlmostEqual(c.R('eV/K'), self.ans.at['test_R', 13])
        self.assertAlmostEqual(c.R('Eh/K'), self.ans.at['test_R', 14])
        self.assertAlmostEqual(c.R('Ha/K'), self.ans.at['test_R', 15])
        # Test R raises an error when an supported unit is passed
        with self.assertRaises(KeyError):
            c.R('arbitrary unit')

    def test_h(self):
        # Test h for all units (bar=False)
        self.assertAlmostEqual(c.h('J s', bar=False), self.ans.at['test_h', 0])
        self.assertAlmostEqual(c.h('kJ s', bar=False), self.ans.at['test_h',
                                                                   1])
        self.assertAlmostEqual(c.h('eV s', bar=False), self.ans.at['test_h',
                                                                   2])
        self.assertAlmostEqual(c.h('Eh s', bar=False), self.ans.at['test_h',
                                                                   3])
        self.assertAlmostEqual(c.h('Ha s', bar=False), self.ans.at['test_h',
                                                                   4])
        # Test h for all units (bar=True)
        self.assertAlmostEqual(c.h('J s', bar=True), self.ans.at['test_h', 5])
        self.assertAlmostEqual(c.h('kJ s', bar=True), self.ans.at['test_h', 6])
        self.assertAlmostEqual(c.h('eV s', bar=True), self.ans.at['test_h', 7])
        self.assertAlmostEqual(c.h('Eh s', bar=True), self.ans.at['test_h', 8])
        self.assertAlmostEqual(c.h('Ha s', bar=True), self.ans.at['test_h', 9])
        # Test h raises an error when an supported unit is passed
        with self.assertRaises(KeyError):
            c.h('arbitrary unit')

    def test_kb(self):
        # Test kb for all units
        self.assertAlmostEqual(c.kb('J/K'), self.ans.at['test_kb', 0])
        self.assertAlmostEqual(c.kb('kJ/K'), self.ans.at['test_kb', 1])
        self.assertAlmostEqual(c.kb('eV/K'), self.ans.at['test_kb', 2])
        self.assertAlmostEqual(c.kb('cal/K'), self.ans.at['test_kb', 3])
        self.assertAlmostEqual(c.kb('kcal/K'), self.ans.at['test_kb', 4])
        self.assertAlmostEqual(c.kb('Eh/K'), self.ans.at['test_kb', 5])
        self.assertAlmostEqual(c.kb('Ha/K'), self.ans.at['test_kb', 6])
        # Test kb raises an error when an unsupported unit is passed
        with self.assertRaises(KeyError):
            c.kb('arbitrary unit')

    def test_c(self):
        # Test c for all units
        self.assertAlmostEqual(c.c('m/s'), self.ans.at['test_c', 0])
        self.assertAlmostEqual(c.c('cm/s'), self.ans.at['test_c', 1])
        # Test c raises an error when an unsupported unit is passed
        with self.assertRaises(KeyError):
            c.c('arbitrary unit')

    def test_m_e(self):
        # Test m_e for all units
        self.assertAlmostEqual(c.m_e('kg'), self.ans.at['test_m_e', 0])
        self.assertAlmostEqual(c.m_e('g'), self.ans.at['test_m_e', 1])
        self.assertAlmostEqual(c.m_e('amu'), self.ans.at['test_m_e', 2])
        # Test m_e raises an error when an unsupported unit is passed
        with self.assertRaises(ValueError):
            c.m_e('arbitrary unit')

    def test_m_p(self):
        # Test m_p for all units
        self.assertAlmostEqual(c.m_p('kg'), self.ans.at['test_m_p', 0])
        self.assertAlmostEqual(c.m_p('g'), self.ans.at['test_m_p', 1])
        self.assertAlmostEqual(c.m_p('amu'), self.ans.at['test_m_p', 2])
        # Test m_p raises an error when an unsupported unit is passed
        with self.assertRaises(ValueError):
            c.m_p('arbitrary unit')

    def test_P0(self):
        # Test P0 for all units
        self.assertAlmostEqual(c.P0('bar'), self.ans.at['test_P0', 0])
        self.assertAlmostEqual(c.P0('atm'), self.ans.at['test_P0', 1])
        self.assertAlmostEqual(c.P0('Pa'), self.ans.at['test_P0', 2])
        self.assertAlmostEqual(c.P0('kPa'), self.ans.at['test_P0', 3])
        self.assertAlmostEqual(c.P0('MPa'), self.ans.at['test_P0', 4])
        self.assertAlmostEqual(c.P0('psi'), self.ans.at['test_P0', 5])
        self.assertAlmostEqual(c.P0('mmHg'), self.ans.at['test_P0', 6])
        self.assertAlmostEqual(c.P0('torr'), self.ans.at['test_P0', 7])
        # Test P0 raises an error when an unsupported unit is passed
        with self.assertRaises(ValueError):
            c.P0('arbitrary unit')

    def test_T0(self):
        # Test T0 for all units
        self.assertAlmostEqual(c.T0('K'), self.ans.at['test_T0', 0])
        self.assertAlmostEqual(c.T0('C'), self.ans.at['test_T0', 1])
        self.assertAlmostEqual(c.T0('R'), self.ans.at['test_T0', 2])
        self.assertAlmostEqual(c.T0('F'), self.ans.at['test_T0', 3])
        # Test T0 raises an error when an unsupported unit is passed
        with self.assertRaises(ValueError):
            c.T0('arbitrary unit')

    def test_V0(self):
        # Test V0 for all units
        self.assertAlmostEqual(c.V0('m3'), self.ans.at['test_V0', 0])
        self.assertAlmostEqual(c.V0('cm3'), self.ans.at['test_V0', 1])
        self.assertAlmostEqual(c.V0('mL'), self.ans.at['test_V0', 2])
        self.assertAlmostEqual(c.V0('L'), self.ans.at['test_V0', 3])
        # Test V0 raises an error when an unsupported unit is passed
        with self.assertRaises(ValueError):
            c.V0('arbitrary unit')

    def test_convert_unit(self):
        # Test all combinations for temperature conversion
        self.assertAlmostEqual(
            c.convert_unit(c.T0('K'), initial='K', final='C'),
            self.ans.at['test_convert_unit', 1])
        self.assertAlmostEqual(
            c.convert_unit(c.T0('K'), initial='K', final='F'),
            self.ans.at['test_convert_unit', 2])
        self.assertAlmostEqual(
            c.convert_unit(c.T0('K'), initial='K', final='R'),
            self.ans.at['test_convert_unit', 3])

        self.assertAlmostEqual(
            c.convert_unit(c.T0('C'), initial='C', final='K'),
            self.ans.at['test_convert_unit', 0])
        self.assertAlmostEqual(
            c.convert_unit(c.T0('C'), initial='C', final='F'),
            self.ans.at['test_convert_unit', 2])
        self.assertAlmostEqual(
            c.convert_unit(c.T0('C'), initial='C', final='R'),
            self.ans.at['test_convert_unit', 3])

        self.assertAlmostEqual(
            c.convert_unit(c.T0('F'), initial='F', final='K'),
            self.ans.at['test_convert_unit', 0])
        self.assertAlmostEqual(
            c.convert_unit(c.T0('F'), initial='F', final='C'),
            self.ans.at['test_convert_unit', 1])
        self.assertAlmostEqual(
            c.convert_unit(c.T0('F'), initial='F', final='R'),
            self.ans.at['test_convert_unit', 3])

        self.assertAlmostEqual(
            c.convert_unit(c.T0('R'), initial='R', final='K'),
            self.ans.at['test_convert_unit', 0])
        self.assertAlmostEqual(
            c.convert_unit(c.T0('R'), initial='R', final='C'),
            self.ans.at['test_convert_unit', 1])
        self.assertAlmostEqual(
            c.convert_unit(c.T0('R'), initial='R', final='F'),
            self.ans.at['test_convert_unit', 2])

        # Test a unit conversion with multiple-based units
        self.assertAlmostEqual(c.convert_unit(initial='m', final='cm'),
                               self.ans.at['test_convert_unit', 4])
        # Test if error raised when units in different set
        with self.assertRaises(ValueError):
            c.convert_unit(initial='cm', final='J')
        # Test if error raised when unaccepted unit inputted
        with self.assertRaises(ValueError):
            c.convert_unit(initial='arbitrary unit', final='J')
        with self.assertRaises(ValueError):
            c.convert_unit(initial='cm', final='arbitrary unit')

    def test_energy_to_freq(self):
        E_J = c.convert_unit(0.1, initial='eV', final='J')
        np.testing.assert_almost_equal(c.energy_to_freq(E_J),
                                       self.ans.at['test_energy_to_freq', 0],
                                       decimal=1)
        # Decimal set to 1 because expected result and actual result differ by
        # 16th decimal place

    def test_energy_to_temp(self):
        E_J = c.convert_unit(0.1, initial='eV', final='J')
        self.assertAlmostEqual(c.energy_to_temp(E_J),
                               self.ans.at['test_energy_to_temp', 0])

    def test_energy_to_wavenumber(self):
        E_J = c.convert_unit(0.1, initial='eV', final='J')
        self.assertAlmostEqual(c.energy_to_wavenumber(E_J),
                               self.ans.at['test_energy_to_wavenumber', 0])

    def test_freq_to_energy(self):
        self.assertAlmostEqual(c.freq_to_energy(2.42E+13),
                               self.ans.at['test_freq_to_energy', 0])

    def test_freq_to_temp(self):
        self.assertAlmostEqual(c.freq_to_temp(2.42E+13),
                               self.ans.at['test_freq_to_temp', 0])

    def test_freq_to_wavenumber(self):
        self.assertAlmostEqual(c.freq_to_wavenumber(2.42E+13),
                               self.ans.at['test_freq_to_wavenumber', 0])

    def test_inertia_to_temp(self):
        self.assertAlmostEqual(c.inertia_to_temp(7.20E-46),
                               self.ans.at['test_inertia_to_temp', 0])

    def test_temp_to_energy(self):
        self.assertAlmostEqual(c.temp_to_energy(1160.),
                               self.ans.at['test_temp_to_energy', 0])

    def test_temp_to_freq(self):
        self.assertAlmostEqual(c.temp_to_freq(1160.),
                               self.ans.at['test_temp_to_freq', 0])

    def test_temp_to_wavenumber(self):
        self.assertAlmostEqual(c.temp_to_wavenumber(1160.),
                               self.ans.at['test_temp_to_wavenumber', 0])

    def test_wavenumber_to_energy(self):
        self.assertAlmostEqual(c.wavenumber_to_energy(810.),
                               self.ans.at['test_wavenumber_to_energy', 0])

    def test_wavenumber_to_freq(self):
        self.assertAlmostEqual(c.wavenumber_to_freq(810.),
                               self.ans.at['test_wavenumber_to_freq', 0])

    def test_wavenumber_to_inertia(self):
        self.assertAlmostEqual(c.wavenumber_to_inertia(810.),
                               self.ans.at['test_wavenumber_to_inertia', 0])

    def test_wavenumber_to_temp(self):
        self.assertAlmostEqual(c.wavenumber_to_temp(810.),
                               self.ans.at['test_wavenumber_to_temp', 0])

    def test_debye_to_einstein(self):
        self.assertAlmostEqual(c.debye_to_einstein(215.),
                               self.ans.at['test_debye_to_einstein', 0])

    def test_einstein_to_debye(self):
        self.assertAlmostEqual(c.einstein_to_debye(175.),
                               self.ans.at['test_einstein_to_debye', 0])

    # def test_wavenumber_to_temp(self):
    #     self.assertAlmostEqual(c.wavenumber_to_temp(1.), 1.4387773538277204)

    # def test_wavenumber_to_energy(self):
    #     self.assertAlmostEqual(c.wavenumber_to_energy(1.), 1.239841974E-04)

    # def test_wavenumber_to_inertia(self):
    #     self.assertTrue(
    #             np.isclose(c.wavenumber_to_inertia(1.), 2.799275137826E-46))

    # def test_debye_to_einstein(self):
    #     self.assertAlmostEqual(c.debye_to_einstein(215.),
    #                            173.28913505677)

    # def test_einstein_to_debye(self):
    #     self.assertAlmostEqual(c.einstein_to_debye(173.28913505677), 215.)


if __name__ == '__main__':
    unittest.main()
