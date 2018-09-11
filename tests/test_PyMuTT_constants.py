# -*- coding: utf-8 -*-
"""
PyMuTT.test_constants
Tests for PyMuTT.constants file
Created on Fri Jul 7 12:31:00 2018
"""
import unittest
from PyMuTT import constants as c
import numpy as np

class TestConstants(unittest.TestCase):

    def test_R(self):
        self.assertEqual(c.R('J/mol/K'), 8.3144598)
        with self.assertRaises(KeyError):
            c.R('arbitrary unit')

    def test_h(self):
        self.assertEqual(c.h('J s', bar = False), 6.626070040e-34)
        self.assertEqual(c.h('J s', bar = True), 6.626070040e-34/(2.*np.pi))
        with self.assertRaises(KeyError):
            c.h('arbitrary unit')

    def test_kb(self):
        self.assertEqual(c.kb('J/K'), 1.38064852e-23)
        with self.assertRaises(KeyError):
            c.kb('arbitrary unit')

    def test_c(self):
        self.assertEqual(c.c('m/s'), 299792458.)
        with self.assertRaises(KeyError):
            c.c('arbitrary unit')

    def test_m_e(self):
        self.assertEqual(c.m_e('amu'), 5.48579909070e-4)
        with self.assertRaises(KeyError):
            c.m_e('arbitrary unit')

    def test_m_p(self):
        self.assertEqual(c.m_p('amu'), 1.007276466879)
        with self.assertRaises(KeyError):
            c.m_p('arbitrary unit')

    def test_P0(self):
        self.assertEqual(c.P0('atm'), 1.)
        with self.assertRaises(KeyError):
            c.P0('arbitrary unit')

    def test_T0(self):
        self.assertEqual(c.T0('K'), 298.15)
        with self.assertRaises(KeyError):
            c.T0('arbitrary unit')

    def test_convert_unit(self):
        self.assertEqual(c.convert_unit(num = 0., from_ = 'C', to = 'K'), 273.15)
        self.assertEqual(c.convert_unit(from_ = 'm', to = 'cm'), 100.)
        with self.assertRaises(ValueError):
            c.convert_unit(from_ = 'cm', to = 'J')
        with self.assertRaises(ValueError):
            c.convert_unit(from_ = 'arbitrary unit', to = 'J')
        with self.assertRaises(ValueError):
            c.convert_unit(from_ = 'cm', to = 'arbitrary unit')

    def test_wavenumber_to_temp(self):
        self.assertAlmostEqual(c.wavenumber_to_temp(1.), 1.4387773538277204)

    def test_wavenumber_to_energy(self):
        self.assertAlmostEqual(c.wavenumber_to_energy(1.), 1.239841974E-04)

    def test_wavenumber_to_inertia(self):
        self.assertTrue( \
                np.isclose(c.wavenumber_to_inertia(1.), 2.799275137826E-46))

if __name__ == '__main__':
    unittest.main()