# -*- coding: utf-8 -*-
"""
pmutt.test_constants
Tests for pmutt.constants file
Created on Fri Jul 7 12:31:00 2018
"""
import unittest
from pmutt import constants as c
import numpy as np


class TestConstants(unittest.TestCase):

    def test_R(self):
        self.assertEqual(c.R('J/mol/K'), 8.3144598)
        with self.assertRaises(KeyError):
            c.R('arbitrary unit')

    def test_h(self):
        self.assertEqual(c.h('J s', bar=False), 6.626070040e-34)
        self.assertEqual(c.h('J s', bar=True), 6.626070040e-34/(2.*np.pi))
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
        with self.assertRaises(ValueError):
            c.m_e('arbitrary unit')

    def test_m_p(self):
        self.assertEqual(c.m_p('amu'), 1.007276466879)
        with self.assertRaises(ValueError):
            c.m_p('arbitrary unit')

    def test_P0(self):
        self.assertEqual(c.P0('bar'), 1.)
        with self.assertRaises(ValueError):
            c.P0('arbitrary unit')

    def test_T0(self):
        self.assertEqual(c.T0('K'), 298.15)
        with self.assertRaises(ValueError):
            c.T0('arbitrary unit')

    def test_convert_unit(self):
        self.assertAlmostEqual(c.convert_unit(num=1., initial='C', final='K'),
                               274.15)
        self.assertAlmostEqual(c.convert_unit(num=1., initial='C', final='F'),
                               33.8)
        self.assertAlmostEqual(c.convert_unit(num=1., initial='C', final='R'), 
                               493.46999999999997)
        self.assertAlmostEqual(c.convert_unit(num=1., initial='K', final='C'),
                               -272.15)
        self.assertAlmostEqual(c.convert_unit(num=1., initial='K', final='F'),
                               -457.87)
        self.assertAlmostEqual(c.convert_unit(num=1., initial='K', final='R'),
                               1.8)
        self.assertAlmostEqual(c.convert_unit(num=1., initial='F', final='C'),
                               -17.22222222222222)
        self.assertAlmostEqual(c.convert_unit(num=1., initial='F', final='K'),
                               255.92777777777778)
        self.assertAlmostEqual(c.convert_unit(num=1., initial='F', final='R'),
                               460.67)

        self.assertAlmostEqual(c.convert_unit(num=1., initial='R', final='C'),
                               -272.59444444444443)
        self.assertAlmostEqual(c.convert_unit(num=1., initial='R', final='K'),
                               0.5555555555555556)
        self.assertAlmostEqual(c.convert_unit(num=1., initial='R', final='F'),
                               -458.67)

        self.assertAlmostEqual(c.convert_unit(initial='m', final='cm'), 100.)
        with self.assertRaises(ValueError):
            c.convert_unit(initial='cm', final='J')
        with self.assertRaises(ValueError):
            c.convert_unit(initial='arbitrary unit', final='J')
        with self.assertRaises(ValueError):
            c.convert_unit(initial='cm', final='arbitrary unit')

    def test_wavenumber_to_temp(self):
        self.assertAlmostEqual(c.wavenumber_to_temp(1.), 1.4387773538277204)

    def test_wavenumber_to_energy(self):
        self.assertAlmostEqual(c.wavenumber_to_energy(1.), 1.239841974E-04)

    def test_wavenumber_to_inertia(self):
        self.assertTrue(
                np.isclose(c.wavenumber_to_inertia(1.), 2.799275137826E-46))

    def test_debye_to_einstein(self):
        self.assertAlmostEqual(c.debye_to_einstein(215.),
                               173.28913505677)

    def test_einstein_to_debye(self):
        self.assertAlmostEqual(c.einstein_to_debye(173.28913505677), 215.)

if __name__ == '__main__':
    unittest.main()
