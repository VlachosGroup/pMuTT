# -*- coding: utf-8 -*-
"""
pmutt.test_pmutt
Tests for pmutt module
Created on Fri Jul 7 12:31:00 2018
"""

import os
import unittest
import pandas as pd
import pmutt


class Testpmutt(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)

        def sum_fn(num1, num2):
            """Function used for testing argument passing"""
            return num1 + num2

        self.sum_fn = sum_fn

        def kwargs_fn(**kwargs):
            return kwargs

        self.kwargs_fn = kwargs_fn

        class sum_class:
            def __init__(self, num1, num2):
                self.num1 = num1
                self.num2 = num2

            def get_sum(self):
                return self.num1 + self.num2

            def get_sum3(self, num3):
                return self.num1 + self.num2 + num3

            def __eq__(self, other):
                return self.__dict__ == other.__dict__

        self.sum_class = sum_class

        class kwargs_class:
            def __init__(self, **kwargs):
                self.kwargs = kwargs

            def get_kwargs(self):
                return self.kwargs

            def __eq__(self, other):
                return self.__dict__ == other.__dict__

        self.kwargs_class = kwargs_class

        os.chdir(os.path.dirname(__file__))
        self.ans = pd.read_excel('test_pmutt.xlsx',
                                 sheet_name='ans',
                                 index_col=0,
                                 header=0)

    def test_parse_formula(self):
        elements_dict = {
            'Ca': 1,
            'Ti': 1,
            'O': 3,
        }
        self.assertEqual(pmutt.parse_formula('CaTiO3'), elements_dict)
        elements_dict = {
            'H': 1,
            'F': 1,
        }
        self.assertEqual(pmutt.parse_formula('HF'), elements_dict)
        elements_dict = {
            'H': 8,
            'C': 3,
        }
        self.assertEqual(pmutt.parse_formula('CH3CH2CH3'), elements_dict)

    def test_get_molecular_weight(self):
        elements_dict = {'Ca': 1, 'Ti': 1, 'O': 3}
        self.assertEqual(pmutt.get_molecular_weight(elements_dict),
                         self.ans.at['test_get_molecular_weight', 0])
        self.assertEqual(pmutt.get_molecular_weight('CaTiO3'),
                         self.ans.at['test_get_molecular_weight', 0])

        elements_dict_error = {'non-existent element': 1, 'O': 1}
        with self.assertRaises(KeyError):
            pmutt.get_molecular_weight(elements_dict_error)

    def test_get_expected_arguments(self):
        self.assertEqual(pmutt._get_expected_arguments(self.sum_fn),
                         ('num1', 'num2'))
        self.assertEqual(pmutt._get_expected_arguments(self.kwargs_fn),
                         tuple())
        self.assertEqual(pmutt._get_expected_arguments(self.sum_class),
                         ('self', 'num1', 'num2'))
        self.assertEqual(pmutt._get_expected_arguments(self.kwargs_class),
                         ('self', ))

    def test_pass_expected_arguments(self):
        self.assertEqual(
            pmutt._pass_expected_arguments(self.sum_fn, **{
                'num1': 1,
                'num2': 2
            }), 3)

        self.assertEqual(
            pmutt._pass_expected_arguments(self.sum_class, **{
                'num1': 1,
                'num2': 2
            }), self.sum_class(num1=1, num2=2))

        self.assertEqual(
            pmutt._pass_expected_arguments(self.sum_class, **{
                'num1': 1,
                'num2': 2,
                'num3': 3
            }), self.sum_class(num1=1, num2=2))

    def test_kwargs_allowed(self):
        self.assertFalse(pmutt._kwargs_allowed(self.sum_fn))
        self.assertTrue(pmutt._kwargs_allowed(self.kwargs_fn))
        self.assertFalse(pmutt._kwargs_allowed(self.sum_class))
        self.assertTrue(pmutt._kwargs_allowed(self.kwargs_class))

    def test_force_pass_arguments(self):
        self.assertEqual(
            pmutt._force_pass_arguments(self.sum_fn, **{
                'num1': 1,
                'num2': 2
            }), 3)

        self.assertEqual(
            pmutt._force_pass_arguments(self.kwargs_fn, **{'num1': 1}),
            {'num1': 1})

        self.assertEqual(
            pmutt._force_pass_arguments(self.sum_class, **{
                'num1': 1,
                'num2': 2
            }), self.sum_class(num1=1, num2=2))

        self.assertEqual(
            pmutt._force_pass_arguments(self.sum_class, **{
                'num1': 1,
                'num2': 2,
                'num3': 3
            }), self.sum_class(num1=1, num2=2))

        self.assertEqual(
            pmutt._force_pass_arguments(self.kwargs_class, **{
                'num1': 1,
                'num2': 2,
                'num3': 3
            }), self.kwargs_class(num1=1, num2=2, num3=3))

    def test_is_iterable(self):
        self.assertTrue(pmutt._is_iterable(list()))
        self.assertTrue(pmutt._is_iterable(tuple()))
        self.assertTrue(pmutt._is_iterable(set()))
        self.assertTrue(pmutt._is_iterable(dict()))
        self.assertFalse(pmutt._is_iterable(''))
        self.assertFalse(pmutt._is_iterable(1))
        self.assertFalse(pmutt._is_iterable(1.))
        self.assertFalse(pmutt._is_iterable(None))

    def test_get_mode_quantity(self):
        sum_obj = self.sum_class(num1=1, num2=2)
        self.assertEqual(
            pmutt._get_mode_quantity(mode=sum_obj,
                                     method_name='get_sum3',
                                     num3=3), 6)

        with self.assertRaises(AttributeError):
            pmutt._get_mode_quantity(mode=sum_obj, method_name='get_prod')

        self.assertEqual(
            pmutt._get_mode_quantity(mode=sum_obj,
                                     method_name='get_prod',
                                     raise_error=False,
                                     raise_warning=False,
                                     default_value=0), 0)

    def test_get_specie_kwargs(self):
        # Tests function does not fail when species' name not present in dict
        in_dict = {'val': 1}
        self.assertDictEqual(pmutt._get_specie_kwargs('test', **in_dict),
                             in_dict)
        # Tests function correctly reads value when species' name is present
        in_dict = {'val1': 1, 'test_kwargs': {'val2': 2}}
        out_dict = {'val1': 1, 'val2': 2}
        self.assertDictEqual(pmutt._get_specie_kwargs('test', **in_dict),
                             out_dict)
        # Tests function correctly reads value when species' name is present
        # and ignores other names
        in_dict = {
            'val1': 1,
            'test_kwargs': {
                'val2': 2
            },
            'other_kwargs': {
                'val3': 3
            }
        }
        out_dict = {'val1': 1, 'val2': 2}
        self.assertDictEqual(pmutt._get_specie_kwargs('test', **in_dict),
                             out_dict)
        # Tests function overwrites default values if species' name is present
        in_dict = {'val1': 1, 'test_kwargs': {'val1': 2}}
        out_dict = {'val1': 2}
        self.assertDictEqual(pmutt._get_specie_kwargs('test', **in_dict),
                             out_dict)

    def test_apply_numpy_operation(self):
        self.assertEqual(
            pmutt._apply_numpy_operation(quantity=[2, 2, 4], operation='sum'),
            8)
        self.assertEqual(
            pmutt._apply_numpy_operation(quantity=[2, 2, 4], operation='prod'),
            16)
        self.assertEqual(
            pmutt._apply_numpy_operation(quantity=[2, 2, 4],
                                         operation='prod',
                                         verbose=True), [2, 2, 4])

    def test_pmutt_list_to_dict(self):
        sum_obj1 = self.sum_class(num1=1, num2=2)
        sum_obj2 = self.sum_class(num1=2, num2=3)
        self.assertDictEqual(
            pmutt.pmutt_list_to_dict([sum_obj1, sum_obj2], key='num1'), {
                1: sum_obj1,
                2: sum_obj2
            })

    def test_format_conditions(self):
        num1 = [1, 2, 3]
        num2 = [2, 4, 6]
        self.assertListEqual(pmutt.format_conditions(num1=num1, num2=num2),
                             [{
                                 'num1': num1[0],
                                 'num2': num2[0]
                             }, {
                                 'num1': num1[1],
                                 'num2': num2[1]
                             }, {
                                 'num1': num1[2],
                                 'num2': num2[2]
                             }])

    def test_get_mass_unit(self):
        self.assertEqual(pmutt._get_mass_unit('kg/s'), 'kg')
        self.assertEqual(pmutt._get_mass_unit('g/s'), 'g')
        self.assertEqual(pmutt._get_mass_unit('cm3/s'), None)

    def test_get_R_adj(self):
        elements = {'H': 2}
        self.assertAlmostEqual(
            pmutt._get_R_adj(units='J/g/K', elements=elements),
            self.ans.at['test_get_R_adj', 0])
        self.assertAlmostEqual(
            pmutt._get_R_adj(units='J/mol/K', elements=elements),
            self.ans.at['test_get_R_adj', 1])

    def test_check_obj(self):
        sum_obj = self.sum_class(num1=1, num2=2)
        self.assertEqual(pmutt._check_obj(sum_obj), sum_obj)
        self.assertEqual(pmutt._check_obj(self.sum_class, num1=1, num2=2),
                         sum_obj)

    def test_check_iterable_attr(self):
        self.assertListEqual(pmutt._check_iterable_attr([1]), [1])
        self.assertTupleEqual(pmutt._check_iterable_attr((1, )), (1, ))
        self.assertDictEqual(pmutt._check_iterable_attr({'1': 1}), {'1': 1})
        self.assertListEqual(pmutt._check_iterable_attr(1), [1])
        self.assertListEqual(pmutt._check_iterable_attr(1.), [1.])
        self.assertEqual(pmutt._check_iterable_attr(None), None)


class Test_pmuttBase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)

        self.pmutt_base = pmutt._pmuttBase()

    def test_eq(self):
        self.assertTrue(self.pmutt_base == self.pmutt_base)
        self.assertFalse(self.pmutt_base == None)

    def test_to_dict(self):
        out_dict = {'class': "<class 'pmutt._pmuttBase'>"}
        self.assertDictEqual(self.pmutt_base.to_dict(), out_dict)

    def test_from_dict(self):
        in_dict = {'class': "<class 'pmutt._pmuttBase'>"}
        self.assertEqual(self.pmutt_base, pmutt._pmuttBase.from_dict(in_dict))


if __name__ == '__main__':
    unittest.main()
