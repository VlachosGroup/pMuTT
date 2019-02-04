# -*- coding: utf-8 -*-
"""
pMuTT.test_pMuTT
Tests for pMuTT module
Created on Fri Jul 7 12:31:00 2018
"""
import unittest
import pMuTT


class TestpMuTT(unittest.TestCase):
    def test_parse_formula(self):
        elements_dict = {'Ca': 1, 'Ti': 1, 'O': 3, }
        self.assertEqual(pMuTT.parse_formula('CaTiO3'), elements_dict)
        elements_dict = {'H': 1, 'F': 1, }
        self.assertEqual(pMuTT.parse_formula('HF'), elements_dict)
        elements_dict = {'H': 8, 'C': 3, }
        self.assertEqual(pMuTT.parse_formula('CH3CH2CH3'), elements_dict)

    def test_get_molecular_weight(self):
        elements_dict = {
            'Ca': 1,
            'Ti': 1,
            'O': 3,
            }
        self.assertEqual(pMuTT.get_molecular_weight(elements_dict), 135.942)
        self.assertEqual(pMuTT.get_molecular_weight('CaTiO3'), 135.942)

        elements_dict_error = {
            'non-existent element': 1,
            'O': 1,
            }
        with self.assertRaises(KeyError):
            pMuTT.get_molecular_weight(elements_dict_error)

    def test_get_expected_arguments(self):
        def sum_fn(num1, num2):
            return num1 + num2
        self.assertEqual(pMuTT._get_expected_arguments(sum_fn),
                         ('num1', 'num2'))

        class sum_class:
            def __init__(self, num1, num2):
                self.num1 = num1
                self.num2 = num2

            def get_sum(self):
                return self.num1 + self.num2
        self.assertEqual(pMuTT._get_expected_arguments(sum_class),
                         ('self', 'num1', 'num2'))

    def test_pass_expected_arguments(self):
        def sum_fn(num1, num2):
            return num1 + num2
        self.assertEqual(pMuTT._pass_expected_arguments
                         (sum_fn, **{'num1': 1, 'num2': 2}), 3)

        class sum_class:
            def __init__(self, num1, num2):
                self.num1 = num1
                self.num2 = num2

            def get_sum(self):
                return self.num1 + self.num2

            def __eq__(self, other):
                return self.__dict__ == other.__dict__
        self.assertEqual(pMuTT._pass_expected_arguments
                         (sum_class, **{'num1': 1, 'num2': 2}),
                         sum_class(num1=1, num2=2))
        self.assertEqual(pMuTT._pass_expected_arguments
                         (sum_class, **{'num1': 1, 'num2': 2, 'num3': 3}),
                         sum_class(num1=1, num2=2))


if __name__ == '__main__':
    unittest.main()
