# -*- coding: utf-8 -*-
"""
Thermochemistry.test_Thermochemistry
Tests for Thermochemistry module
Created on Fri Jul 7 12:31:00 2018
"""
import unittest
import Thermochemistry
import numpy as np
from ase.build import molecule
from ase.thermochemistry import IdealGasThermo

class TestThermochemistry(unittest.TestCase):
	def test_parse_formula(self):
		elements_dict = {
			'Ca': 1,
			'Ti': 1,
			'O': 3,
			}
		self.assertEqual(Thermochemistry.parse_formula('CaTiO3'), elements_dict)

	def test_get_molecular_weight(self):
		elements_dict = {
			'Ca': 1,
			'Ti': 1,
			'O': 3,
			}
		self.assertEqual(Thermochemistry.get_molecular_weight(elements_dict), 135.942)
		self.assertEqual(Thermochemistry.get_molecular_weight('CaTiO3'), 135.942)

		elements_dict_error = {
			'non-existent element': 1,
			'O': 1,
			}
		with self.assertRaises(KeyError):
			Thermochemistry.get_molecular_weight(elements_dict_error)

	def test_get_expected_arguments(self):
		def sum_fn(num1, num2):
			return num1 + num2
		self.assertEqual(Thermochemistry._get_expected_arguments(sum_fn), ('num1', 'num2'))

		class sum_class:
			def __init__(self, num1, num2):
				self.num1 = num1
				self.num2 = num2
			def get_sum(self):
				return self.num1 + self.num2
		self.assertEqual(Thermochemistry._get_expected_arguments(sum_class), ('self', 'num1', 'num2'))

	def test_pass_expected_arguments(self):
		def sum_fn(num1, num2):
			return num1 + num2
		self.assertEqual(Thermochemistry._pass_expected_arguments(sum_fn, **{'num1': 1, 'num2': 2}), 3)

		class sum_class:
			def __init__(self, num1, num2):
				self.num1 = num1
				self.num2 = num2
			def get_sum(self):
				return self.num1 + self.num2
			def __eq__(self, other):
				return self.__dict__ == other.__dict__
		self.assertEqual(Thermochemistry._pass_expected_arguments(sum_class, **{'num1': 1, 'num2': 2}), sum_class(num1 = 1, num2 = 2))

	def setUp(self):
		unittest.TestCase.setUp(self)
		self.BaseThermo = Thermochemistry.BaseThermo(
			name = 'H2O',
			elements = {'H': 2, 'O': 1},
			phase = 'g',
			thermo_model = IdealGasThermo,
			potentialenergy = -14.2209,
			geometry = 'nonlinear',
			atoms = molecule('H2O'),
			symmetrynumber = 2,
			spin = 0,
			vib_energies = np.array([0.47462, 0.46033, 0.19633])
			)

	# def test_BaseThermo(unittest.TestCase):
	# 	pass


if __name__ == '__main__':
	unittest.main()