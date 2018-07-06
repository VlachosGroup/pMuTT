# -*- coding: utf-8 -*-
"""
Thermochemistry.test_Thermochemistry
Tests for Thermochemistry module
Created on Fri Jul 7 12:31:00 2018
"""
import unittest
import Thermochemistry

class TestThermochemistry(unittest.TestCase):
	def test_parse_formula(self):
		elements_dict = {
			'Ca': 1,
			'Ti': 1,
			'O': 3,
			}
		self.assertEqual(Thermochemistry.parse_formula('CaTiO3'), elements_dict)

if __name__ == '__main__':
	unittest.main()