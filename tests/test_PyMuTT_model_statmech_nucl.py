# -*- coding: utf-8 -*-
"""
PyMuTT.test_PyMuTT_model_statmech_nucl
Tests for PyMuTT module
"""
import unittest
from PyMuTT.models.statmech import nucl

class TestIdealNucl(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.nuclear = nucl.IdealNucl()
        self.nuclear_dict = {
            'class': "<class 'PyMuTT.models.statmech.nucl.IdealNucl'>"
        }

    def test_get_q(self):
        self.assertEqual(self.nuclear.get_q(), 1.)

    def test_get_CvoR(self):
        self.assertEqual(self.nuclear.get_CvoR(), 0.)

    def test_get_CpoR(self):
        self.assertEqual(self.nuclear.get_CpoR(), 0.)

    def test_get_UoRT(self):
        self.assertEqual(self.nuclear.get_UoRT(), 0.)

    def test_get_HoRT(self):
        self.assertEqual(self.nuclear.get_HoRT(), 0.)

    def test_get_SoR(self):
        self.assertEqual(self.nuclear.get_SoR(), 0.)

    def test_get_AoRT(self):
        self.assertEqual(self.nuclear.get_AoRT(), 0.)

    def test_get_GoRT(self):
        self.assertEqual(self.nuclear.get_GoRT(), 0.)

    def test_to_dict(self):
        self.assertEqual(self.nuclear.to_dict(), self.nuclear_dict)

    def test_from_dict(self):
        self.assertEqual(nucl.IdealNucl.from_dict(self.nuclear_dict), 
                self.nuclear)

if __name__ == '__main__':
    unittest.main()