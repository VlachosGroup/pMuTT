# -*- coding: utf-8 -*-
"""
PyMuTT.test_PyMuTT_model_statmech_trans
Tests for PyMuTT module
"""
import unittest
import numpy as np
from PyMuTT.models.statmech import trans

class TestIdealTrans(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)

        #Using Cl2 as an example
        molecular_weight = 71.
        self.trans_1D = trans.IdealTrans(molecular_weight=molecular_weight,
                                         n_degrees=1)
        self.trans_2D = trans.IdealTrans(molecular_weight=molecular_weight,
                                         n_degrees=2)
        self.trans_3D = trans.IdealTrans(molecular_weight=molecular_weight,
                                         n_degrees=3)
        self.T = 300 # K
        self.P = 0.99768 # bar

    def test_get_q(self):
        # Using np.isclose instead of self.assertAlmostEqual since the latter
        # does not compare large floats very well
        self.assertTrue(np.isclose(self.trans_1D.get_q(T=self.T, P=self.P), 
                                   2090036406.020292))
        self.assertTrue(np.isclose(self.trans_2D.get_q(T=self.T, P=self.P), 
                                   1.747204243477979e+20))
        self.assertTrue(np.isclose(self.trans_3D.get_q(T=self.T, P=self.P), 
                                   1.4606074131695383e+31))        

    def test_get_CvoR(self):
        self.assertEqual(self.trans_1D.get_CvoR(), 0.5)
        self.assertEqual(self.trans_2D.get_CvoR(), 1.)
        self.assertEqual(self.trans_3D.get_CvoR(), 1.5)

    def test_get_CpoR(self):
        self.assertEqual(self.trans_1D.get_CpoR(), 1.5)
        self.assertEqual(self.trans_2D.get_CpoR(), 2.)
        self.assertEqual(self.trans_3D.get_CpoR(), 2.5)

    def test_get_UoRT(self):
        self.assertEqual(self.trans_1D.get_UoRT(), 0.5)
        self.assertEqual(self.trans_2D.get_UoRT(), 1.)
        self.assertEqual(self.trans_3D.get_UoRT(), 1.5)

    def test_get_HoRT(self):
        self.assertEqual(self.trans_1D.get_HoRT(), 1.5)
        self.assertEqual(self.trans_2D.get_HoRT(), 2.)
        self.assertEqual(self.trans_3D.get_HoRT(), 2.5)

    def test_get_SoR(self):
        # Using np.isclose instead of self.assertAlmostEqual since the latter
        # does not compare large floats very well
        self.assertTrue(np.isclose(self.trans_1D.get_SoR(T=self.T, P=self.P), 
                                   -3.1794507940E+01))
        self.assertTrue(np.isclose(self.trans_2D.get_SoR(T=self.T, P=self.P), 
                                   -6.1452364662E+00))
        self.assertTrue(np.isclose(self.trans_3D.get_SoR(T=self.T, P=self.P), 
                                   1.9504035007E+01))        

    def test_get_AoRT(self):
        # Using np.isclose instead of self.assertAlmostEqual since the latter
        # does not compare large floats very well
        self.assertTrue(np.isclose(self.trans_1D.get_AoRT(T=self.T, P=self.P), 
                                   3.2294507940E+01))
        self.assertTrue(np.isclose(self.trans_2D.get_AoRT(T=self.T, P=self.P), 
                                   7.1452364662E+00))
        self.assertTrue(np.isclose(self.trans_3D.get_AoRT(T=self.T, P=self.P), 
                                   -1.8004035007E+01))        

    def test_get_GoRT(self):
        # Using np.isclose instead of self.assertAlmostEqual since the latter
        # does not compare large floats very well
        self.assertTrue(np.isclose(self.trans_1D.get_GoRT(T=self.T, P=self.P), 
                                   3.3294507940E+01))
        self.assertTrue(np.isclose(self.trans_2D.get_GoRT(T=self.T, P=self.P), 
                                   8.1452364662E+00))
        self.assertTrue(np.isclose(self.trans_3D.get_GoRT(T=self.T, P=self.P), 
                                   -1.7004035007E+01))        

if __name__ == '__main__':
    unittest.main()