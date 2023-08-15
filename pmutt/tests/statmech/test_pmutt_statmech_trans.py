# -*- coding: utf-8 -*-
"""
pmutt.test_pmutt_model_statmech_trans
Tests for pmutt module
"""
import unittest
import numpy as np
from pmutt.statmech import trans


class TestFreeTrans(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)

        # Using Cl2 as an example
        molecular_weight = 71.
        self.trans_1D = trans.FreeTrans(molecular_weight=molecular_weight,
                                        n_degrees=1)
        self.trans_2D = trans.FreeTrans(molecular_weight=molecular_weight,
                                        n_degrees=2)
        self.trans_3D = trans.FreeTrans(molecular_weight=molecular_weight,
                                        n_degrees=3)
        self.trans_3D_dict = {
            'class': "<class 'pmutt.statmech.trans.FreeTrans'>",
            'molecular_weight': 71.,
            'n_degrees': 3
        }

        self.T = 300.  # K
        self.P = 0.99768  # bar

    def test_get_q(self):
        np.testing.assert_allclose(self.trans_1D.get_q(T=self.T, P=self.P),
                                   3.470587046381861E-15)
        np.testing.assert_allclose(self.trans_2D.get_q(T=self.T, P=self.P),
                                   2.901300856449876E-04)
        np.testing.assert_allclose(self.trans_3D.get_q(T=self.T, P=self.P),
                                   2.425395631097109E+07)

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
        self.assertTrue(
            np.isclose(self.trans_1D.get_SoR(T=self.T, P=self.P),
                       -3.1794507940E+01))
        self.assertTrue(
            np.isclose(self.trans_2D.get_SoR(T=self.T, P=self.P),
                       -6.1452364662E+00))
        self.assertTrue(
            np.isclose(self.trans_3D.get_SoR(T=self.T, P=self.P),
                       1.9504035007E+01))

    def test_get_FoRT(self):
        # Using np.isclose instead of self.assertAlmostEqual since the latter
        # does not compare large floats very well
        self.assertTrue(
            np.isclose(self.trans_1D.get_FoRT(T=self.T, P=self.P),
                       3.2294507940E+01))
        self.assertTrue(
            np.isclose(self.trans_2D.get_FoRT(T=self.T, P=self.P),
                       7.1452364662E+00))
        self.assertTrue(
            np.isclose(self.trans_3D.get_FoRT(T=self.T, P=self.P),
                       -1.8004035007E+01))

    def test_get_GoRT(self):
        # Using np.isclose instead of self.assertAlmostEqual since the latter
        # does not compare large floats very well
        self.assertTrue(
            np.isclose(self.trans_1D.get_GoRT(T=self.T, P=self.P),
                       3.3294507940E+01))
        self.assertTrue(
            np.isclose(self.trans_2D.get_GoRT(T=self.T, P=self.P),
                       8.1452364662E+00))
        self.assertTrue(
            np.isclose(self.trans_3D.get_GoRT(T=self.T, P=self.P),
                       -1.7004035007E+01))

    def test_to_dict(self):
        self.assertEqual(self.trans_3D.to_dict(), self.trans_3D_dict)

    def test_from_dict(self):
        self.assertEqual(trans.FreeTrans.from_dict(self.trans_3D_dict),
                         self.trans_3D)


if __name__ == '__main__':
    unittest.main()
