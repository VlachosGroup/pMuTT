import unittest
import warnings
import numpy as np
from ase.build import molecule
from PyMuTT.models.statmech import presets, StatMech, trans, rot, vib, elec
from PyMuTT.models.empirical.nasa import Nasa

class TestNasa(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.Nasa = Nasa(
            name = 'H2O',
            elements = {'H': 2, 'O': 1},
            phase = 'g',
            thermo_model = StatMech,
            trans_model = trans.IdealTrans,
            n_degrees = 3,
            vib_model = vib.HarmonicVib,
            elec_model = elec.IdealElec,
            rot_model = rot.RigidRotor,
            potentialenergy = -14.2209,
            geometry = 'nonlinear',
            atoms = molecule('H2O'),
            symmetrynumber = 2,
            spin = 0,
            vib_wavenumbers = np.array([0.47462, 0.46033, 0.19633]),
            a_low = np.array([4.04618796E+00, -6.87238823E-04, 2.79722240E-06, -1.42318006E-09, 2.34551159E-13, -3.02826236E+04, -2.50036531E-01]),
            a_high = np.array([2.41854323E+00, 3.35448922E-03, -9.66398101E-07, 1.34441829E-10, -7.18940063E-15, -2.97582484E+04, 8.37839787E+00]),
            T_low = 100.,
            T_mid = 1610.97, 
            T_high = 5000.
            )

    def test_get_a(self):
        np.testing.assert_array_equal(self.Nasa.get_a(T=self.Nasa.T_mid-500), self.Nasa.a_low)
        np.testing.assert_array_equal(self.Nasa.get_a(T=self.Nasa.T_mid+500), self.Nasa.a_high)
        with self.assertWarns(RuntimeWarning):
            np.testing.assert_array_equal(self.Nasa.get_a(T=self.Nasa.T_high+100), self.Nasa.a_high)
        with self.assertWarns(RuntimeWarning):
            np.testing.assert_array_equal(self.Nasa.get_a(T=self.Nasa.T_low-100), self.Nasa.a_low)

    def test_get_CpoR(self):
        Ts = np.array([500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200])
        CpoR_expected = np.array([4.238636088, 4.363835667, 4.503924733, 4.654023202, 4.809813915, 4.967542636, 5.124018051, 5.276611768, 5.423258319, 5.56245516, 5.693262665, 5.815304137, 5.928750505, 6.034087273, 6.131819121, 6.222433488, 6.306400563, 6.384173277])
        np.testing.assert_almost_equal(self.Nasa.get_CpoR(Ts=Ts[0]), CpoR_expected[0])
        np.testing.assert_array_almost_equal(self.Nasa.get_CpoR(Ts=Ts), CpoR_expected)

    def test_get_HoRT(self):
        Ts = np.array([500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200])
        HoRT_expected = np.array([-56.49930957, -46.36612849, -39.10913137, -33.64819891, -29.38377578, -25.95653237, -23.13812007, -20.77654898, -18.76677584, -17.03389718, -15.52306522, -14.19318522, -13.01283758, -11.95756475, -11.00803153, -10.14874498, -9.367140366, -8.652916499])
        np.testing.assert_almost_equal(self.Nasa.get_HoRT(Ts=Ts[0]), HoRT_expected[0])
        np.testing.assert_array_almost_equal(self.Nasa.get_HoRT(Ts=Ts), HoRT_expected)

    def test_get_SoR(self):
        Ts = np.array([500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200])
        SoR_expected = np.array([24.84583501, 25.62943045, 26.31248017, 26.92360771, 27.48073089, 27.99565652, 28.47647349, 28.92890014, 29.35709049, 29.76414079, 30.15242096, 30.52379873, 30.87979567, 31.22169282, 31.55059054, 31.86744523, 32.17309661, 32.46828858])
        np.testing.assert_almost_equal(self.Nasa.get_SoR(Ts=Ts[0]), SoR_expected[0])
        np.testing.assert_array_almost_equal(self.Nasa.get_SoR(Ts=Ts), SoR_expected)

    def test_get_GoRT(self):
        Ts = np.array([500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200])
        HoRT_expected = np.array([-56.49930957, -46.36612849, -39.10913137, -33.64819891, -29.38377578, -25.95653237, -23.13812007, -20.77654898, -18.76677584, -17.03389718, -15.52306522, -14.19318522, -13.01283758, -11.95756475, -11.00803153, -10.14874498, -9.367140366, -8.652916499])
        SoR_expected = np.array([24.84583501, 25.62943045, 26.31248017, 26.92360771, 27.48073089, 27.99565652, 28.47647349, 28.92890014, 29.35709049, 29.76414079, 30.15242096, 30.52379873, 30.87979567, 31.22169282, 31.55059054, 31.86744523, 32.17309661, 32.46828858])
        GoRT_expected = HoRT_expected - SoR_expected
        np.testing.assert_almost_equal(self.Nasa.get_GoRT(Ts=Ts[0]), GoRT_expected[0])
        np.testing.assert_array_almost_equal(self.Nasa.get_GoRT(Ts=Ts), GoRT_expected)


if __name__ == '__main__':
    unittest.main()