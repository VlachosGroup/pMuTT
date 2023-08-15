import unittest
import numpy as np
from ase.build import molecule
from pmutt import constants as c
from pmutt import get_molecular_weight
from pmutt.statmech import StatMech, trans, rot, vib, elec
from pmutt.empirical.nasa import Nasa


class TestNasa(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.Nasa_direct = Nasa(name='H2O',
                                elements={
                                    'H': 2,
                                    'O': 1
                                },
                                phase='g',
                                a_low=np.array([
                                    4.04618796E+00, -6.87238823E-04,
                                    2.79722240E-06, -1.42318006E-09,
                                    2.34551159E-13, -3.02826236E+04,
                                    -2.50036531E-01
                                ]),
                                a_high=np.array([
                                    2.41854323E+00, 3.35448922E-03,
                                    -9.66398101E-07, 1.34441829E-10,
                                    -7.18940063E-15, -2.97582484E+04,
                                    8.37839787E+00
                                ]),
                                T_low=100.,
                                T_mid=1610.97,
                                T_high=5000.)

        self.Nasa_direct_dict = {
            'class':
            "<class 'pmutt.empirical.nasa.Nasa'>",
            'name':
            'H2O',
            'elements': {
                'H': 2,
                'O': 1
            },
            'phase':
            'g',
            'a_low': [
                4.04618796E+00, -6.87238823E-04, 2.79722240E-06,
                -1.42318006E-09, 2.34551159E-13, -3.02826236E+04,
                -2.50036531E-01
            ],
            'a_high': [
                2.41854323E+00, 3.35448922E-03, -9.66398101E-07,
                1.34441829E-10, -7.18940063E-15, -2.97582484E+04,
                8.37839787E+00
            ],
            'T_low':
            100.,
            'T_mid':
            1610.97,
            'T_high':
            5000.,
            'notes':
            None,
            'model':
            None,
            'misc_models': [{
                'class': "<class 'pmutt.empirical.GasPressureAdj'>"
            }],
            'cat_site':
            None,
            'n_sites':
            None,
            'smiles':
            None,
            'type':
            'nasa'
        }

        self.Nasa_data = Nasa.from_data(
            name='H2O',
            elements={
                'H': 2,
                'O': 1
            },
            phase='g',
            T=np.array([
                500., 600., 700., 800., 900., 1000., 1100., 1200., 1300.,
                1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200
            ]),
            CpoR=np.array([
                4.238636088, 4.363835667, 4.503924733, 4.654023202,
                4.809813915, 4.967542636, 5.124018051, 5.276611768,
                5.423258319, 5.56245516, 5.693262665, 5.815304137, 5.928750505,
                6.034087273, 6.131819121, 6.222433488, 6.306400563, 6.384173277
            ]),
            T_ref=500.,
            HoRT_ref=-56.49930957,
            SoR_ref=24.84583501)

        self.Nasa_statmech = Nasa.from_model(name='H2O',
                                             elements={
                                                 'H': 2,
                                                 'O': 1
                                             },
                                             phase='g',
                                             model=StatMech,
                                             trans_model=trans.FreeTrans,
                                             n_degrees=3,
                                             vib_model=vib.HarmonicVib,
                                             elec_model=elec.GroundStateElec,
                                             rot_model=rot.RigidRotor,
                                             potentialenergy=-14.2209,
                                             atoms=molecule('H2O'),
                                             symmetrynumber=2,
                                             spin=0,
                                             vib_wavenumbers=np.array(
                                                 [0.47462, 0.46033, 0.19633]),
                                             T_low=100.,
                                             T_mid=1610.97,
                                             T_high=5000.)
        self.mw = get_molecular_weight({'H': 2, 'O': 1})  # g/mol

    def test_get_a(self):
        np.testing.assert_array_equal(self.Nasa_direct.get_a(T=300.),
                                      self.Nasa_direct.a_low)
        np.testing.assert_array_equal(self.Nasa_direct.get_a(T=2000.),
                                      self.Nasa_direct.a_high)
        with self.assertWarns(RuntimeWarning):
            np.testing.assert_array_equal(self.Nasa_direct.get_a(T=6000.),
                                          self.Nasa_direct.a_high)
        with self.assertWarns(RuntimeWarning):
            np.testing.assert_array_equal(self.Nasa_direct.get_a(T=50.),
                                          self.Nasa_direct.a_low)

    def test_get_CpoR(self):
        T = np.array([
            500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400.,
            1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200
        ])
        CpoR_expected = np.array([
            4.238636088, 4.363835667, 4.503924733, 4.654023202, 4.809813915,
            4.967542636, 5.124018051, 5.276611768, 5.423258319, 5.56245516,
            5.693262665, 5.815304137, 5.928750505, 6.034087273, 6.131819121,
            6.222433488, 6.306400563, 6.384173277
        ])
        np.testing.assert_almost_equal(self.Nasa_direct.get_CpoR(T=T[0]),
                                       CpoR_expected[0])
        np.testing.assert_array_almost_equal(self.Nasa_direct.get_CpoR(T=T),
                                             CpoR_expected)

    def test_get_Cp(self):
        T = np.array([
            500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400.,
            1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200
        ])
        Cp_expected = c.R('J/mol/K') *\
            np.array([4.238636088, 4.363835667,
                      4.503924733, 4.654023202, 4.809813915, 4.967542636,
                      5.124018051, 5.276611768, 5.423258319, 5.56245516,
                      5.693262665, 5.815304137, 5.928750505, 6.034087273,
                      6.131819121, 6.222433488, 6.306400563, 6.384173277])
        np.testing.assert_almost_equal(
            self.Nasa_direct.get_Cp(T=T[0], units='J/mol/K'), Cp_expected[0])
        np.testing.assert_array_almost_equal(
            self.Nasa_direct.get_Cp(T=T, units='J/mol/K'), Cp_expected)
        np.testing.assert_almost_equal(
            self.Nasa_direct.get_Cp(T=T[0], units='J/g/K'),
            Cp_expected[0] / self.mw)
        np.testing.assert_array_almost_equal(
            self.Nasa_direct.get_Cp(T=T, units='J/g/K'), Cp_expected / self.mw)

    def test_get_HoRT(self):
        T = np.array([
            500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400.,
            1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200
        ])
        HoRT_expected = np.array([
            -56.49930957, -46.36612849, -39.10913137, -33.64819891,
            -29.38377578, -25.95653237, -23.13812007, -20.77654898,
            -18.76677584, -17.03389718, -15.52306522, -14.19318522,
            -13.01283758, -11.95756475, -11.00803153, -10.14874498,
            -9.367140366, -8.652916499
        ])
        np.testing.assert_almost_equal(self.Nasa_direct.get_HoRT(T=T[0]),
                                       HoRT_expected[0])
        np.testing.assert_array_almost_equal(self.Nasa_direct.get_HoRT(T=T),
                                             HoRT_expected)

    def test_get_H(self):
        T = np.array([
            500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400.,
            1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200
        ])
        H_expected = c.R('J/mol/K')*T *\
            np.array([-56.49930957, -46.36612849, -39.10913137, -33.64819891,
                      -29.38377578, -25.95653237, -23.13812007, -20.77654898,
                      -18.76677584, -17.03389718, -15.52306522, -14.19318522,
                      -13.01283758, -11.95756475, -11.00803153, -10.14874498,
                      -9.367140366, -8.652916499])
        np.testing.assert_almost_equal(self.Nasa_direct.get_H(T=T[0],
                                                              units='J/mol'),
                                       H_expected[0],
                                       decimal=4)
        np.testing.assert_array_almost_equal(self.Nasa_direct.get_H(
            T=T, units='J/mol'),
                                             H_expected,
                                             decimal=4)
        np.testing.assert_almost_equal(self.Nasa_direct.get_H(T=T[0],
                                                              units='J/g'),
                                       H_expected[0] / self.mw,
                                       decimal=4)
        np.testing.assert_array_almost_equal(self.Nasa_direct.get_H(
            T=T, units='J/g'),
                                             H_expected / self.mw,
                                             decimal=4)

    def test_get_SoR(self):
        T = np.array([
            500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400.,
            1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200
        ])
        SoR_expected = np.array([
            24.84583501, 25.62943045, 26.31248017, 26.92360771, 27.48073089,
            27.99565652, 28.47647349, 28.92890014, 29.35709049, 29.76414079,
            30.15242096, 30.52379873, 30.87979567, 31.22169282, 31.55059054,
            31.86744523, 32.17309661, 32.46828858
        ])
        np.testing.assert_almost_equal(self.Nasa_direct.get_SoR(T=T[0]),
                                       SoR_expected[0])
        np.testing.assert_array_almost_equal(self.Nasa_direct.get_SoR(T=T),
                                             SoR_expected)

    def test_get_S(self):
        T = np.array([
            500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400.,
            1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200
        ])
        S_expected = c.R('J/mol/K') *\
            np.array([24.84583501, 25.62943045, 26.31248017, 26.92360771,
                      27.48073089, 27.99565652, 28.47647349, 28.92890014,
                      29.35709049, 29.76414079, 30.15242096, 30.52379873,
                      30.87979567, 31.22169282, 31.55059054, 31.86744523,
                      32.17309661, 32.46828858])
        np.testing.assert_almost_equal(
            self.Nasa_direct.get_S(T=T[0], units='J/mol/K'), S_expected[0])
        np.testing.assert_array_almost_equal(
            self.Nasa_direct.get_S(T=T, units='J/mol/K'), S_expected)
        np.testing.assert_almost_equal(
            self.Nasa_direct.get_S(T=T[0], units='J/g/K'),
            S_expected[0] / self.mw)
        np.testing.assert_array_almost_equal(
            self.Nasa_direct.get_S(T=T, units='J/g/K'), S_expected / self.mw)

    def test_get_GoRT(self):
        T = np.array([
            500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400.,
            1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200
        ])
        HoRT_expected = np.array([
            -56.49930957, -46.36612849, -39.10913137, -33.64819891,
            -29.38377578, -25.95653237, -23.13812007, -20.77654898,
            -18.76677584, -17.03389718, -15.52306522, -14.19318522,
            -13.01283758, -11.95756475, -11.00803153, -10.14874498,
            -9.367140366, -8.652916499
        ])
        SoR_expected = np.array([
            24.84583501, 25.62943045, 26.31248017, 26.92360771, 27.48073089,
            27.99565652, 28.47647349, 28.92890014, 29.35709049, 29.76414079,
            30.15242096, 30.52379873, 30.87979567, 31.22169282, 31.55059054,
            31.86744523, 32.17309661, 32.46828858
        ])
        GoRT_expected = HoRT_expected - SoR_expected
        np.testing.assert_almost_equal(self.Nasa_direct.get_GoRT(T=T[0]),
                                       GoRT_expected[0])
        np.testing.assert_array_almost_equal(self.Nasa_direct.get_GoRT(T=T),
                                             GoRT_expected)

    def test_get_G(self):
        T = np.array([
            500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400.,
            1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200
        ])
        HoRT_expected = np.array([
            -56.49930957, -46.36612849, -39.10913137, -33.64819891,
            -29.38377578, -25.95653237, -23.13812007, -20.77654898,
            -18.76677584, -17.03389718, -15.52306522, -14.19318522,
            -13.01283758, -11.95756475, -11.00803153, -10.14874498,
            -9.367140366, -8.652916499
        ])
        SoR_expected = np.array([
            24.84583501, 25.62943045, 26.31248017, 26.92360771, 27.48073089,
            27.99565652, 28.47647349, 28.92890014, 29.35709049, 29.76414079,
            30.15242096, 30.52379873, 30.87979567, 31.22169282, 31.55059054,
            31.86744523, 32.17309661, 32.46828858
        ])
        GoRT_expected = HoRT_expected - SoR_expected
        np.testing.assert_almost_equal(self.Nasa_direct.get_GoRT(T=T[0]),
                                       GoRT_expected[0])
        np.testing.assert_array_almost_equal(self.Nasa_direct.get_GoRT(T=T),
                                             GoRT_expected)

    def test_get_GoRT_Selements(self):
        T = np.array([
            500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400.,
            1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200
        ])
        HoRT_expected = np.array([
            -56.49930957, -46.36612849, -39.10913137, -33.64819891,
            -29.38377578, -25.95653237, -23.13812007, -20.77654898,
            -18.76677584, -17.03389718, -15.52306522, -14.19318522,
            -13.01283758, -11.95756475, -11.00803153, -10.14874498,
            -9.367140366, -8.652916499
        ])
        SoR_Selements_expected = np.array([
            -3.20842299, -2.42482755, -1.74177783, -1.13065029, -0.57352711,
            -0.05860148,  0.42221549,  0.87464214,  1.30283249,  1.70988279,
            2.09816296,  2.46954073,  2.82553767,  3.16743482,  3.49633254,
            3.81318723,  4.11883861,  4.41403058
        ])
        GoRT_expected = HoRT_expected - SoR_Selements_expected
        np.testing.assert_almost_equal(self.Nasa_direct.
                                       get_GoRT(T=T[0], S_elements=True),
                                       GoRT_expected[0])
        np.testing.assert_array_almost_equal(self.Nasa_direct.
                                             get_GoRT(T=T, S_elements=True),
                                             GoRT_expected)

    def test_to_dict(self):
        self.maxDiff = None
        self.assertEqual(self.Nasa_direct.to_dict(), self.Nasa_direct_dict)

    def test_from_dict(self):
        self.assertEqual(Nasa.from_dict(self.Nasa_direct_dict),
                         self.Nasa_direct)


if __name__ == '__main__':
    unittest.main()
