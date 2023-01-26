import unittest
import numpy as np
from ase.build import molecule
from pmutt.statmech import StatMech, presets
from pmutt.empirical.nasa import Nasa9, SingleNasa9


class TestNasa(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.Nasa9_direct = Nasa9(
            name='CO2',
            elements={
                'C': 1,
                'O': 2
            },
            phase='g',
            nasas=[
                SingleNasa9(T_low=200.,
                            T_high=1000.,
                            a=np.array([
                                4.943650540E+04, -6.264116010E+02,
                                5.301725240E+00, 2.503813816E-03,
                                -2.127308728E-07, -7.689988780E-10,
                                2.849677801E-13, -4.528198460E+04,
                                -7.048279440E+00
                            ])),
                SingleNasa9(T_low=1000.,
                            T_high=6000.,
                            a=np.array([
                                1.176962419E+05, -1.788791477E+03,
                                8.291523190E+00, -9.223156780E-05,
                                4.863676880E-09, -1.891053312E-12,
                                6.330036590E-16, -3.908350590E+04,
                                -2.652669281E+01
                            ])),
                SingleNasa9(T_low=6000.,
                            T_high=20000.,
                            a=np.array([
                                -1.544423287E+09, 1.016847056E+06,
                                -2.561405230E+02, 3.369401080E-02,
                                -2.181184337E-06, 6.991420840E-11,
                                -8.842351500E-16, -8.043214510E+06,
                                2.254177493E+03
                            ]))
            ])

        H2O_statmech = StatMech(name='H2O',
                                spin=0,
                                symmetrynumber=2,
                                atoms=molecule('H2O'),
                                potentialenergy=-14.2209,
                                vib_wavenumbers=np.array(
                                    [3825.434, 3710.2642, 1582.432]),
                                **presets['idealgas'])

        self.Nasa9_statmech = Nasa9.from_model(name='H2O',
                                               elements={
                                                   'H': 2,
                                                   'O': 1
                                               },
                                               phase='g',
                                               T_low=100.,
                                               T_high=5000.,
                                               model=H2O_statmech)

    def test_get_GoRT_Selements(self):
        T = np.array([
            500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400.,
            1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200
        ])
        HoRT_expected = np.array([
            -313.39206362, -260.30706692, -222.32876415, -193.80083208,
            -171.57810622, -153.77147961, -139.17740223, -126.99263078,
            -116.66063231, -107.78381497, -100.07086638,  -93.30384015,
            -87.31685891,  -81.98191634,  -77.19916276,  -72.89010504,
            -68.99274941,  -65.45806892
        ])
        SoR_Selements_expected = np.array([
            -4.91521161, -3.98339921, -3.13000180, -2.34336165, -1.61301028,
            -0.92967107, -0.28526057,  0.32711071,  0.91313782,  1.47732114,
            2.02296614,  2.55218329,  3.06588796,  3.56380033,  4.04444539,
            4.50515288,  4.94205724,  5.35009766
        ])
        GoRT_expected = HoRT_expected - SoR_Selements_expected
        np.testing.assert_almost_equal(self.Nasa9_statmech.
                                       get_GoRT(T=T[0], S_elements=True),
                                       GoRT_expected[0])
        np.testing.assert_array_almost_equal(self.Nasa9_statmech.
                                             get_GoRT(T=T, S_elements=True),
                                             GoRT_expected)


if __name__ == '__main__':
    unittest.main()
