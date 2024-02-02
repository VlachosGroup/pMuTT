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
            -312.8864984, -260.02368702, -222.24560239, -193.89448982,
            -171.82759584, -154.15958325, -139.69092339, -127.62210266,
            -117.39973902, -108.62865137, -101.01910275, -94.35382838,
            -88.46670135, -83.22851013, -78.53722772, -74.31120008,
            -70.48428094, -67.00230448
        ])
        SoR_Selements_expected = np.array([
            -3.21032649, -2.42845119, -1.746929, -1.13779778, -0.58349529,
            -0.07236386, 0.40366771, 0.85037837, 1.27202368, 1.67181601,
            2.05223696, 2.41524774, 2.7624351, 3.09511512, 3.41440855,
            3.7212968, 4.01666398, 4.30131811
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
