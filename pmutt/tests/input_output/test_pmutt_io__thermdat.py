# -*- coding: utf-8 -*-
"""
Tests for pmutt.io.thermdat module
Created on Fri Jul 7 12:31:00 2018
"""
import unittest
import pmutt.io.thermdat as thermdat
import numpy as np


class TestThermdat(unittest.TestCase):
    def test__get_fields(self):
        lines = [
            '    test 200 ## %^& test100', 'test 200 ## %^& test100    ',
            'test 200 ## %^& test100  \n', '\n  test 200 ## %^& test100',
            '\n test 200 ## %^& test100 \n'
        ]

        # Default delimiter and remove fields
        expected_fields = ['test', '200', '##', '%^&', 'test100']
        for line in lines:
            self.assertListEqual(thermdat._get_fields(line), expected_fields)

        # Changed delimiter to '#'
        expected_fields = ['test200', '%^&test100']
        for line in lines:
            self.assertListEqual(
                thermdat._get_fields(line,
                                     delimiter='#',
                                     remove_fields=['', '\n', ' ']),
                expected_fields)

        # Added 'test' to remove_fields
        expected_fields = ['200', '##', '%^&', '100']
        for line in lines:
            self.assertListEqual(
                thermdat._get_fields(line, remove_fields=['', '\n', 'test']),
                expected_fields)

    def test__is_temperature_header(self):
        false_lines = [
            'THERMO ALL',
            'H2O                     H   2O   1          G200.0     1100.0    493.9         1',
            ' 3.65264072E+00 1.06108515E-03 3.83455580E-08 3.84923664E-10-2.13953966E-13    2',
            '-3.02204928E+04 1.60236266E+00 3.99524709E+00 5.18551442E-04-5.53026360E-06    3',
            ' 1.85895538E-08-1.55138452E-11-3.02807840E+04-7.89384507E-02                   4',
            '100 200 300 400'
        ]
        for line in false_lines:
            self.assertFalse(thermdat._is_temperature_header(line))

        true_lines = [
            '       100       500      1500   \n', '100 500  1500\n',
            '100 500  1500 \n    '
        ]
        for line in true_lines:
            self.assertTrue(thermdat._is_temperature_header(line))

    def test__read_line_num(self):
        lines = [
            'a b c d  1\n', 'a        2 \n', 'a   c    3   ', '         4'
        ]
        expected_values = [1, 2, 3, 4]
        for line, expected_value in zip(lines, expected_values):
            self.assertEqual(thermdat._read_line_num(line), expected_value)

        with self.assertRaises(ValueError):
            thermdat._read_line_num('line ending with letter')

    def test__read_line1(self):
        lines = [
            'H2                      H   2               G200.0     1100.0    493.9         1',
            'H2O                     H   2O   1          G200.0     1100.0    493.9         1',
            'CH3OH                   H   4O   1C   1     G200.0     1100.0    493.9         1',
            'Oxazirene               H   1O   1C   1N   1L200.0     1100.0    493.9         1',
            'Pt(B)   bulk species    Pt  1               S200.0     1100.0    493.9         1'
        ]

        expected_values = [{
            'name': 'H2',
            'elements': {
                'H': 2
            },
            'phase': 'G',
            'T_low': 200.,
            'T_high': 1100.,
            'T_mid': 493.9,
        }, {
            'name': 'H2O',
            'elements': {
                'H': 2,
                'O': 1
            },
            'phase': 'G',
            'T_low': 200.,
            'T_high': 1100.,
            'T_mid': 493.9,
        }, {
            'name': 'CH3OH',
            'elements': {
                'H': 4,
                'O': 1,
                'C': 1
            },
            'phase': 'G',
            'T_low': 200.,
            'T_high': 1100.,
            'T_mid': 493.9,
        }, {
            'name': 'Oxazirene',
            'elements': {
                'H': 1,
                'O': 1,
                'C': 1,
                'N': 1
            },
            'phase': 'L',
            'T_low': 200.,
            'T_high': 1100.,
            'T_mid': 493.9,
        }, {
            'name': 'Pt(B)',
            'elements': {
                'Pt': 1
            },
            'phase': 'S',
            'T_low': 200.,
            'T_high': 1100.,
            'T_mid': 493.9,
            'notes': 'bulk species'
        }]

        for line, expected_value in zip(lines, expected_values):
            self.assertDictEqual(thermdat._read_line1(line), expected_value)

    def test__read_line2(self):
        lines = [
            ' 3.65264072E+00 1.06108515E-03 3.83455580E-08 3.84923664E-10-2.13953966E-13    2',
            '  3.652641E+00   1.061085E-03   3.834556E-08   3.849237E-10  -2.139540E-13     2'
        ]
        expected_values = [
            np.array([
                3.65264072E+00, 1.06108515E-03, 3.83455580E-08, 3.84923664E-10,
                -2.13953966E-13, 0., 0.
            ]),
            np.array([
                3.652641E+00, 1.061085E-03, 3.834556E-08, 3.849237E-10,
                -2.139540E-13, 0., 0.
            ])
        ]

        for line, expected_value in zip(lines, expected_values):
            data = thermdat._read_line2(line, nasa_data={})
            np.testing.assert_allclose(data['a_high'], expected_value)

    def test__read_line3(self):
        lines = [
            ' 3.65264072E+00 1.06108515E-03 3.83455580E-08 3.84923664E-10-2.13953966E-13    3',
            '  3.652641E+00   1.061085E-03   3.834556E-08   3.849237E-10  -2.139540E-13     3'
        ]
        expected_a_high_values = [
            np.array([0., 0., 0., 0., 0., 3.65264072E+00, 1.06108515E-03]),
            np.array([0., 0., 0., 0., 0., 3.652641E+00, 1.061085E-03])
        ]
        expected_a_low_values = [
            np.array([
                3.83455580E-08, 3.84923664E-10, -2.13953966E-13, 0., 0., 0., 0.
            ]),
            np.array(
                [3.834556E-08, 3.849237E-10, -2.139540E-13, 0., 0., 0., 0.])
        ]

        for line, expected_a_high_value, expected_a_low_value in zip(
                lines, expected_a_high_values, expected_a_low_values):
            data = thermdat._read_line3(line,
                                        nasa_data={'a_high': np.zeros(7)})
            np.testing.assert_allclose(data['a_high'], expected_a_high_value)
            np.testing.assert_allclose(data['a_low'], expected_a_low_value)

    def test__read_line4(self):
        lines = [
            ' 3.65264072E+00 1.06108515E-03 3.83455580E-08 3.84923664E-10                   4',
            '  3.652641E+00   1.061085E-03   3.834556E-08   3.849237E-10                    4'
        ]
        expected_values = [
            np.array([
                0., 0., 0., 3.65264072E+00, 1.06108515E-03, 3.83455580E-08,
                3.84923664E-10
            ]),
            np.array([
                0., 0., 0., 3.652641E+00, 1.061085E-03, 3.834556E-08,
                3.849237E-10
            ])
        ]

        for line, expected_value in zip(lines, expected_values):
            data = thermdat._read_line4(line, nasa_data={'a_low': np.zeros(7)})
            np.testing.assert_allclose(data['a_low'], expected_value)


if __name__ == '__main__':
    unittest.main()
