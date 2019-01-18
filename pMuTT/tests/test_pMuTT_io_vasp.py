import os
import unittest
from pMuTT.io_ import vasp


class TestVasp(unittest.TestCase):
    def test_set_vib_wavenumbers_from_outcar(self):
        test_file = os.path.join(os.path.dirname(__file__), 'test_OUTCAR')
        incorrect_file = 'this_is_not_a_file.txt'
        out_dict = dict()
        expected_wavenumber_nocutoff = {'vib_wavenumbers':
                                        [3821.717493, 3703.479948,
                                         1535.727129, 115.153397,
                                         105.380772, 64.404843]}
        cutoff = 100
        expected_wavenumber_cutoff = {'vib_wavenumbers':
                                      [3821.717493, 3703.479948,
                                       1535.727129, 115.153397,
                                       105.380772]}
        # check parsing with no min_frequency_cutoff
        self.assertEqual(vasp.set_vib_wavenumbers_from_outcar(test_file,
                                                              out_dict, 0),
                         expected_wavenumber_nocutoff)
        # check min_frequency_cutoff working
        self.assertEqual(vasp.set_vib_wavenumbers_from_outcar(test_file,
                                                              out_dict,
                                                              cutoff),
                         expected_wavenumber_cutoff)
        # check file opening exception handling
        with self.assertRaises(FileNotFoundError):
            vasp.set_vib_wavenumbers_from_outcar(incorrect_file, out_dict,
                                                 cutoff)

    def test_get_vib_wavenumber_from_line(self):
        expected_frequency = 3821.717493
        sample_line = '   1 f  =  114.572212 THz   719.878437 2PiTHz ' \
                      '3821.717493 cm-1   473.832750 meV'
        self.assertEqual(vasp.get_vib_wavenumber_from_line(sample_line),
                         expected_frequency)


if __name__ == '__main__':
    unittest.main()
