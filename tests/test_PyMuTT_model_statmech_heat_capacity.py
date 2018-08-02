import unittest
import numpy as np
from PyMuTT import constants as c
from PyMuTT.models.statmech import heat_capacity

class TestHeatCapacity(unittest.TestCase):
	def test_get_CvoR_trans(self):
		self.assertEqual(heat_capacity.get_CvoR_trans(degrees=1), 0.5)
		self.assertEqual(heat_capacity.get_CvoR_trans(), 1.5)

	def test_get_CvoR_rot(self):
		self.assertEqual(heat_capacity.get_CvoR_rot(geometry='monatomic'), 0.)
		self.assertEqual(heat_capacity.get_CvoR_rot(geometry='linear'), 1.)
		self.assertEqual(heat_capacity.get_CvoR_rot(geometry='nonlinear'), 1.5)		
		with self.assertRaises(ValueError):
			heat_capacity.get_CvoR_rot(geometry='unsupported option')

	def test_get_CvoR_vib(self):
		#Vibrational frequencies correspond to H2O
		vib_energies = np.array([5360., 5160., 2290.]) * c.kb('eV/K')
		self.assertAlmostEqual(heat_capacity._get_single_CvoR_vib(vib_energies=vib_energies, T=300.), 0.02824711469596)
		self.assertAlmostEqual(heat_capacity.get_CvoR_vib(vib_energies=vib_energies, Ts=300.), 0.02824711469596)
		np.testing.assert_almost_equal(heat_capacity.get_CvoR_vib(vib_energies=vib_energies, Ts=np.array([300., 400., 500.])), np.array([0.02824711469596, 0.10834772440392, 0.22564243671432]))

if __name__ == '__main__':
	unittest.main()