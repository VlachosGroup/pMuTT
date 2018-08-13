# -*- coding: utf-8 -*-
"""
PyMuTT.models.statmech.heat_capacity

Common heat capacity calculations
"""

import numpy as np
from PyMuTT import constants as c

def get_CvoR_trans(degrees=3):
	"""Calculates the dimensionless translational heat capacity
	
	:math:`\\frac {Cv_{trans}} {R} = \\frac {N} {2}` where N is the degrees of freedom

	Parameters
	----------
		degrees : int
			Degrees of freedom the specie has. Gas-phase species will be 3 degrees whereas surface species may have less depending on the thermodynamic model
	Returns
	-------
		CvoR_trans : float
			Dimensionless translational heat capacity
	Raises
	------
		ValueError
			If degrees > 3 since this is not possible
	"""
	if degrees > 3:
		raise ValueError('Not possible to have more than translational 3 degrees of freedom.')
	return 0.5 * degrees

def get_CvoR_rot(geometry):
	"""Calculates the dimensionless rotational heat capacity

	:math:`\\frac {Cv_{rot}} {R} = 0` if monatomic

	:math:`\\frac {Cv_{rot}} {R} = 1` if linear

	:math:`\\frac {Cv_{rot}} {R} = 1.5` if nonlinear

	Parameters
	----------
		geometry : str
			Geometry of the specie. Currently accepts:

			- monatomic
			- linear
			- nonlinear
	Returns
	-------
		CvoR_rot : float
			Dimensionless rotational heat capacity
	Raises
	------
		ValueError
			If none of the supported units are provided
	"""
	if geometry == 'monatomic':
		return 0.
	elif geometry == 'linear':
		return 1.
	elif geometry == 'nonlinear':
		return 1.5
	else:
		raise ValueError('Geometry {} not supported.'.format(geometry))

def get_CvoR_vib(vib_energies, Ts):
	"""Calculates the dimensionless vibrational heat capacity

	:math:`\\frac {Cv_{vib}} {R} = \\sum_{i=1}^{3n-6} \\bigg(\\frac {\\Theta_{Vi}} {T} \\frac {e^{\\Theta_{Vi}}} {1-e^{\\frac {\\Theta_{Vi}} {T}}} - \\ln ({1 - e^{\\frac {-\\Theta_{Vi}}{T}}})\\bigg)`

	Parameters
	----------
		vib_energies : array-like
			Vibrational energies in eV
		Ts : array-like or scalar
			Temperatures in K to calculate heat capacity
	Returns
	-------
		CvoR_vib : float
			Dimensionless vibrational heat capacity
	"""
	#Check if T is scalar or array-like
	try:
		iter(Ts)
	except TypeError:
		#If scalar
		CvoR_vib = _get_single_CvoR_vib(vib_energies=vib_energies, T=Ts)
	else:
		#If array-like
		CvoR_vib = np.zeros_like(Ts)
		for i, T in enumerate(Ts):
			CvoR_vib[i] = _get_single_CvoR_vib(vib_energies=vib_energies, T=T)
	return CvoR_vib

def _get_single_CvoR_vib(vib_energies, T):
	"""Calculates the dimensionless vibrational heat capacity for a single temperature

	Parameters
	----------
		vib_energies : array-like
			Vibrational energies in eV
		T : float
			Temperature in K to calculate heat capacity
	Returns
	-------
		CvoR_vib : float
			Dimensionless vibrational heat capacity at a single temperature
	"""
	dimensionless_vibs = vib_energies/c.kb('eV/K')/T
	CvoR_vib = np.sum((0.5 * dimensionless_vibs)**2 * (1./np.sinh(0.5 * dimensionless_vibs))**2)
	return CvoR_vib