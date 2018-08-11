# -*- coding: utf-8 -*-
"""
PyMuTT.models.statmech.harmonicthermo
Vlachos group code for Harmonic approximation.
Created on Fri Jul 7 12:40:00 2018
"""

import numpy as np
from ase import thermochemistry
from PyMuTT import constants as c
from PyMuTT.models.statmech.heat_capacity import get_CvoR_trans, get_CvoR_vib, get_CvoR_rot

class HarmonicThermo:
	"""Treats all degrees of freedom harmonically. Uses ase.thermochemistry.HarmonicThermo to calculate	enthalpy and entropy

	Attributes
		model (ase.thermochemistry.HarmonicThermo): The HarmonicThermo object has the following attributes.
			vib_energies ((3N,) ndarray where N is number of atoms): Vibrational energies in eV.
			potentialenergy (float): Potential energy in eV
	"""
	def __init__(self, vib_energies, potentialenergy=0.0):
		self.model = thermochemistry.HarmonicThermo(
			vib_energies = vib_energies,
			potentialenergy = potentialenergy)

	def get_CpoR(self, Ts):
		"""Calculates the dimensionless heat capacity (Cp/R) at a given temperature.

		Args
			Ts (float or (N,) ndarray): Temperature(s) in K
		Returns
			float or (N,) ndarray: Dimensionless heat capacity (Cp/R)
		""" 
		return get_CvoR_vib(vib_energies=self.model.vib_energies, Ts=Ts)

	def get_HoRT(self, Ts, verbose = False):
		"""Returns the dimensionless enthalpy at a given temperature

		Args
			Ts (float or (N,) ndarray): Temperature(s) in K
			verbose (bool): Whether a table breaking down each contribution should be printed
		Returns
			float or (N,) ndarray: Dimensionless heat capacity (H/RT) at the specified temperature
		"""
		try:
			iter(Ts)
		except TypeError:
				HoRT = self.model.get_internal_energy(temperature=Ts, verbose=verbose)/(c.kb('eV/K') * Ts)
		else:
			HoRT = np.zeros_like(Ts)
			for i, T in enumerate(Ts):
				HoRT[i] = self.model.get_internal_energy(temperature=T, verbose=verbose)/(c.kb('eV/K') * T)
		return HoRT

	def get_SoR(self, Ts, verbose=False):
		"""Returns the dimensionless entropy at a given temperature and pressure

		Args
			Ts (float or (N,) ndarray): Temperature(s) in K
			verbose (bool): Whether a table breaking down each contribution should be printed
		Returns
			float or (N,) ndarray: Dimensionless entropy (S/R) at the specified temperature and pressure
		"""
		try:
			iter(Ts)
		except TypeError:
			SoR = self.model.get_entropy(temperature=Ts, verbose=verbose)/c.R('eV/K')
		else:
			SoR = np.zeros_like(Ts)
			for i, T in enumerate(Ts):
				SoR[i] = self.model.get_entropy(temperature=T, verbose=verbose)/c.R('eV/K')
		return SoR

	def get_GoRT(self, Ts, verbose=False):
		"""Returns the dimensionless Gibbs energy at a given temperature

		Args
			Ts (float or (N,) ndarray): Temperature(s) in K
			verbose (bool): Whether a table breaking down each contribution should be printed
		Returns
			float or (N,) ndarray: Dimensionless heat capacity (G/RT) at the specified temperature
		"""
		try:
			iter(Ts)
		except TypeError:
			GoRT = self.model.get_helmholtz_energy(temperature=Ts, verbose=verbose)/(c.kb('eV/K') * Ts)
		else:
			GoRT = np.zeros_like(Ts)
			for i, T in enumerate(Ts):
				GoRT[i] = self.model.get_helmholtz_energy(temperature=T, verbose=verbose)/(c.kb('eV/K') * T)
		return GoRT