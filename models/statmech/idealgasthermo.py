# -*- coding: utf-8 -*-
"""
Thermochemistry.thermo_model
Vlachos group code for thermodynamic models.
Created on Tues Jul 10 12:40:00 2018
"""

from ase import thermochemistry
from Thermochemistry import constants as c
from Thermochemistry.models.statmech.heat_capacity import get_CvoR_trans, get_CvoR_vib, get_CvoR_rot

class IdealGasThermo:
	"""
	Ideal Gas Model. Uses ase.thermochemistry.IdealGasThermo to calculate
	enthalpy and entropy
	Attributes
		model - ase.thermochemistry.IdealGasThermo
			Ideal Gas Model
	"""
	def __init__(self, vib_energies, geometry, potentialenergy=0.0, atoms=None, symmetrynumber=None, spin=None, natoms=None):
		self.model = thermochemistry.IdealGasThermo(
			vib_energies = vib_energies,
			geometry = geometry,
			potentialenergy = potentialenergy,
			atoms = atoms, 
			symmetrynumber = symmetrynumber,
			spin = spin,
			natoms = natoms)

	def get_CpoR(Ts):
		"""
		Calculates the dimensionless heat capacity (Cp/R) at a given temperature.
		If you would like to use different behavior from the default, the
		thermo_model used must have the function 'get_CpoR'.

		Parameters
			Ts - array-like or scalar
				Temperature(s) in K
		Returns
			float or (N,) ndarray
				Dimensionless heat capacity (Cp/R)
		"""
		CvoR_trans = get_CvoR_trans(degrees=3.)
		CvoR_vib = get_CvoR_vib(vib_energies=self.model.vib_energies, Ts=Ts)
		CvoR_rot = get_CvoR_rot(geometry=self.model.geometry)
		CvoR_to_CpoR = 1.
		CpoR = CvoR_trans + CvoR_vib + CvoR_rot + CvoR_to_CpoR
		return CpoR

	def get_HoRT(self, T, verbose=False):
		"""
		Returns the dimensionless enthalpy at a given temperature

		Parameters
			T - float
				Temperature in K
			verbose - bool
				Whether a table breaking down each contribution should be printed
		Returns
			float
				Dimensionless heat capacity (H/RT) at the specified temperature
		"""
		return self.model.get_enthalpy(temperature=T, verbose=verbose)/(c.kb('eV/K') * T)

	def get_SoR(self, T, P=1., verbose=False):
		"""
		Returns the dimensionless entropy at a given temperature and pressure

		Parameters
			T - float
				Temperature in K
			P - float
				Pressure in bar
			verbose - bool
				Whether a table breaking down each contribution should be printed
		Returns
			float
				Dimensionless entropy (S/R) at the specified temperature and pressure
		"""
		return self.model.get_entropy(temperature=T, pressure=P*c.convert_unit(from_='bar', to='Pa'), verbose=verbose)

	def get_GoRT(self, T, P=1., verbose=False):
		"""
		Returns the dimensionless Gibbs energy at a given temperature

		Parameters
			T - float
				Temperature in K
			P - float
				Pressure in bar
			verbose - bool
				Whether a table breaking down each contribution should be printed
		Returns
			float
				Dimensionless heat capacity (G/RT) at the specified temperature
		"""
		return self.model.get_gibbs_energy(temperature=T, pressure=P*c.convert_unit(from_='bar', to='Pa'), verbose=verbose)/(c.kb('eV/K') * T)