"""
Thermochemistry.models.statmech.basethermo
Vlachos group code for thermodynamic models.
Created on Tues Jul 10 12:40:00 2018
"""

import inspect
from Thermochemistry import _pass_expected_arguments

class BaseThermo:
	"""
	The Thermodynamic Parent class. Holds properties of a specie, the 
	statistical-mechanical thermodynamic model.

	Attributes
		name - str
			Name of the specie
		phase - str
			Phase of the specie
				G - gas
				S - surface
		elements - dict
			Composition of the species. Keys of dictionary are elements, 
			values are stoichiometric values in a formula unit
			e.g. CH3OH can be represented as:
			{
				'C': 1,
				'H': 4,
				'O': 1,
			}
		thermo_model - Thermochemistry.thermo_model class or custom class
			Class should have the following methods:
				get_CpoR
				get_HoRT
				get_SoR
				get_GoRT
		T_ref - float
			Reference temperature. Only used for reference species.
		HoRT_ref - float
			Reference dimensionless enthalpy corresponding to T_ref. Only used 
			for reference species.
	"""

	def __init__(self, name, phase, elements, thermo_model, T_ref = None, HoRT_ref = None, **kwargs):
		self.name = name
		self.phase = phase
		self.elements = elements
		self.T_ref = T_ref
		self.HoRT_ref = HoRT_ref
		if inspect.isclass(thermo_model):
			self.thermo_model = _pass_expected_arguments(thermo_model, **kwargs)
		else:
			#If it's an object that has already been initialized
			self.thermo_model = thermo_model

