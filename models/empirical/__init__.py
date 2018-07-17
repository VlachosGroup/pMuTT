"""
Thermochemistry.models.statmech.basethermo
Vlachos group code for thermodynamic models.
Created on Tues Jul 10 12:40:00 2018
"""

import inspect
from matplotlib import pyplot as plt
import numpy as np
from Thermochemistry import _pass_expected_arguments
from Thermochemistry import constants as c

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
		notes - str
			Any additional details you would like to include such as computational set up.
	"""

	def __init__(self, name, phase, elements, thermo_model = None, T_ref = None, HoRT_ref = None, notes = None, **kwargs):
		self.name = name
		self.phase = phase
		self.elements = elements
		self.T_ref = T_ref
		self.HoRT_ref = HoRT_ref
		self.notes = notes
		if inspect.isclass(thermo_model):
			#If you're passing a class. Note that the required arguments will be guessed.
			self.thermo_model = _pass_expected_arguments(thermo_model, **kwargs)
		else:
			#If it's an object that has already been initialized
			self.thermo_model = thermo_model

	def plot_empirical(self, T_low = None, T_high = None, Cp_units = None, H_units = None, S_units = None, G_units = None):
		"""
		Plots the thermodynamic profiles between T_low and T_high using empirical relationship
		Parameters
			T_low - float
				Lower temperature in K. If not specified, T_low attribute used
			T_high - float
				Upper temperature in K. If not specified, T_high attribute used
			Cp_units - str
				Units to plot heat capacity. See Thermochemistry.constants.R for accepted units. 
				If not specified, dimensionless units used.
			H_units - str
				Units to plot enthalpy. See Thermochemistry.constants.R for accepted units but omit the '/K' (e.g. J/mol).
				If not specified, dimensionless units used.
			S_units - str
				Units to plot entropy. See Thermochemistry.constants.R for accepted units. 
				If not specified, dimensionless units used.
			G_units - str
				Units to plot Gibbs free energy. See Thermochemistry.constants.R for accepted units but omit the '/K' (e.g. J/mol).
				If not specified, dimensionless units used.
		"""
		if T_low is None:
			T_low = self.T_low
		if T_high  is None:
			T_high = self.T_high
		Ts = np.linspace(T_low, T_high)

		plt.figure()
		'''
		Heat Capacity
		'''
		plt.subplot(411)
		plt.title('Thermdat Specie: {}'.format(self.name))
		plt.xlabel('Temperature (K)')
		CpoR = self.get_CpoR(Ts=Ts)
		if Cp_units is None:
			plt.ylabel('Cp/R')
			y1 = CpoR
		else:
			plt.ylabel('Cp ({})'.format(Cp_units))
			y1 = CpoR * c.R(Cp_units)
		plt.plot(Ts, y1, 'r-')

		'''
		Enthalpy
		'''
		plt.subplot(412)
		plt.xlabel('Temperature (K)')
		HoRT = self.get_HoRT(Ts=Ts)
		if H_units is None:
			plt.ylabel('H/RT')
			y2 = HoRT
		else:
			plt.ylabel('H ({})'.format(H_units))
			y2 = HoRT * c.R('{}/K'.format(H_units)) * Ts
		plt.plot(Ts, y2, 'g-')

		'''
		Entropy
		'''
		plt.subplot(413)
		plt.xlabel('Temperature (K)')
		SoR = self.get_SoR(Ts=Ts)
		if S_units is None:
			plt.ylabel('S/R')
			y3 = SoR
		else:
			plt.ylabel('S ({})'.format(S_units))
			y3 = SoR * c.R(S_units)
		plt.plot(Ts, y3, 'b-')

		'''
		Gibbs energy
		'''
		plt.subplot(414)
		plt.xlabel('Temperature (K)')
		GoRT = self.get_GoRT(Ts=Ts)
		if G_units is None:
			plt.ylabel('G/RT')
			y4 = GoRT
		else:
			plt.ylabel('G ({})'.format(H_units))
			y4 = GoRT * c.R('{}/K'.format(H_units)) * Ts
		plt.plot(Ts, y4, 'k-')

	def plot_thermo_model(self, T_low = None, T_high = None, Cp_units = None, H_units = None, S_units = None, G_units = None):
		"""
		Plots the thermodynamic profiles between T_low and T_high using statistical mechanics thermodynamic model
		Parameters
			T_low - float
				Lower temperature in K. If not specified, T_low attribute used
			T_high - float
				Upper temperature in K. If not specified, T_high attribute used
			Cp_units - str
				Units to plot heat capacity. See Thermochemistry.constants.R for accepted units. 
				If not specified, dimensionless units used.
			H_units - str
				Units to plot enthalpy. See Thermochemistry.constants.R for accepted units but omit the '/K' (e.g. J/mol).
				If not specified, dimensionless units used.
			S_units - str
				Units to plot entropy. See Thermochemistry.constants.R for accepted units. 
				If not specified, dimensionless units used.
			G_units - str
				Units to plot Gibbs free energy. See Thermochemistry.constants.R for accepted units but omit the '/K' (e.g. J/mol).
				If not specified, dimensionless units used.
		"""
		if T_low is None:
			T_low = self.T_low
		if T_high  is None:
			T_high = self.T_high
		Ts = np.linspace(T_low, T_high)

		plt.figure()
		'''
		Heat Capacity
		'''
		plt.subplot(411)
		plt.title('Thermdat Specie: {}'.format(self.name))
		plt.xlabel('Temperature (K)')
		CpoR = self.thermo_model.get_CpoR(Ts=Ts)
		if Cp_units is None:
			plt.ylabel('Cp/R')
			y1 = CpoR
		else:
			plt.ylabel('Cp ({})'.format(Cp_units))
			y1 = CpoR * c.R(Cp_units)
		plt.plot(Ts, y1, 'r-')

		'''
		Enthalpy
		'''
		plt.subplot(412)
		plt.xlabel('Temperature (K)')
		HoRT = self.thermo_model.get_HoRT(Ts=Ts)
		if H_units is None:
			plt.ylabel('H/RT')
			y2 = HoRT
		else:
			plt.ylabel('H ({})'.format(H_units))
			y2 = HoRT * c.R('{}/K'.format(H_units)) * Ts
		plt.plot(Ts, y2, 'g-')

		'''
		Entropy
		'''
		plt.subplot(413)
		plt.xlabel('Temperature (K)')
		SoR = self.thermo_model.get_SoR(Ts=Ts)
		if S_units is None:
			plt.ylabel('S/R')
			y3 = SoR
		else:
			plt.ylabel('S ({})'.format(S_units))
			y3 = SoR * c.R(S_units)
		plt.plot(Ts, y3, 'b-')

		'''
		Gibbs energy
		'''
		plt.subplot(414)
		plt.xlabel('Temperature (K)')
		GoRT = self.thermo_model.get_GoRT(Ts=Ts)
		if G_units is None:
			plt.ylabel('G/RT')
			y4 = GoRT
		else:
			plt.ylabel('G ({})'.format(H_units))
			y4 = GoRT * c.R('{}/K'.format(H_units)) * Ts
		plt.plot(Ts, y4, 'k-')

	def plot_thermo_model_and_empirical(self, T_low = None, T_high = None, Cp_units = None, H_units = None, S_units = None, G_units = None):
		"""
		Plots the thermodynamic profiles between T_low and T_high
		Parameters
			T_low - float
				Lower temperature in K. If not specified, T_low attribute used
			T_high - float
				Upper temperature in K. If not specified, T_high attribute used
			Cp_units - str
				Units to plot heat capacity. See Thermochemistry.constants.R for accepted units. 
				If not specified, dimensionless units used.
			H_units - str
				Units to plot enthalpy. See Thermochemistry.constants.R for accepted units but omit the '/K' (e.g. J/mol).
				If not specified, dimensionless units used.
			S_units - str
				Units to plot entropy. See Thermochemistry.constants.R for accepted units. 
				If not specified, dimensionless units used.
			G_units - str
				Units to plot Gibbs free energy. See Thermochemistry.constants.R for accepted units but omit the '/K' (e.g. J/mol).
				If not specified, dimensionless units used.
		"""
		if T_low is None:
			T_low = self.T_low
		if T_high  is None:
			T_high = self.T_high
		Ts = np.linspace(T_low, T_high)

		plt.figure()
		'''
		Heat Capacity
		'''
		plt.subplot(411)
		plt.title('Thermdat Specie: {}'.format(self.name))
		plt.xlabel('Temperature (K)')
		CpoR_thermo_model = self.thermo_model.get_CpoR(Ts=Ts)
		CpoR_empirical = self.get_CpoR(Ts=Ts)
		if Cp_units is None:
			plt.ylabel('Cp/R')
			y1_thermo_model = CpoR_thermo_model
			y1_empirical = CpoR_empirical
		else:
			plt.ylabel('Cp ({})'.format(Cp_units))
			y1_thermo_model = CpoR_thermo_model * c.R(Cp_units)
			y1_empirical = CpoR_empirical * c.R(Cp_units)
		plt.plot(Ts, y1_thermo_model, 'r-', label = 'Stat Mech Model')
		plt.plot(Ts, y1_empirical, 'b-', label = 'Empirical Model')
		plt.legend()
		
		'''
		Enthalpy
		'''
		plt.subplot(412)
		plt.xlabel('Temperature (K)')
		HoRT_thermo_model = self.thermo_model.get_HoRT(Ts=Ts)
		HoRT_empirical = self.get_HoRT(Ts=Ts)
		if H_units is None:
			plt.ylabel('H/RT')
			y2_thermo_model = HoRT_thermo_model
			y2_empirical = HoRT_empirical
		else:
			plt.ylabel('H ({})'.format(H_units))
			y2_thermo_model = HoRT_thermo_model * c.R('{}/K'.format(H_units)) * Ts
			y2_empirical = HoRT_empirical * c.R('{}/K'.format(H_units)) * Ts
		plt.plot(Ts, y2_thermo_model, 'r-')
		plt.plot(Ts, y2_empirical, 'b-')

		'''
		Entropy
		'''
		plt.subplot(413)
		plt.xlabel('Temperature (K)')
		SoR_thermo_model = self.thermo_model.get_SoR(Ts=Ts)
		SoR_empirical = self.get_SoR(Ts=Ts)
		if S_units is None:
			plt.ylabel('S/R')
			y3_thermo_model = SoR_thermo_model
			y3_empirical = SoR_empirical
		else:
			plt.ylabel('S ({})'.format(S_units))
			y3_thermo_model = SoR_thermo_model * c.R(S_units)
			y3_empirical = SoR_empirical * c.R(S_units)
		plt.plot(Ts, y3_thermo_model, 'r-')
		plt.plot(Ts, y3_empirical, 'b-')

		'''
		Gibbs energy
		'''
		plt.subplot(414)
		plt.xlabel('Temperature (K)')
		GoRT_thermo_model = self.thermo_model.get_GoRT(Ts=Ts)
		GoRT_empirical = self.get_GoRT(Ts=Ts)
		if G_units is None:
			plt.ylabel('G/RT')
			y4_thermo_model = GoRT_thermo_model
			y4_empirical = GoRT_empirical
		else:
			plt.ylabel('G ({})'.format(H_units))
			y4_thermo_model = GoRT_thermo_model * c.R('{}/K'.format(H_units)) * Ts
			y4_empirical = GoRT_empirical * c.R('{}/K'.format(H_units)) * Ts
		plt.plot(Ts, y4_thermo_model, 'r-')
		plt.plot(Ts, y4_empirical, 'b-')
