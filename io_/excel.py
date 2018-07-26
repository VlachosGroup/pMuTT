# -*- coding: utf-8 -*-
"""
Thermochemistry.io
Vlachos group code to read from/write to csv files of particular format.
Created on Fri Jul 7 12:40:00 2018
"""

import numpy as np
import pandas as pd
import os
from ase.io import read
from Thermochemistry import constants as c
from Thermochemistry import parse_formula
from Thermochemistry.models.empirical import BaseThermo

def read_excel(io, skiprows=[1], header=0, delimiter='~', **kwargs):
	"""
	Reads an excel file and returns it as a list of BaseThermo objects

	Parameters
		io - str
			Name of the Excel spreadsheet
		skiprows - list-like
			Rows to skip at the beginning (0-indexed). Default is [1] so comments can be put in that row
		header - int
			Location to find header names (0-index). Default is 0
		delimiter - str
			Delimiter to parse column names. Default is '~'
		**kwargs
			Parameters used by pandas.read_excel
			Not required but some potentially useful parameters include:
				sheet_name - Specify the name of the sheet you're reading.
				converters - Specify how to process certain columns
				dtype - Expected data type. Will be guessed if not specified
				na_values - What strings to interpret as NaN
				convert_float - Converts integral floats to int (i.e. 1.0 --> 1)
	Returns
		list of dict objects
			dict objects that can be used to initialize BaseThermo or BaseThermo-derived object
	"""
	input_data = pd.read_excel(io=io, skiprows=skiprows, header=header, **kwargs)
	excel_path = os.path.dirname(io)
	thermos_out = []
	for row, row_data in input_data.iterrows():
		thermo_data = {'vib_energies': []}
		for col, cell_data in row_data.iteritems():
			#Special parsing instructions
			if pd.isnull(cell_data):
			 	continue
			elif 'element' in col:
				thermo_data = set_element(header=col, value=cell_data, output_structure=thermo_data, delimiter=delimiter)
			elif 'formula' in col:
				thermo_data = set_formula(formula=cell_data, output_structure=thermo_data)
			elif 'atoms' in col:
				thermo_data = set_atoms(path=cell_data, excel_path=excel_path, output_structure=thermo_data)
			elif 'thermo_model' in col:
				thermo_data = set_thermo_model(model=cell_data, output_structure=thermo_data)
			elif 'vib_wavenumber' in col:
				thermo_data = set_vib_wavenumber(value=cell_data, output_structure=thermo_data)				
			else:
				thermo_data[col] = cell_data
		thermos_out.append(thermo_data)
	return thermos_out

def set_element(header, value, output_structure, delimiter = '~'):
	"""
	Parses element header and assigns to output_structure['elements']
	
	Parameters
		header - str
			String containing the element name. Element symbol should be at the end
			e.g. 'element~O'
		value - int
			Amount found in formula
		output_structure - dict
			Structure to assign value. Will assign to output_structure['elements'][element]
		delimiter - str
			Delimiter for element. Element symbol should be at the end
	Returns
		dict
			output_structure with new element added
	"""
	element = header.split(delimiter)[-1]
	try:
		output_structure['elements'][element] = value
	except (NameError, KeyError):
		output_structure['elements'] = {element: value}
	return output_structure

def set_formula(formula, output_structure):
	"""
	Parses stoichiometric formula unit and assigns to output_structure

	Parameters
		formula - str
			Stoichiometric formula unit. e.g. H2O
			Note that an element cannot be specified multiple times. e.g. CH3OH is not supported
		output_structure - dict
			Structure to assign value. Will assign to output_structure['elements']
	Returns
		dict
			output_structure with new elements added
	"""
	elements = parse_formula(formula=formula)
	output_structure['elements'] = elements
	return output_structure

def set_atoms(path, output_structure, excel_path=None):
	"""
	Reads the atoms object ans assigns to output_structure

	Parameters
		path - str
			Location to import atoms object. If relative references used, the path should be relative to excel_path.
			See ase.read for supported formats
		excel_path - str
			Location where excel path is located
		output_structure - dict
			Structure to assign value. Will assign to output_structure['atoms']
	Returns
		output_structure - dict
			output_structure with new thermo model added
	"""
	try:
		output_structure['atoms'] = read(path)
	except FileNotFoundError:
		try:
			output_structure['atoms'] = read(os.path.join(excel_path, path))
		except FileNotFoundError:
			raise FileNotFoundError('If using relative references for atoms files, use a path relative to the spreadsheet imported.')

	return output_structure		

def set_thermo_model(model, output_structure):
	"""
	Imports module and assigns the class to output_structure

	Parameters
		model - str
			Thermodynamic model to import.
			Supported Options:
				IdealGas
				Harmonic
				HinderedRotor
		output_structure - dict
			Structure to assign value. Will assign to output_structure['thermo_model']
	Returns
		output_structure - dict
			output_structure with new thermo model added
	"""
	model = model.lower()
	if model == 'idealgas':
		import Thermochemistry.models.statmech.idealgasthermo as idealgasthermo
		output_structure['thermo_model'] = idealgasthermo.IdealGasThermo
	elif model == 'harmonic':
		import Thermochemistry.models.statmech.harmonicthermo as harmonicthermo
		output_structure['thermo_model'] = harmonicthermo.HarmonicThermo
	elif model == 'hinderedrotor':
		import Thermochemistry.models.statmech.hinderedrotor as hinderedrotor
		output_structure['thermo_model'] = hinderedrotor.HinderedRotor	
	else:
		raise ValueError('Unsupported thermodynamic model, {}. See docstring of Thermochemistry.io_.excel.set_thermo_model for supported models.'.format(model))
	# module = '.'.join(path.split('.')[:-1])
	# thermo_model = path.split('.')[-1]
	# exec('from {} import {}'.format(module, thermo_model))
	# exec("output_structure['thermo_model'] = {}".format(thermo_model))
	return output_structure

def set_vib_wavenumber(value, output_structure):
	"""
	Parses element header and assigns to output_structure['vib_energies']
	
	Parameters
		value - float
			Vibrational frequency in 1/cm
		output_structure - dict
			Structure to assign value. Will assign to output_structure['elements'][element]
	Returns
		dict
			output_structure with new vibration added
	"""
	vib_energy = value*c.c('cm/s')*c.h('eV s')
	try:
		output_structure['vib_energies'].append(vib_energy)
	except (NameError, KeyError):
		output_structure['vib_energies'] = [vib_energy]
	return output_structure