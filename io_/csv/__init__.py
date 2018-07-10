# -*- coding: utf-8 -*-
"""
Thermochemistry.io
Vlachos group code to read from/write to csv files of particular format.
Created on Fri Jul 7 12:40:00 2018
"""

import numpy as np
import pandas as pd
from ase.io import read
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
	thermos_out = []
	for row, row_data in input_data.iterrows():
		thermo_data = {}
		for col, cell_data in row_data.iteritems():
			#Special parsing instructions
			if 'element' in col:
				thermo_data = set_element(header=col, value=cell_data, output_structure=thermo_data, delimiter=delimiter)
			elif 'formula' in col:
				thermo_data = set_formula(formula=cell_data, output_structure=thermo_data)
			elif 'atoms' in col:
				thermo_data = set_atoms(path=cell_data, output_structure=thermo_data)
			elif 'thermo_model' in col:
				thermo_data = set_thermo_model(path=cell_data, output_structure=thermo_data)
			elif 'vib_energy' in col:
				thermo_data = set_vib_energy(value=cell_data, output_structure=thermo_data)
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

def set_atoms(path, output_structure):
	"""
	Reads the atoms object ans assigns to output_structure

	Parameters
		path - str
			Location to import atoms object. See ase.read for supported formats
		output_structure - dict
			Structure to assign value. Will assign to output_structure['atoms']
	Returns
		output_structure - dict
			output_structure with new thermo model added
	"""
	output_structure['atoms'] = read(path)
	return output_structure		

def set_thermo_model(path, output_structure):
	"""
	Imports module and assigns the class to output_structure

	Parameters
		path - str
			Module to be imported.
			e.g. Thermochemistry.thermo_model.IdealGasThermo
		output_structure - dict
			Structure to assign value. Will assign to output_structure['thermo_model']
	Returns
		output_structure - dict
			output_structure with new thermo model added
	"""
	module = '.'.join(path.split('.')[:-1])
	thermo_model = path.split('.')[-1]
	exec('from {} import {}'.format(module, thermo_model))
	exec("output_structure['thermo_model'] = {}".format(thermo_model))
	return output_structure

def set_vib_energy(value, output_structure):
	"""
	Parses element header and assigns to output_structure['vib_energies']
	
	Parameters
		value - int
			Amount found in formula
		output_structure - dict
			Structure to assign value. Will assign to output_structure['elements'][element]
	Returns
		dict
			output_structure with new vibration added
	"""
	if np.isnan(value):
		return output_structure
	try:
		output_structure['vib_energies'].append(value)
	except (NameError, KeyError):
		output_structure['vib_energies'] = [value]
	return output_structure
