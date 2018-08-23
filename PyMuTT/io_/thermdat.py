# -*- coding: utf-8 -*-
"""
PyMuTT.io_.thermdat

Read from/write to thermdat files.
"""

import numpy as np
from datetime import datetime
from PyMuTT import constants as c
from PyMuTT.models.empirical.nasa import Nasa

def read_thermdat(filename):
	"""Directly read thermdat file that is in the Chemkin format
	
	Parameters
	----------
		filename : str
			Input filename
	Returns
	-------
		Nasas : list of ``PyMuTT.models.empirical.nasa.Nasa``
	Raises
	------
		FileNotFoundError
			If the file isn't found.
		IOError
			Invalid line number found.
	"""

	species = []
	with open(filename, 'r') as f_ptr:
		for line in f_ptr:
			'''
			Lines to skip
			'''
			#Skip the header line
			if 'THERMO' in line:
				continue
			#Skip the end line
			if 'END' in line:
				continue
			#Skip comment lines
			if line[0] == '!':
				continue
			#Skip header temperatures
			if _is_temperature_header(line):
				continue

			'''
			Parse lines
			'''
			line_num = _read_line_num(line)
			if line_num == 1:
				nasa_data = _read_line1(line)
			elif line_num == 2:
				nasa_data = _read_line2(line, nasa_data)
			elif line_num == 3:
				nasa_data = _read_line3(line, nasa_data)
			elif line_num == 4:
				nasa_data = _read_line4(line, nasa_data)
				species.append(Nasa(**nasa_data))
			else:
				raise IOError('Invalid line number, {}, in thermdat file: {}'.format(line_num, filename))
	return species

def _get_fields(line, delimiter = ' ', remove_fields = ['', '\n']):
	"""Gets the fields from a line delimited by delimiter and without entries in remove_fields

	Parameters
	----------
		line : str
			Line to find the fields
		delimiter : str
			Text separating fields in string
		remove_fields : list of str
			Fields to delete
	Returns
	-------
		fields : list of str
			Fields of the line
	"""
	for remove_field in remove_fields:
		line = line.replace(remove_field, '')
	all_fields = line.split(delimiter)
	fields = []
	for field in all_fields:
		if all([field != remove_field for remove_field in remove_fields]):
			#Add field if it does not match any of the remove_fields
			fields.append(field)
	return fields

def _is_temperature_header(line):
	"""Determines if the line if the temperature header by seeing if the line only contains three numbers.

	Parameters
	----------
		line : str
			Line to test
	Returns
	-------
		temperature_header : bool
			True if the line is the temperature header. False otherwise.
	"""

	fields = _get_fields(line)
	n_num = 0
	for field in fields:
		#See if the field is a float
		try:
			float(field)
		except ValueError:
			#Contains text field
			return False
		else:
			n_num += 1
	#Temperature header contains 3 floats
	if n_num == 3:
		return True
	else:
		return False

def _read_line_num(line):
	"""Reads the line number. Assumes the line number is the last character

	Parameters
	----------
		line : str
			Line to be read
	Returns
	-------
		line_num : int
			Line number.
	"""

	fields = _get_fields(line)
	return int(fields[-1])

def _read_line1(line):
	"""Reads the first line of a thermdat specie

	Parameters
	----------
		line : str
			Line 1 of thermdat specie
	Returns
	-------
		nasa_data : dict
			Nasa input fields
	"""
	nasa_data = {}
	ref_pos = 24
	ref_offset = 5
	max_elements = 4
	phase_pos = 44

	#Store the name
	blank_pos = line.find(' ')
	nasa_data['name'] = line[:blank_pos]

	#Store the notes if any
	notes = line[blank_pos:ref_pos].strip()
	if len(notes) > 0:
		nasa_data['notes'] = notes

	#Store the elements
	nasa_data['elements'] = {}
	for i in range(max_elements):
		blank_pos = line.find(' ', ref_pos)
		#All the elements have been assigned
		if blank_pos == ref_pos:
			break

		element = line[ref_pos:blank_pos]
		ref_pos += ref_offset
		coeff = int(line[blank_pos:ref_pos])

		nasa_data['elements'][element] = coeff

	#Store the phase
	nasa_data['phase'] = line[phase_pos]

	#Store the temperatures
	fields = _get_fields(line[phase_pos+1:])
	nasa_data['T_low'] = float(fields[0])
	nasa_data['T_high'] = float(fields[1])
	nasa_data['T_mid'] = float(fields[2])
	return nasa_data

def _read_line2(line, nasa_data):
	"""Reads the second line of a thermdat specie

	Parameters
	----------
		line : str
			Line 2 of thermdat specie
		nasa_data : dict
			Pre-filled Nasa input fields
	Returns
		nasa_data : dict
			Nasa input fields
	"""
	#Locations to find a values
	positions = [0, 15, 30, 45, 60]
	offset = 15

	nasa_data['a_high'] = np.zeros(7)

	for i, position in enumerate(positions):
		nasa_data['a_high'][i] = float(line[position:position+offset])
	return nasa_data

def _read_line3(line, nasa_data):
	"""Reads the third line of a thermdat specie

	Parameters
	----------
		line : str
			Line 3 of thermdat specie
		nasa_data : dict
			Pre-filled Nasa input fields
	Returns
	-------
		nasa_data : dict
			Nasa input fields
	"""
	#Locations to find a values
	positions = [0, 15, 30, 45, 60]
	offset = 15

	nasa_data['a_low'] = np.zeros(7)

	j = 5 #Counter for a_high
	k = 0 #Counter for a_low
	for i, position in enumerate(positions):
		if i < 2:
			nasa_data['a_high'][j] = float(line[position:position+offset])
			j += 1
		else:
			nasa_data['a_low'][k] = float(line[position:position+offset])
			k += 1
	return nasa_data

def _read_line4(line, nasa_data):
	"""Reads the third line of a thermdat specie

	Parameters
	----------
		line : str
			Line 3 of thermdat specie
		nasa_data : dict
			Pre-filled Nasa input fields
	Returns
	-------
		nasa_data : dict
			Nasa input fields
	"""
	#Locations to find a values
	positions = [0, 15, 30, 45]
	offset = 15

	j = 3
	for position in positions:
		nasa_data['a_low'][j] = float(line[position:position+offset])
		j += 1
	return nasa_data

def write_thermdat(filename, nasa_species, write_date=True, newline='\n'):
	"""Writes thermdats in the Chemkin format

	Parameters
	----------
		filename : str
			Output file name
		nasa_species : list of ``PyMuTT.models.empirical.nasa.Nasa``
		write_date : bool, optional
			Whether or not the date should be written. If False, writes the first 8 characters of ``notes`` attribute. Defaults to True			
		newline : str, optional
			Newline character to use. Default is the Unix convention (\\n)
	"""
	with open(filename, 'w', newline=newline) as f_ptr:
		f_ptr.write('THERMO ALL\n       100       500      1500\n')
		
		float_string = '%.8E'
		for nasa_specie in nasa_species:
			_write_line1(f_ptr, nasa_specie, write_date)
			_write_line2(f_ptr, nasa_specie, float_string)
			_write_line3(f_ptr, nasa_specie, float_string)
			_write_line4(f_ptr, nasa_specie, float_string)
		f_ptr.write('END')

def _write_line1(thermdat_file, nasa_specie, write_date=True):
	"""Writes the first line of the thermdat file, which contains information on the composition, phase, and temperature ranges

	Parameters
	----------
		thermdat_file : file object
			Thermdat file that is being written to
		nasa_specie : ``PyMuTT.models.empirical.thermdat.Thermdat``
			Nasa specie to take information from
		write_date : bool, optional
			Whether or not the date should be written. If False, writes the first 8 characters of ``notes`` attribute. Defaults to True
	"""
	element_pos = [
		24, #Element 1
		28, #Element 1#
		29, #Element 2
		33, #Element 2#
		34, #Element 3
		38, #Element 3#
		39, #Element 4
		43] #Element 4#
	temperature_pos = [
		44, # Phase
		45, # T_low
		55, # T_high
		65, # T_mid
		79] # Line num

	#Adjusts the position based on the number of elements
	line1_pos = [16]
	for element, val in nasa_specie.elements.items():
		if val > 0.:
			two_digit = len(str(val)) - 1
			line1_pos.append(element_pos.pop(0))
			line1_pos.append(element_pos.pop(0) - two_digit)
	line1_pos.extend(temperature_pos)

	#Creating a list of the text to insert
	if write_date:
		now = datetime.now()
		notes = now.strftime('%Y%m%d')
	elif (nasa_specie.notes is None) or (nasa_specie.notes == ''):
		notes = ''
	else:
		notes = nasa_specie.notes[:8]

	line1_fields = [nasa_specie.name, notes]
	for element, val in nasa_specie.elements.items():
		if val > 0.:
			line1_fields.extend([element, '%d' % val])
	line1_fields.extend([nasa_specie.phase, '%.1f' % nasa_specie.T_low, '%.1f' % nasa_specie.T_high, '%.1f' % nasa_specie.T_mid])

	#Write the content with appropriate spacing
	line = ''
	for pos, field in zip(line1_pos, line1_fields):
		line += field
		line = _insert_space(pos, line)
	line += '1\n'
	thermdat_file.write(line)

def _write_line2(thermdat_file, nasa_specie, float_string):
	"""Writes the second line of the thermdat file

	Parameters
	----------
		thermdat_file : file object
			Thermdat file that is being written to
		nasa_specie : ``PyMuTT.models.empirical.thermdat.Thermdat``
			Nasa specie to take information from
		float_string : str
			float format
	"""
	line = ''
	for i in range(5):
		a = nasa_specie.a_high[i]
		if a >= 0:
			line += ' '
		line += float_string % a
	line += '    2\n'
	thermdat_file.write(line)

def _write_line3(thermdat_file, nasa_specie, float_string):
	"""Writes the third line of the thermdat file

	Parameters
	----------
		thermdat_file : file object
			Thermdat file that is being written to
		nasa_specie : ``PyMuTT.models.empirical.thermdat.Thermdat``
			Nasa specie to take information from
		float_string : str
			float format
	"""
	line = ''
	for i in range(5):
		if i < 2:
			a = nasa_specie.a_high[i+5]
		else:
			a = nasa_specie.a_low[i-2]
		if a >= 0:
			line += ' '
		line += float_string % a
	line += '    3\n'
	thermdat_file.write(line)

def _write_line4(thermdat_file, nasa_specie, float_string):
	"""Writes the fourth line of the thermdat file

	Parameters
	----------
		thermdat_file : file object
			Thermdat file that is being written to
		nasa_specie : ``PyMuTT.models.empirical.thermdat.Thermdat``
			Nasa specie to take information from
		float_string : str
			float format
	"""
	line = ''
	for i in range(3,7):
		a = nasa_specie.a_low[i]
		if a >= 0:
			line += ' '
		line += float_string % a
	line += '                   4\n'
	thermdat_file.write(line)

def _insert_space(end_index, string):
	"""Inserts the number of spaces required given the string and the position of the next non-blank field.

	Parameters:
		end_index : int
			Expected string length
		string : str
			String to add spaces to
	Returns:
		string_with_blanks : str
			String with spaces padded on the end
	"""
	string += ' ' * (end_index - len(string))
	return string

def thermdat_to_json(filename, thermdats):
	raise(NotImplementedError)

def json_to_thermdats(filename):
	raise(NotImplementedError)