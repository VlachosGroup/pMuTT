# -*- coding: utf-8 -*-
"""
PyMuTT
Vlachos group code to read, write, and generate thermochemical data.
Created on Fri Jul 7 12:40:00 2018
"""

import re
import inspect
from PyMuTT import constants as c

def _get_expected_arguments(fn):
	"""Returns the arguments expected by a function. Useful for determining
	where to assign **kwargs parameters.
	
	Args
		fn (function or class): Function or class you would like to find the expected arguments.
	Returns
		tuple of str: Expected arguments. If a class is specified, returns the expected arguments of __init__
	"""

	#If class passed, use __init__ to find expected arguments
	if inspect.isclass(fn):
		fn = fn.__init__

	fn_code = fn.__code__
	arg_count = fn_code.co_argcount
	args = fn_code.co_varnames[:arg_count]
	return args

def _pass_expected_arguments(fn, **kwargs):
	"""Finds expected values from a function or class and passes the appropriate arguments.

	Args
		fn (Function or Class): Function or class you would like to find the expected arguments.
		**kwargs: Keyword arguments that contain parameters to pass to fn
	Returns
		Output of fn that has been fed the expected arguments.
	"""
	expected_args = _get_expected_arguments(fn)
	expected_arg_val = {}
	for arg in expected_args:
		try:
			expected_arg_val[arg] = kwargs[arg]
		except KeyError:
			continue
	return fn(**expected_arg_val)

def parse_formula(formula):
	"""Parses chemical formula into its elements and returns it as a dictionary.

	Args	
		formula (str): Chemical formula e.g. Al2O3
	Returns
		elements (dict): ELement composition of formula e.g. {'Al': 2, 'O': 3}
	"""
	elements_tuples = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
	elements = {}
	for (element, coefficient) in elements_tuples:
		if coefficient == '':
			elements[element] = 1
		else:
			elements[element] = int(coefficient)
	return elements

def get_molecular_weight(elements):
	"""Molecular mass (in g/mol) given the elemental composition.
	Data taken from: https://en.wikipedia.org/wiki/Standard_atomic_weight

	Parameters
		elements (:obj: `dict` , :obj: `str`): Elemental composition of species.
			If a dictionary is passed, the keys are the element symbol, atomic number, 
			or element name and the value is the stoichiometric coefficient.
			If a string is passed, the formula will be guessed using PyMuTT.parse_formula

	Returns
		molecular_weight (float): Molecular weight as float in kg/mol
	"""
	if isinstance(elements, str):
		elements = parse_formula(elements)

	molecular_weight = 0.
	for element, coefficient in elements.items():
		molecular_weight += c.atomic_weight[element] * coefficient
	return molecular_weight