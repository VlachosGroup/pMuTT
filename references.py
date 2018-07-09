# -*- coding: utf-8 -*-
"""
Thermochemistry.references
Vlachos group code for reference species.
Created on Sun Jul 8 13:50:00 2018
"""

import collections
from Thermochemistry import BaseThermo

class References:
	"""
	Holds reference species to adjust DFT energies to experimental data.

	Attributes
		_references - list
			Reference species. Each member of the list should have the attributes:
				T_ref
				HoRT_ref
		element_offset - dict
			Enthalpy offset for each element
	"""
	def __init__(self, references):
		self._references = references
		self._element_offset = {}

	def __iter__(self):
		"""
		Iterates over references attribute
		"""
		for reference in self._references:
			yield reference

	def append(self, obj):
		"""
		Operates the same as list.append using references attribute
		"""
		self._references.append(obj)

	def extend(self, seq):
		"""
		Operates the same as list.extend using references attribute
		"""
		self._references.extend(seq)

	def insert(self, obj):
		"""
		Operates the same as list.insert using references attribute
		"""
		self._references.insert(obj)

	def pop(self, obj=-1):
		"""
		Operates the same as list.pop using references attribute
		"""
		self._references.pop(obj)

	def remove(self, obj):
		"""
		Operates the same as list.remove using references attribute
		"""
		self._references.remove(obj)

	def index(self, obj):
		"""
		Operates the same as list.index using references attribute
		"""
        for i, reference in enumerate(self):
            if obj == reference.name:
                return i
        else:
            return None

    def clear_offset(self):
    	"""
    	Removes all entries from element offset dictionary.
    	"""
    	self._element_offset.clear()

    def calc_offset(self):
    	pass

    def get_specie_offset(self, elements):
    	"""
    	Returns the offset due to the element composition of a specie.
    	calc_offset must be run before meaningful offset values can be returned.

    	Parameters
    		elements - dict
    			Dictionary where the keys are elements and the values are
    			the number of each element in a formula unit
    	Returns
    		float
    			Offset to add to potentialenergy (in eV) to adjust to References
		Raises
			KeyError
				If any element in elements is not found in _element_offset
    	"""
    	offset = 0.
    	for element, val in elements.items():
    		offset += self._element_offset[element] * val
    	return offset