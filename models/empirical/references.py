# -*- coding: utf-8 -*-
"""
PyMuTT.references
Vlachos group code for reference species.
Created on Sun Jul 8 13:50:00 2018
"""
from warnings import warn
import numpy as np
from PyMuTT.models.empirical import BaseThermo

class References:
	"""
	Holds reference species to adjust DFT energies to experimental data.

	Attributes
		_references - list
			Reference species. Each member of the list should have the attributes:
				T_ref
				HoRT_ref
		HoRT_element_offset - dict
			Dimensionless enthalpy offset for each element
		T_ref - float
			Reference temperature in K
	"""
	def __init__(self, references):
		self._references = references
		self.HoRT_element_offset = {}
		self.fit_HoRT_offset()

	def __iter__(self):
		"""
		Iterates over references attribute
		"""
		for reference in self._references:
			yield reference

	def __len__(self):
		"""
		Returns length of references attribute
		"""
		return len(self._references)

	def __setitem__(self, index, reference):
		self._references[index] = reference

	def __getitem__(self, index):
		return self._references[index]

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
		self.HoRT_element_offset.clear()

	def get_elements(self):
		"""
		Returns the elements in references.

		Returns
			tuple
				Unique elements in reference species
		"""
		unique_elements = []
		for reference in self:
			for element in reference.elements.keys():
				if element not in unique_elements:
					unique_elements.append(element)
		return tuple(sorted(unique_elements))

	def get_elements_matrix(self):
		"""
		Creates the elements matrix required for calculating the offset. The elements are sorted
		in alphabetical order.

		Returns
			(M,N) ndarray
				Rows correspond to reference species. Columns correspond to elements
		"""
		elements = self.get_elements()
		elements_mat = np.zeros((len(self), len(elements)))
		for i, reference in enumerate(self):
			for j, element in enumerate(elements):
				try:
					elements_mat[i, j] = reference.elements[element]
				except KeyError:
					#If element not in dictionary
					elements_mat[i, j] = 0.
		return elements_mat

	def fit_HoRT_offset(self):
		"""
		Calculate the elemental offset between DFT and formation energies using reference species.
		"""
		elements = self.get_elements()
		elements_mat = self.get_elements_matrix()

		T_refs = np.array([reference.T_ref for reference in self])
		#If any of the T_ref values are not close to the others.
		if any([not np.isclose(T_refs[0], T_ref) for T_ref in T_refs]):
			warn('All the reference temperatures are not the same. May cause error in referencing. Using mean temperature.')
		self.T_ref = np.mean(T_refs)

		HoRT_ref_dft = np.array([reference.thermo_model.get_HoRT(Ts=reference.T_ref) for reference in self])
		HoRT_ref_exp = np.array([reference.HoRT_ref for reference in self])
		ref_offset = HoRT_ref_dft - HoRT_ref_exp #Offset between the DFT energies and experimental values for reference species
		HoRT_element_offset = np.linalg.lstsq(elements_mat, ref_offset, rcond=None)[0] #Offset between the DFT energies and experimental values for each element
		#Convert HoRT_element_offset to a dictionary
		self.HoRT_element_offset = {element: offset for element, offset in zip(elements, HoRT_element_offset)}

	def get_HoRT_offset(self, elements, Ts=None):
		"""
		Returns the offset due to the element composition of a specie. The offset is defined as follows:
		HoRT_exp = HoRT_dft + offset
		calc_offset must be run before meaningful offset values can be returned.

		Parameters
			elements - dict
				Dictionary where the keys are elements and the values are
				the number of each element in a formula unit
			Ts - float or (N,) ndarray
				Temperatures in K. If not specified, adjusts using T_ref
		Returns
			float
				Offset to add to potentialenergy (in eV) to adjust to References
		"""
		HoRT_offset = 0.
		for element, coefficient in elements.items():
			try:
				HoRT_offset -= self.HoRT_element_offset[element] * coefficient
			except KeyError:
				warn('References does not have offset value for the element: {}.'.format(element), RuntimeWarning)
		if Ts is None:
			return HoRT_offset
		else:
			#Adjust for the temperature
			return HoRT_offset * self.T_ref / Ts