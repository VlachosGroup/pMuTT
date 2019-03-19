# -*- coding: utf-8 -*-
"""
pMuTT.empirical.references

Operations related to referencing DFT energy to enthalpies of experimental
references.
"""

from warnings import warn
import numpy as np
from pMuTT.empirical import EmpiricalBase
from pMuTT.io.json import json_to_pMuTT, remove_class


class Reference(EmpiricalBase):
    """Single reference specie used to adjust DFT energies to
    experimental data.

    Attributes
    ----------
        T_ref : float
            Temperature reference in K
        HoRT_ref : float
            Dimensionless enthalpy corresponding to T_ref
    """

    def __init__(self, T_ref, HoRT_ref, **kwargs):
        super().__init__(**kwargs)
        self.T_ref = T_ref
        self.HoRT_ref = HoRT_ref

    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = super().to_dict()
        obj_dict['T_ref'] = self.T_ref
        obj_dict['HoRT_ref'] = self.HoRT_ref
        return obj_dict

    @classmethod
    def from_dict(cls, json_obj):
        """Recreate an object from the JSON representation.

        Parameters
        ----------
            json_obj : dict
                JSON representation
        Returns
        -------
            Reference : Reference object
        """
        json_obj = remove_class(json_obj)
        # Reconstruct statmech model
        return cls(**json_obj)


class References:
    """Holds reference species to adjust DFT energies to experimental data.

    Attributes
    ----------
        offset : dict
            Dimensionless enthalpy offset for each descriptor
        references : list of :class:`~pMuTT.empirical.references.Reference`,
        optional
            Reference species. Each member of the list should have the
            attributes ``T_ref`` and ``HoRT_ref``
        T_ref : float
            Reference temperature in K
    """
    def __init__(self, offset=None, references=None, descriptor='elements'):
        self.offset = offset
        self.references = references
        self.descriptor = descriptor
        # If offset not specified but references is specified
        if self.offset is None and self.references is not None:
            self.fit_HoRT_offset()

    def __iter__(self):
        """Iterates over references attribute

        Yields
        ------
            reference : :class:`~pMuTT.empirical.references.Reference`
        """
        for reference in self.references:
            yield reference

    def __len__(self):
        return len(self.references)

    def __setitem__(self, index, reference):
        self.references[index] = reference

    def __getitem__(self, index):
        return self.references[index]

    def append(self, obj):
        self.references.append(obj)

    def extend(self, seq):
        self.references.extend(seq)

    def insert(self, obj):
        self.references.insert(obj)

    def pop(self, obj=-1):
        self.references.pop(obj)

    def remove(self, obj):
        self.references.remove(obj)

    def index(self, name):
        for i, reference in enumerate(self):
            if name == reference.name:
                return i
        else:
            return None

    def clear_offset(self):
        """Removes all entries from descriptor offset dictionary."""
        self.offset.clear()

    def get_descriptors(self):
        """Returns the descriptors in references.

        Returns
        -------
            descriptors : tuple
                Unique descriptors in reference species
        """
        unique_descriptors = []
        for reference in self.references:
            for desc in getattr(reference, self.descriptor).keys():
                if desc not in unique_descriptors:
                    unique_descriptors.append(desc)
        return tuple(sorted(unique_descriptors))

    def get_descriptors_matrix(self):
        """Creates the descriptors matrix required for calculating the offset.
        The descriptors are sorted in alphabetical order.

        Returns
        -------
            descriptor matrix : (M,N) `numpy.ndarray`_
                Rows correspond to reference species. Columns correspond to
                descriptors
        """
        descriptors = self.get_descriptors()
        descriptors_mat = np.zeros((len(self.references), len(descriptors)))
        for i, reference in enumerate(self):
            for j, descriptor_name in enumerate(descriptors):
                try:
                    descriptors_mat[i, j] = \
                        getattr(reference,
                                self.descriptor)[descriptor_name]
                except KeyError:
                    # If descriptor not in dictionary
                    descriptors_mat[i, j] = 0.
        return descriptors_mat

    def fit_HoRT_offset(self):
        """Calculate the descriptoral offset between DFT and formation energies
        using reference species."""
        descriptors = self.get_descriptors()
        descriptors_mat = self.get_descriptors_matrix()

        T_refs = np.array([reference.T_ref for reference in self])
        # If any of the T_ref values are not close to the others.
        if any([not np.isclose(T_refs[0], T_ref) for T_ref in T_refs]):
            warn('All the reference temperatures are not the same. May cause '
                 'error in referencing. Using mean temperature.')
        self.T_ref = np.mean(T_refs)

        HoRT_ref_dft = np.array(
                [reference.statmech_model.get_HoRT(T=reference.T_ref)
                    for reference in self])
        HoRT_ref_exp = np.array([reference.HoRT_ref for reference in self])
        # Offset between the DFT energies and experimentalvalues
        # for reference species
        ref_offset = HoRT_ref_dft - HoRT_ref_exp
        # Offset between the DFT energies and experimental values for
        # each descriptor
        offset = np.linalg.lstsq(descriptors_mat, ref_offset,
                                 rcond=None)[0]
        # Convert offset to a dictionary
        self.offset = {descriptor: val for descriptor, val
                       in zip(descriptors, offset)}

    def get_HoRT_offset(self, descriptors, T=None):
        """Returns the offset due to the descriptor composition of a specie.
        The offset is defined as follows:

        :math:`HoRT_{exp} = HoRT_{dft} + offset`

        Parameters
        ----------
            descriptors : dict
                Dictionary where the keys are decriptors and the values are the
                number of each descriptor in a formula unit
            Ts : float or (N,) numpy.ndarray_
                Temperatures in K. If not specified, adjusts using ``T_ref``
        Returns
        -------
            HoRT_offset : float
                Offset to add to potentialenergy (in eV) to adjust to
                References

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        HoRT_offset = 0.
        for descriptor, coefficient in descriptors.items():
            try:
                HoRT_offset -= self.offset[descriptor]*coefficient
            except KeyError:
                warn('References does not have offset value for the '
                     'descriptor: {}.'.format(descriptor), RuntimeWarning)
        if T is None:
            return HoRT_offset
        else:
            # Adjust for the temperature
            return HoRT_offset * self.T_ref/T

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = {'class': str(self.__class__),
                    'offset': self.offset,
                    'descriptor': self.descriptor}
        try:
            obj_dict['references'] = [ref.to_dict() for ref in self.references]
        except (AttributeError, TypeError):
            obj_dict['references'] = self.references
        return obj_dict

    @classmethod
    def from_dict(cls, json_obj):
        """Recreate an object from the JSON representation.

        Parameters
        ----------
            json_obj : dict
                JSON representation
        Returns
        -------
            References : References object
        """
        json_obj = remove_class(json_obj)
        json_obj['offset'] = np.array(json_obj['offset'])
        json_obj['references'] = [json_to_pMuTT(ref_dict)
                                  for ref_dict in json_obj['references']]
        return cls(**json_obj)
