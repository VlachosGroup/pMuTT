# -*- coding: utf-8 -*-

import numpy as np

from pmutt import _ModelBase
from pmutt import constants as c
from pmutt import get_molecular_weight
from pmutt.io.json import remove_class


class FreeTrans(_ModelBase):
    """Translational mode using ideal gas assumption. Equations sourced from:
    
    * Sandler, S. I. An Introduction to Applied Statistical Thermodynamics;
      John Wiley & Sons, 2010.

    Attributes
    ----------
        n_degrees : int, optional
            Number of degrees of freedom. Default is 3
        molecular_weight : float
            Molecular weight of the molecule in g/mol.
        atoms : `ase.Atoms`_ object, optional
            An atoms object can be used to calculate molecular weight. Not
            stored by FreeTrans

    .. _`ase.Atoms`: https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms
    """
    def __init__(self, n_degrees=3, molecular_weight=None, atoms=None):
        self.n_degrees = n_degrees
        if molecular_weight is None and atoms is not None:
            self.molecular_weight = get_molecular_weight(
                atoms.get_chemical_formula(mode='hill'))
        else:
            self.molecular_weight = molecular_weight

    def get_V(self, T, P):
        """Calculates the molar volume of an ideal gas at T and P

        :math:`V_m=\\frac{RT}{P}`

        Parameters
        ----------
            T : float
                Temperature in K
            P : float
                Pressure in bar
        Returns
        -------
            V : float
                Molar volume in m3
        """
        return T * c.R('J/mol/K') / (P *
                                     c.convert_unit(initial='bar', final='Pa'))

    def get_q(self, T, P=c.P0('bar')):
        """Calculates the partition function

        :math:`q_{trans} = \\bigg(\\frac{2\\pi \\sum_{i}^{atoms}m_ikT}{h^2}
        \\bigg)^\\frac {n_{degrees}} {2}V`

        Parameters
        ----------
            T : float
                Temperature in K
            P : float, optional
                Pressure (bar) or pressure-like quantity.
                Default is atmospheric pressure
        Returns
        -------
            q_trans : float
                Translational partition function
        """
        molecular_volume = self.get_V(T=T, P=P)/c.Na
        unit_mass = self.molecular_weight *\
            c.convert_unit(initial='g', final='kg')/c.Na
        return molecular_volume*(2.*np.pi*c.kb('J/K')*T*unit_mass/(c.h('J s')**2.0)) \
            ** (float(self.n_degrees)/2.)

    def get_CvoR(self):
        """Calculates the dimensionless heat capacity at constant volume

        :math:`\\frac{Cv^{trans}}{R}=\\frac{n_{degrees}}{2}`

        Returns
        -------
            CvoR_trans : float
                Translational dimensionless heat capacity at constant V
        """
        return float(self.n_degrees) / 2.

    def get_CpoR(self):
        """Calculates the dimensionless heat capacity at constant pressure

        :math:`\\frac{Cp^{trans}}{R}=\\frac{Cv^{trans}}{R} + 1`

        Returns
        -------
            CpoR_trans : float
                Translational dimensionless heat capacity at constant P
        """
        return self.get_CvoR() + 1.

    def get_UoRT(self):
        """Calculates the dimensionless internal energy

        :math:`\\frac{U^{trans}}{RT}=\\frac{n_{degrees}}{2}`

        Returns
        -------
            UoRT_trans : float
                Translational internal energy
        """
        return float(self.n_degrees) / 2.

    def get_HoRT(self):
        """Calculates the dimensionless enthalpy

        :math:`\\frac{H^{trans}}{RT}=\\frac{U^{trans}}{RT} + 1`

        Returns
        -------
            HoRT_trans : float
                Translational enthalpy
        """
        return self.get_UoRT() + 1.

    def get_SoR(self, T, P=c.P0('bar')):
        """Calculates the dimensionless entropy

        :math:`\\frac{S^{trans}}{R}=1+\\frac{n_{degrees}}{2}+\\log\\bigg(\\big(
        \\frac{2\\pi mk_bT}{h^2})^\\frac{n_{degrees}}{2}\\frac{RT}{PN_a}\\bigg)`

        Parameters
        ----------
            T : float
                Temperature in K
            P : float, optional
                Pressure (bar) or pressure-like quantity.
                Default is atmospheric pressure

        Returns
        -------
            SoR_trans : float
                Translational dimensionless entropy
        """
        V = self.get_V(T=T, P=P)
        unit_mass = self.molecular_weight *\
            c.convert_unit(initial='g', final='kg')/c.Na
        return 1. + float(self.n_degrees)/2. \
            + np.log((2.*np.pi*unit_mass*c.kb('J/K')*T/c.h('J s')**2)
                     ** (float(self.n_degrees)/2.)*V/c.Na)

    def get_FoRT(self, T, P=c.P0('bar')):
        """Calculates the dimensionless Helmholtz energy

        :math:`\\frac{A^{trans}}{RT}=\\frac{U^{trans}}{RT}-\\frac{S^{trans}}{R}`

        Parameters
        ----------
            T : float
                Temperature in K
            P : float, optional
                Pressure (bar) or pressure-like quantity.
                Default is atmospheric pressure
        Returns
        -------
            FoRT_trans : float
                Translational dimensionless Helmholtz energy
        """
        return self.get_UoRT() - self.get_SoR(T=T, P=P)

    def get_GoRT(self, T, P=c.P0('bar')):
        """Calculates the dimensionless Gibbs energy

        :math:`\\frac{G^{trans}}{RT}=\\frac{H^{trans}}{RT}-\\frac{S^{trans}}{R}`       

        Parameters
        ----------
            T : float
                Temperature in K
            P : float, optional
                Pressure (bar) or pressure-like quantity.
                Default is atmospheric pressure
        Returns
        -------
            GoR_trans : float
                Translational dimensionless Gibbs energy
        """
        return self.get_HoRT() - self.get_SoR(T=T, P=P)

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        return {
            'class': str(self.__class__),
            'n_degrees': self.n_degrees,
            'molecular_weight': self.molecular_weight
        }
