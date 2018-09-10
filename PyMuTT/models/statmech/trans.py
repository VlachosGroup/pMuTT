# -*- coding: utf-8 -*-

import numpy as np
from PyMuTT import constants as c

class IdealTrans:
    """Translational mode using ideal gas assumption

    Attributes
    ----------
        n_degrees : int, optional
            Number of degrees of freedom. Default is 3
        molecular_weight : float
            Molecular weight of the molecule in g/mol.
    """

    def __init__(self, molecular_weight, n_degrees=3):
        self.molecular_weight = molecular_weight
        self.n_degrees = n_degrees

    def get_q(self, T, V=c.V0):
        """Calculates the partition function

        :math:`q_{trans} = \\bigg(\\frac{2\\pi \\sum_{i}^{atoms}m_ikT}{h^2}
        \\bigg)^\\frac {n_{degrees}} {2}V`

        Parameters
        ----------
            T : float
                Temperature in K
            V : float, optional
                Volume or volume-like (if n_degrees < 3) quantity. 
                Default is molar volume at standard conditions

                +-----------+--------------+----------------+
                | n_degrees | Meaning of V | Expected units |
                +===========+==============+================+
                | 1         | Length       | m              |
                +-----------+--------------+----------------+
                | 2         | Area         | m2             |
                +-----------+--------------+----------------+
                | 3         | Volume       | m3             |
                +-----------+--------------+----------------+

        Returns
        -------
            q_trans : float
                Translational partition function
        """
        unit_mass = self.molecular_weight*c.convert_unit(from_='g', to='kg')/c.Na
        return V*(2*np.pi*c.kb('J/K')*T*unit_mass/c.h('J s')**2) \
            **(float(self.n_degrees)/2.)

    def get_CvoR(self):
        """Calculates the dimensionless heat capacity at constant volume

        :math:`\\frac{Cv^{trans}}{R}=\\frac{n_{degrees}}{2}`
        
        Returns
        -------
            CvoR_trans : float
                Translational dimensionless heat capacity at constant V
        """
        return float(self.n_degrees)/2.

    def get_CpoR(self):
        """Calculates the dimensionless heat capacity at constant pressure

        :math:`\\frac{Cv^{trans}}{R}=\\frac{Cv^{trans}}{R} + 1`
        
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
        return float(self.n_degrees)/2.

    def get_HoRT(self):
        """Calculates the dimensionless enthalpy

        :math:`\\frac{H^{trans}}{RT}=\\frac{U^{trans}}{RT} + 1`

        Returns
        -------
            HoRT_trans : float
                Translational enthalpy
        """
        return self.get_UoRT() + 1.

    def get_SoR(self, T, V=c.V0):
        """Calculates the dimensionless entropy

        :math:`\\frac{S^{trans}}{R} = 1 + \\ln {\\frac {q_{trans}}{N_A}} 
        + \\frac{U}{RT}`

        Parameters
        ----------
            T : float
                Temperature in K
            V : float, optional
                Volume or volume-like (if n_degrees < 3) quantity. 
                Default is molar volume at standard conditions

                +-----------+--------------+----------------+
                | n_degrees | Meaning of V | Expected units |
                +===========+==============+================+
                | 1         | Length       | m              |
                +-----------+--------------+----------------+
                | 2         | Area         | m2             |
                +-----------+--------------+----------------+
                | 3         | Volume       | m3             |
                +-----------+--------------+----------------+

        Returns
        -------
            SoR_trans : float
                Translational dimensionless entropy 
        """
        return 1. + self.get_UoRT() + np.log(self.get_q(T=T, V=V)/c.Na)

    def get_AoRT(self, T, V=c.V0):
        """Calculates the dimensionless Helmholtz energy

        :math:`\\frac{A^{trans}}{RT}=\\frac{U^{trans}}{RT}-\\frac{S^{trans}}{R}`

        Parameters
        ----------
            T : float
                Temperature in K
            V : float, optional
                Volume or volume-like (if n_degrees < 3) quantity. 
                Default is molar volume at standard conditions

                +-----------+--------------+----------------+
                | n_degrees | Meaning of V | Expected units |
                +===========+==============+================+
                | 1         | Length       | m              |
                +-----------+--------------+----------------+
                | 2         | Area         | m2             |
                +-----------+--------------+----------------+
                | 3         | Volume       | m3             |
                +-----------+--------------+----------------+

        Returns
        -------
            AoRT_trans : float
                Translational dimensionless Helmholtz energy 
        """
        return self.get_UoRT()-self.get_SoR(T=T, V=V)

    def get_GoRT(self, T, V=c.V0):
        """Calculates the dimensionless Gibbs energy

        :math:`\\frac{G^{trans}}{RT}=\\frac{H^{trans}}{RT}-\\frac{S^{trans}}{R}`       
        
        Parameters
        ----------
            T : float
                Temperature in K
            V : float, optional
                Volume or volume-like (if n_degrees < 3) quantity. 
                Default is molar volume at standard conditions

                +-----------+--------------+----------------+
                | n_degrees | Meaning of V | Expected units |
                +===========+==============+================+
                | 1         | Length       | m              |
                +-----------+--------------+----------------+
                | 2         | Area         | m2             |
                +-----------+--------------+----------------+
                | 3         | Volume       | m3             |
                +-----------+--------------+----------------+

        Returns
        -------
            GoR_trans : float
                Translational dimensionless Gibbs energy 
        """
        return self.get_HoRT()-self.get_SoR(T=T, V=V)