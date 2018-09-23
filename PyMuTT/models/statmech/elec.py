# -*- coding: utf-8 -*-

import numpy as np
from PyMuTT import constants as c
from PyMuTT.io_.jsonio import remove_class

class IdealElec:
    """Electronic modes using the ideal gas assumption. Equations found in
    Sandler, S. I. An Introduction to Applied Statistical Thermodynamics; 
    John Wiley & Sons, 2010.
    
    Attributes
    ----------
        potentialenergy : float, optional
            Potential energy in eV. Default is 0
        spin : float, optional
            The total electron spin. Default is 0
            
            - 0 for molecules in which all electrons are paired
            - 0.5 for a free radical with a single unpaired electron
            - 1.0 for a triplet with two unpaired electrons, such as O :sub:`2`
    """

    def __init__(self, potentialenergy=0., spin=0.):
        self.potentialenergy = potentialenergy
        self.spin = spin
        self._degeneracy = 2.*self.spin + 1.

    def get_q(self, T, ignore_q_elec=False):
        """Calculates the partition function

        :math:`q^{elec}=\\omega_i \\exp\\bigg(-\\frac{E}{RT}\\bigg)`

        Parameters
        ----------
            T : float
                Temperature in K
            ignore_q_elec : bool, optional
                Ignore the contribution of electronic mode to partition function
                . Often necessary since DFT's value for potentialenergy is
                very negative causing q_elec to go to infinity.
        Returns
        -------
            q_elec : float
                Electronic partition function
        """
        if ignore_q_elec:
            return 1.
        else:
            return self._degeneracy*np.exp(-self.get_UoRT(T=T))

    def get_CvoR(self):
        """Calculates the dimensionless heat capacity at constant volume

        :math:`\\frac{C_V^{elec}}{R}=0`

        Returns
        -------
            CvoR_elec : float
                electronic dimensionless heat capacity at constant volume
        """
        return 0.

    def get_CpoR(self):
        """Calculates the dimensionless heat capacity at constant pressure

        :math:`\\frac{C_V^{elec}}{R}=\\frac{C_P^{elec}}{R}`

        Returns
        -------
            CpoR_elec : float
                Electronic dimensionless heat capacity at constant pressure
        """
        return self.get_CvoR()
    
    def get_UoRT(self, T):
        """Calculates the imensionless internal energy

        :math:`\\frac{U^{elec}}{RT}=\\frac{E}{RT}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            UoRT_elec : float
                Electronic dimensionless internal energy
        """
        return self.potentialenergy/c.kb('eV/K')/T

    def get_HoRT(self, T):
        """Calculates the dimensionless enthalpy

        :math:`\\frac{H^{elec}}{RT}=\\frac{U^{elec}}{RT}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            HoRT_elec : float
                Electronic dimensionless enthalpy
        """
        return self.get_UoRT(T=T)

    def get_SoR(self):
        """Calculates the dimensionless entropy

        :math:`\\frac{S^{elec}}{R}=\\log \\omega`

        Returns
        -------
            SoR_elec : float
                Electronic dimensionless entropy
        """
        return np.log(self._degeneracy)

    def get_AoRT(self, T):
        """Calculates the dimensionless Helmholtz energy

        :math:`\\frac{A^{elec}}{RT}=\\frac{U^{elec}}{RT}-\\frac{S^{elec}}{R}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            AoRT_elec : float
                Electronic dimensionless Helmholtz energy
        """
        return self.get_UoRT(T=T) - self.get_SoR()

    def get_GoRT(self, T):
        """Calculates the dimensionless Gibbs energy

        :math:`\\frac{G^{elec}}{RT}=\\frac{H^{elec}}{RT}-\\frac{S^{elec}}{R}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            GoRT_elec : float
                Electronic dimensionless Gibbs energy
        """
        return self.get_HoRT(T=T) - self.get_SoR()

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes
        
        Returns
        -------
            obj_dict : dict
        """
        return {'class': str(self.__class__),
                'potentialenergy': self.potentialenergy, 
                'spin': self.spin}

    @classmethod
    def from_dict(cls, json_obj):
        """Recreate an object from the JSON representation.

        Parameters
        ----------
            json_obj : dict
                JSON representation
        Returns
        -------
            IdealElec : IdealElec object
        """
        json_obj = remove_class(json_obj)
        return cls(**json_obj)