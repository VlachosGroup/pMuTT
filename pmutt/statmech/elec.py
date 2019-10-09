# -*- coding: utf-8 -*-

import numpy as np

from pmutt import _ModelBase
from pmutt import constants as c
from pmutt.io.json import remove_class


class GroundStateElec(_ModelBase):
    """Electronic modes using the ideal gas assumption. Equations sourced from:
    
    * Sandler, S. I. An Introduction to Applied Statistical Thermodynamics;
      John Wiley & Sons, 2010.

    Attributes
    ----------
        potentialenergy : float, optional
            Potential energy in eV. If this value and `atoms` are not
            specified, defaults to 0
        spin : float, optional
            The total electron spin. Default is 0

            - 0 for molecules in which all electrons are paired
            - 0.5 for a free radical with a single unpaired electron
            - 1.0 for a triplet with two unpaired electrons, such as O :sub:`2`
        atoms : `ase.Atoms`_ object, optional
            If `potentialenergy` is not specified, the energy will be obtained
            using `atoms.get_potential_energy()`.
        D0 : float, optional
            Bond strength in eV. Preferentially used when calculating partition
            coefficient.

    .. _`ase.Atoms`: https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms
    """

    def __init__(self, potentialenergy=None, spin=0., atoms=None, D0=None):
        # If the potentialenergy was not specified, calculate from atoms object
        if potentialenergy is None:
            # If atoms was not specified, assign default value
            if atoms is None:
                potentialenergy = 0.
            else:
                # Try to read the potentialenergy from the atoms object. If it
                # fails, use the default value
                try:
                    potentialenergy = atoms.get_potential_energy()
                except (RuntimeError, ValueError):
                    potentialenergy = 0. 
        self.potentialenergy = potentialenergy
        self.D0 = D0
        self.spin = spin

    @property
    def spin(self):
        return self._spin
    
    @spin.setter
    def spin(self, val):
        self._spin = val
        self._degeneracy = 2.*val + 1.

    def get_q(self, T, ignore_q_elec=True):
        """Calculates the partition function

        :math:`q^{elec}=1 + \\omega_i \\exp\\bigg(-\\frac{E}{RT}\\bigg)`

        Parameters
        ----------
            T : float
                Temperature in K
            ignore_q_elec : bool, optional
                Ignore contribution of electronic mode to partition function
                . Often necessary since DFT's value for potentialenergy is
                very negative causing q_elec to go to infinity. Default is True
        Returns
        -------
            q_elec : float
                Electronic partition function
        """
        if ignore_q_elec:
            return 1.
        else:
            if self.D0 is not None:
                Epsilon = self.D0/c.kb('eV/K')/T
            else:
                Epsilon = self.get_UoRT(T=T)
            return self._degeneracy*(1 + np.exp(-Epsilon))

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
        return (self.potentialenergy)/c.kb('eV/K')/T

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

    def get_FoRT(self, T):
        """Calculates the dimensionless Helmholtz energy

        :math:`\\frac{A^{elec}}{RT}=\\frac{U^{elec}}{RT}-\\frac{S^{elec}}{R}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            FoRT_elec : float
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
