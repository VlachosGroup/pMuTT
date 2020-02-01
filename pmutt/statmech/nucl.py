# -*- coding: utf-8 -*-
from pmutt import _ModelBase
from pmutt.io.json import remove_class


class EmptyNucl(_ModelBase):
    """Nuclear modes. Assumes no change in any chemical reaction and hence
    does not affect thermodynamic quantities."""
    def __init__(self):
        pass

    def get_q(self):
        """Calculates the partition function

        Returns
        -------
            q_nucl : float
                Nuclear partition function
        """
        return 1.

    def get_CvoR(self):
        """Calculates the dimensionless heat capacity at constant volume

        Returns
        -------
            CvoR_nucl : float
                Nuclear dimensionless heat capacity at constant volume
        """
        return 0.

    def get_CpoR(self):
        """Calculates the dimensionless heat capacity at constant pressure

        Returns
        -------
            CpoR_nucl : float
                Nuclear dimensionless heat capacity at constant pressure
        """
        return 0.

    def get_UoRT(self):
        """Calculates the imensionless internal energy

        Returns
        -------
            UoRT_nucl : float
                Nuclear dimensionless internal energy
        """
        return 0.

    def get_HoRT(self):
        """Calculates the dimensionless enthalpy

        Returns
        -------
            HoRT_nucl : float
                Nuclear dimensionless enthalpy
        """
        return 0.

    def get_SoR(self):
        """Calculates the dimensionless entropy

        Returns
        -------
            SoR_nucl : float
                Nuclear dimensionless entropy
        """
        return 0.

    def get_FoRT(self):
        """Calculates the dimensionless Helmholtz energy

        Returns
        -------
            FoRT_nucl : float
                Nuclear dimensionless Helmholtz energy
        """
        return 0.

    def get_GoRT(self):
        """Calculates the dimensionless Gibbs energy

        Returns
        -------
            GoRT_nucl : float
                Nuclear dimensionless Gibbs energy
        """
        return 0.
