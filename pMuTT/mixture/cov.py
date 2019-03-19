# -*- coding: utf-8 -*-
"""
pMuTT.mixture.cov
Vlachos group code for coverage effects.
Created on Thurs Feb 7 10:15:00 2019
"""

import numpy as np
from pMuTT import constants as c
from pMuTT.io.jsonio import remove_class


class CovEffect:
    """Models enthalpic contribution to coverage effect between two species

    Attributes
    ----------
        name_i : str
            Name of specie's properties affected by coverage of specie_j
        name_j : str
            Name of specie affecting the properties of specie_i
        intervals : list (length N) of floats
            Intervals (in ML) to change slopes. The first element of the list
            should be 0 and the list should be sorted in ascending order to
            correctly specify a piecewise function with mole fraction domain
            [0, 1]
        slopes : list (length N) of float
            Slopes (in kcal/mol) to use between the intervals
    """

    def __init__(self, name_i, name_j, intervals, slopes):
        self.name_i = name_i
        self.name_j = name_j
        self.intervals = intervals
        self.slopes = slopes
        self._set_intercepts()

    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict

    def insert(self, interval, slope):
        """Inserts the a new interval and slope for the piecewise function

        Parameters
        ----------
            interval : float
                Interval to insert. It will automatically be inserted in the
                intervals list in the position so that piecewise function
                remains continuous
            slope : float
                Slope to insert. Inserted into the same position as interval
        """
        i = np.argmax(interval < np.array(self.intervals))
        self.intervals.insert(i, interval)
        self.slopes.insert(i, slope)
        self._set_intercepts()

    def pop(self, i):
        """Removes the interval and slope specified by an index

        Parameters
        ----------
            i : float
                Index to remove from ``self.intervals`` and ``self.slopes``.
                Value cannot be 0 since removing this index would not allow
                mole fractions to span from 0 to 1
        Raises
        ------
            ValueError
                Raised when i = 0
        """
        if i == 0:
            raise ValueError('First index cannot be removed')
        self.intervals.pop(i)
        self.slopes.pop(i)
        self._set_intercepts()

    def _set_intercepts(self):
        """Calculates the intercepts using the intervals and slopes"""
        self._intercepts = []
        for i, (interval, slope) in enumerate(zip(self.intervals,
                                                  self.slopes)):
            if i == 0:
                self._intercepts.append(0.)
            else:
                # Calculate H value at interval
                prev_intercept = self._intercepts[-1]
                prev_slope = self.slopes[i-1]
                H = prev_slope*interval + prev_intercept
                # Calculate intercept of new area of curve
                self._intercepts.append(H - slope*interval)

    def get_UoRT(self, x=0., T=c.T0('K')):
        """Calculates the excess internal energy

        Parameters
        ----------
            x : float, optional
                Coverage (in ML) of species j. Default is 0
            T : float, optional
                Temperature in K. Default is 298.15 K
        Returns
        -------
            UoRT : float
                Dimensionless internal energy
        """
        i = np.argmax(x < np.array(self.intervals))-1
        UoRT = (self.slopes[i]*x + self._intercepts[i])/(c.R('kcal/mol/K')*T)
        return UoRT

    def get_HoRT(self, x=0., T=c.T0('K')):
        """Calculates the excess enthalpy

        Parameters
        ----------
            x : float, optional
                Coverage (in ML) of species j. Default is 0
            T : float, optional
                Temperature in K. Default is 298.15 K
        Returns
        -------
            HoRT : float
                Dimensionless excess enthalpy
        """
        return self.get_UoRT(x=x, T=T)

    def get_FoRT(self, x=0., T=c.T0('K')):
        """Calculates the excess Helmholtz energy

        Parameters
        ----------
            x : float, optional
                Coverage (in ML) of species j. Default is 0
            T : float, optional
                Temperature in K. Default is 298.15 K
        Returns
        -------
            FoRT : float
                Dimensionless excess Helmholtz energy
        """
        return self.get_UoRT(x=x, T=T) - self.get_SoR()

    def get_GoRT(self, x=0., T=c.T0('K')):
        """Calculates the excess Gibbs energy

        Parameters
        ----------
            x : float, optional
                Coverage (in ML) of species j. Default is 0
            T : float, optional
                Temperature in K. Default is 298.15 K
        Returns
        -------
            GoRT : float
                Dimensionless excess Gibbs energy
        """
        return self.get_HoRT(x=x, T=T) - self.get_SoR()

    def get_q(self):
        return 1.

    def get_CvoR(self):
        return 0.

    def get_CpoR(self):
        return 0.

    def get_SoR(self):
        return 0.

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = {'class': str(self.__class__),
                    'name_i': self.name_i,
                    'name_j': self.name_j,
                    'intervals': list(self.intervals),
                    'slopes': list(self.slopes),
                    'intercepts': list(self._intercepts)}
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
            EmpiricalBase : EmpiricalBase object
        """
        json_obj = remove_class(json_obj)
        # Recalculate the intercepts to ensure range is smooth
        json_obj.pop('intercepts', None)
        return cls(**json_obj)
