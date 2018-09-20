# -*- coding: utf-8 -*-
"""
PyMuTT.models.statmech.harmonicthermo

Vlachos group code for Harmonic approximation.
"""

from warnings import warn
import numpy as np
from ase import thermochemistry
from PyMuTT import constants as c
from PyMuTT.models.statmech.heat_capacity import get_CvoR_vib


class HarmonicThermo:
    """Treats all degrees of freedom harmonically.
    Uses ase.thermochemistry.HarmonicThermo to calculate enthalpy and entropy

    Attributes
    ----------
        model : `ase.thermochemistry.HarmonicThermo`_

    Below are parameters accepted by `ase.thermochemistry.HarmonicThermo`_

    ================= =======================================================
    Parameter         Description
    ================= =======================================================
    vib_energies      ((3N,) `numpy.ndarray`_ where N is number of atoms)
                      Vibrational energies in eV
    potentialenergy   (float) Potential energy in eV
    ================= =======================================================

    .. _`ase.thermochemistry.HarmonicThermo`: https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html#ase.thermochemistry.HarmonicThermo
    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
    """
    def __init__(self, vib_energies, potentialenergy=0.0):  
        warn('Class replaced with StatMech class. To replicate HarmonicThermo '
             'behavior, use PyMuTT.models.statmech.presents.',
             DeprecationWarning)
        self.model = thermochemistry.HarmonicThermo(
                vib_energies=vib_energies,
                potentialenergy=potentialenergy)

    def get_CpoR(self, Ts):
        """Calculates the dimensionless heat capacity (Cp/R) at a given
        temperature.

        Parameters
        ----------
            Ts : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
        Returns
        -------
            CpoR : float or (N,) `numpy.ndarray`_
                Dimensionless heat capacity (Cp/R)

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return get_CvoR_vib(vib_energies=self.model.vib_energies, Ts=Ts)

    def get_HoRT(self, Ts, verbose=False):
        """Calculates the dimensionless enthalpy at a given temperature

        Parameters
        ----------
            Ts : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
            verbose : bool
                Whether a table breaking down each contribution should
                be printed
        Returns
        -------
            HoRT : float or (N,) `numpy.ndarray`_
                Dimensionless heat capacity (H/RT) at the specified temperature

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        try:
            iter(Ts)
        except TypeError:
                HoRT = self.model.get_internal_energy(
                        temperature=Ts, verbose=verbose)/(c.kb('eV/K') * Ts)
        else:
            HoRT = np.zeros_like(Ts)
            for i, T in enumerate(Ts):
                HoRT[i] = self.model.get_internal_energy(
                        temperature=T, verbose=verbose)/(c.kb('eV/K') * T)
        return HoRT

    def get_SoR(self, Ts, verbose=False):
        """Returns the dimensionless entropy at a given temperature and
        pressure

        Parameters
        ----------
            Ts : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
            verbose : bool
                Whether a table breaking down each contribution should be
                printed
        Returns
        -------
            SoR : float or (N,) `numpy.ndarray`_
                Dimensionless entropy (S/R) at the specified temperature and
                pressure

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        try:
            iter(Ts)
        except TypeError:
            SoR = self.model.get_entropy(temperature=Ts,
                                         verbose=verbose)/c.R('eV/K')
        else:
            SoR = np.zeros_like(Ts)
            for i, T in enumerate(Ts):
                SoR[i] = self.model.get_entropy(temperature=T,
                                                verbose=verbose)/c.R('eV/K')
        return SoR

    def get_GoRT(self, Ts, verbose=False):
        """Returns the dimensionless Gibbs energy at a given temperature

        Parameters
        ----------
            Ts : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
            verbose : bool
                Whether a table breaking down each contribution should be
                printed
        Returns
        -------
            GoRT : float or (N,) `numpy.ndarray`_
                Dimensionless heat capacity (G/RT) at the specified temperature

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        try:
            iter(Ts)
        except TypeError:
            GoRT = self.model.get_helmholtz_energy(temperature=Ts,
                                                   verbose=verbose) /\
                                                   (c.kb('eV/K') * Ts)
        else:
            GoRT = np.zeros_like(Ts)
            for i, T in enumerate(Ts):
                GoRT[i] = self.model.get_helmholtz_energy(temperature=T,
                                                          verbose=verbose) /\
                                                          (c.kb('eV/K') * T)
        return GoRT
