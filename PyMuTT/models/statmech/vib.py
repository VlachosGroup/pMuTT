# -*- coding: utf-8 -*-

import numpy as np
from PyMuTT import constants as c

class HarmonicVib:
    """Vibrational modes using the harmonic approximation.
    
    Attributes
    ----------
        vib_wavenumbers : list of float
            Vibrational wavenumbers in 1/cm
    """

    def __init__(self, vib_wavenumbers):
        self.vib_wavenumbers = vib_wavenumbers
        self._vib_temperatures = [c.wavenumber_to_temp(wavenumber) 
            for wavenumber in vib_wavenumbers]

    def get_q(self, T):
        """Calculates the partition function

        :math:`q^{vib}=\\prod_i \\frac{\\exp({-\\frac{\\Theta_{V,i}}{2T}})}
        {1-\\exp({-\\frac{\\Theta_{V,i}}{T}})}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            q_vib : float
                Vibrational partition function
        """
        qs = []
        for vib_temperature in self._vib_temperatures:
            vib_dimless = vib_temperature/T            
            qs.append(np.exp(-vib_dimless/2.)/(1. - np.exp(-vib_dimless)))
        return np.prod(qs)

    def get_CvoR(self, T):
        """Calculates the dimensionless heat capacity at constant volume

        :math:`\\frac{C_V^{vib}}{R}=\\sum_i \\bigg(\\frac{\\Theta_{V,i}}{2T}
        \\bigg)^2 \\frac{1}{\\big(\\sinh{\\frac{\\Theta_{V,i}}{2T}}\\big)^2}`
        
        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            CvoR_vib : float
                Vibrational dimensionless heat capacity at constant volume
        """
        CvoR = []
        for vib_temperature in self._vib_temperatures:
            vib_dimless = vib_temperature/T
            CvoR.append((0.5*vib_dimless)**2*(1./np.sinh(vib_dimless/2.))**2)
        return np.sum(CvoR)

    def get_CpoR(self, T):
        """Calculates the dimensionless heat capacity at constant pressure

        :math:`\\frac{C_P^{vib}}{R}=\\frac{C_V^{vib}}{R}=\\sum_i \\bigg(\\frac{
        \\Theta_{V,i}}{2T}\\bigg)^2 \\frac{1}{\\big(\\sinh{\\frac{\\Theta_{V,i}}
        {2T}}\\big)^2}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            CpoR_vib : float
                Vibrational dimensionless heat capacity at constant pressure
        """
        return self.get_CvoR(T=T)
    
    def get_ZPE(self):
        """Calculates the zero point energy

        Returns
        -------
            zpe : float
                Zero point energy in eV
        """
        return 0.5*c.kb('eV/K')*np.sum(self._vib_temperatures)

    def get_UoRT(self, T):
        """Calculates the imensionless internal energy

        :math:`\\frac{U^{vib}}{RT}=\\sum_i \\bigg(\\frac{\\Theta_{V,i}}{2T}+
        \\frac{\\Theta_{V,i}}{T}\\frac{\\exp\\big(-\\frac{\\Theta_{V,i}}{T}
        \\big)}{1-\\exp\\big(-\\frac{\\Theta_{V_i}}{T}\\big)}\\bigg)`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            UoRT_vib : float
                Vibrational dimensionless internal energy
        """
        UoRT = []
        for vib_temperature in self._vib_temperatures:
            vib_dimless = vib_temperature/T
            UoRT.append(
                vib_dimless/2. 
                + vib_dimless*np.exp(-vib_dimless)/(1.-np.exp(-vib_dimless)))
        return np.sum(UoRT)

    def get_HoRT(self, T):
        """Calculates the dimensionless enthalpy

        :math:`\\frac{H^{vib}}{RT}=\\frac{U^{vib}}{RT}=\\sum_i \\bigg(\\frac{
        \\Theta_{V,i}}{2T}+\\frac{\\Theta_{V,i}}{T}\\frac{\\exp\\big(-\\frac{
        \\Theta_{V,i}}{T}\\big)}{1-\\exp\\big(-\\frac{\\Theta_{V_i}}{T}\\big)}
        \\bigg)`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            HoRT_vib : float
                Vibrational dimensionless enthalpy
        """
        return self.get_UoRT(T=T)

    def get_SoR(self, T):
        """Calculates the dimensionless entropy

        :math:`\\frac{S^{vib}}{R}=\\sum_i \\frac{\\Theta_{V,i}}{T}\\frac{\\exp
        \\big(-\\frac{\\Theta_{V,i}}{T}\\big)}{1-\\exp\\big(-\\frac{
        \\Theta_{V,i}}{T}\\big)}-\\ln \\bigg(1-\\exp\\big(-\\frac{
        \\Theta_{V,i}}{T}\\big)\\bigg)`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            SoR_vib : float
                Vibrational dimensionless entropy
        """
        SoR = []
        for vib_temperature in self._vib_temperatures:
            vib_dimless = vib_temperature/T
            SoR.append(
                vib_dimless*np.exp(-vib_dimless)/(1.-np.exp(-vib_dimless))
                - np.log(1. - np.exp(-vib_dimless)))
        return np.sum(SoR)

    def get_AoRT(self, T):
        """Calculates the dimensionless Helmholtz energy

        :math:`\\frac{A^{vib}}{RT}=\\frac{U^{vib}}{RT}-\\frac{S^{vib}}{R}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            AoRT_vib : float
                Vibrational dimensionless Helmholtz energy
        """
        return self.get_UoRT(T=T) - self.get_SoR(T=T)

    def get_GoRT(self, T):
        """Calculates the dimensionless Gibbs energy

        :math:`\\frac{G^{vib}}{RT}=\\frac{H^{vib}}{RT}-\\frac{S^{vib}}{R}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            GoRT_vib : float
                Vibrational dimensionless Gibbs energy
        """
        return self.get_HoRT(T=T) - self.get_SoR(T=T)