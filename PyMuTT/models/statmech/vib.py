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
        self._vib_temperatures = [
            c.wavenumber_to_temp(wavenumber) for wavenumber in vib_wavenumbers]

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
                vib_dimless/2. \
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
                vib_dimless*np.exp(-vib_dimless)/(1.-np.exp(-vib_dimless)) \
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

class QRRHOVib:
    """Vibrational modes using the Quasi Rigid Rotor Harmonic Oscillator
    approximation.
    
    Attributes
    ----------
        vib_wavenumber : list of float
            Vibrational wavenumbers in 1/cm
        Bav : float, optional
            Average molecular moment of inertia as a limiting value of small
            wavenumbers. Default is 1.e-44 kg m2
        v0 : float
        alpha : int, optional
            Power to raise ratio of wavenumbers. Default is 4
        harmonic_model : `PyMuTT.models.statmech.vib.HarmonicVib` object
            Harmonic approximation contribution
        rigidrotor_model: `PyMuTT.models.statmech.vib.RigidRotor` object
            Rigid rotor approximation contribution
    """

    def __init__(self, vib_wavenumbers, Bav=1.e-44, v0=100., alpha=4):
        self.Bav = Bav
        self.v0 = v0
        self.vib_wavenumbers = vib_wavenumbers
        self.alpha = alpha
        self._vib_temperatures = [c.wavenumber_to_temp(wavenumber) \
                                  for wavenumber in vib_wavenumbers]
        self._scaled_wavenumbers = [self._get_scaled_wavenumber(vi) \
                                    for vi in vib_wavenumbers]
        self._inertia = [c.wavenumber_to_inertia(wavenumber) \
                         for wavenumber in vib_wavenumbers]
        self._scaled_inertia = [self._get_scaled_inertia(vi) \
                                for vi in vib_wavenumbers]

    def _get_scaled_wavenumber(self, wavenumber):
        """Calculates the scaled wavenumber determining mixture of RRHO to
        add.

        Parameters
        ----------
            wavenumber : float
                Vibrational wavenumber in 1/cm
        Returns
        -------
            scaled_wavenumber : float
                Scaled wavenumber
        """
        return 1./(1. + (self.v0/wavenumber)**self.alpha)
    
    def _get_scaled_inertia(self, mu):
        """Calculates the scaled moment of inertia.

        Parameters
        ----------
            mu : float
                Moment of inertia in kg*m2
        Returns
        -------
            mu1 : float
                Scaled moment of inertia in kg*m2
        """
        return mu*self.Bav/(mu + self.Bav)

    def get_q(self):
        """Calculates the partition function

        Returns
        -------
            q_vib : float
                Vibrational partition function
        """
        raise NotImplementedError()

    def get_CvoR(self, T):
        """Calculates the dimensionless heat capacity at constant volume

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
        for theta_i, w_i in zip(self._vib_temperatures, self._scaled_wavenumbers):
            CvoR_RRHO = (theta_i/T/(1. - np.exp(-theta_i/T)))**2
            CvoR.append(w_i*CvoR_RRHO + 0.5*(1.-w_i))
        return np.sum(CvoR)

    def get_CpoR(self, T):
        """Calculates the dimensionless heat capacity at constant pressure

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
    
    def _get_UoRT_RRHO(self, T, vib_temperature):
        """Calculates the dimensionless RRHO contribution to internal energy
        
        Parameters
        ----------
            T : float
                Temperature in K
            vib_temperature : float
                Vibrational temperature in K
        Returns
        -------
            UoRT_RRHO : float
                Dimensionless internal energy of Rigid Rotor Harmonic Oscillator
        """
        vib_dimless = vib_temperature/T
        return vib_dimless*(0.5+np.exp(-vib_dimless)/(1.-np.exp(-vib_dimless)))

    def get_UoRT(self, T):
        """Calculates the imensionless internal energy

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            UoRT_vib : float
                Vibrational dimensionless internal energy
        """
        UoRT_QRRHO = []
        for theta_i, w_i in \
                zip(self._vib_temperatures, self._scaled_wavenumbers):
            UoRT_RRHO = self._get_UoRT_RRHO(T=T, vib_temperature=theta_i)
            UoRT_QRRHO.append(w_i*UoRT_RRHO + (1.-w_i)*0.5)
        return np.sum(UoRT_QRRHO)

    def get_HoRT(self, T):
        """Calculates the dimensionless enthalpy

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

    def _get_SoR_H(self, T, vib_temperature):
        """Calculates the dimensionless harmonic osccilator contribution to 
        entropy

        Parameters
        ----------
            T : float
                Temperature in K
            vib_temperature : float
                Vibrational temperature in K
        Returns
        -------
            SoR_RHHO : float
                Dimensionless entropy of Rigid Rotor Harmonic Oscillator
        """
        return vib_temperature/T/(np.exp(vib_temperature/T)-1) \
              -np.log(1-np.exp(-vib_temperature/T))

    def _get_SoR_RRHO(self, T, vib_inertia):
        """Calculates the dimensionless RRHO contribution to entropy

        Parameters
        ----------
            T : float
                Temperature in K
            vib_inertia : float
                Vibrational inertia in kg m2
        Returns
        -------
            SoR_RHHO : float
                Dimensionless entropy of Rigid Rotor Harmonic Oscillator
        """
        return 0.5 + np.log((8.*np.pi**3*vib_inertia*c.kb('J/K')*T
                     /c.h('J s')**2)**0.5)

    def get_SoR(self, T):
        """Calculates the dimensionless entropy

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            SoR_vib : float
                Vibrational dimensionless entropy
        """
        SoR_QRRHO = []
        for theta_i, mu_i, w_i in zip(self._vib_temperatures, self._scaled_inertia, self._scaled_wavenumbers):
            SoR_H = self._get_SoR_H(T=T, vib_temperature=theta_i)
            SoR_RRHO = self._get_SoR_RRHO(T=T, vib_inertia=mu_i)
            SoR_QRRHO.append(w_i*SoR_H + (1.-w_i)*SoR_RRHO)
        return np.sum(SoR_QRRHO)

    def get_AoRT(self, T):
        """Calculates the dimensionless Helmholtz energy

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