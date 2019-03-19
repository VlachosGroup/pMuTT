# -*- coding: utf-8 -*-

import numpy as np
from pMuTT import constants as c
from pMuTT.io.jsonio import remove_class


class HarmonicVib:
    """Vibrational modes using the harmonic approximation. Equations found in
    Sandler, S. I. An Introduction to Applied Statistical Thermodynamics;
    John Wiley & Sons, 2010.

    Attributes
    ----------
        vib_wavenumbers : list of float
            Vibrational wavenumbers in 1/cm
        imaginary_substitute : float, optional
            If this value is set, imaginary frequencies are substituted with
            this value for calculations. Otherwise, imaginary frequencies are
            ignored. Default is None
    """

    def __init__(self, vib_wavenumbers=[], imaginary_substitute=None):
        self.vib_wavenumbers = np.array(vib_wavenumbers)
        self.imaginary_substitute = imaginary_substitute

    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict

    def get_q(self, T, include_ZPE=True):
        """Calculates the partition function

        :math:`q^{vib}=\\prod_i \\frac{\\exp({-\\frac{\\Theta_{V,i}}{2T}})}
        {1-\\exp({-\\frac{\\Theta_{V,i}}{T}})}` if include_ZPE = True
        :math:`q^{vib}=\\prod_i \\frac{1}
        {1-\\exp({-\\frac{\\Theta_{V,i}}{T}})}` if include_ZPE = False            

        Parameters
        ----------
            T : float
                Temperature in K
            include_ZPE : bool, optional
                If True, includes the zero-point energy term
        Returns
        -------
            q_vib : float
                Vibrational partition function
        """
        vib_dimless = _get_vib_dimless(T=T,
                                       wavenumbers=self.vib_wavenumbers,
                                       substitute=self.imaginary_substitute)
        if include_ZPE:
            qs = np.array(np.exp(-vib_dimless/2.)/(1. - np.exp(-vib_dimless)))
        else:
            qs = np.array(1./(1. - np.exp(-vib_dimless)))
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
        vib_dimless = _get_vib_dimless(T=T,
                                       wavenumbers=self.vib_wavenumbers,
                                       substitute=self.imaginary_substitute)
        CvoRs = np.array([(0.5*vib_dimless)**2 *
                          (1./np.sinh(vib_dimless/2.))**2])
        return np.sum(CvoRs)

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

        :math:`ZPE=\\frac{1}{2}k_b\\sum_i \\Theta_{V,i}`

        Returns
        -------
            zpe : float
                Zero point energy in eV
        """
        valid_wavenumbers = _get_valid_vib_wavenumbers(
                wavenumbers=self.vib_wavenumbers,
                substitute=self.imaginary_substitute)
        vib_temperatures = c.wavenumber_to_temp(valid_wavenumbers)
        return 0.5*c.kb('eV/K')*np.sum(vib_temperatures)

    def get_UoRT(self, T):
        """Calculates the dimensionless internal energy

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
        vib_dimless = _get_vib_dimless(T=T,
                                       wavenumbers=self.vib_wavenumbers,
                                       substitute=self.imaginary_substitute)
        UoRT = np.array([vib_dimless/2. + vib_dimless*np.exp(-vib_dimless)
                        / (1.-np.exp(-vib_dimless))])
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
        vib_dimless = _get_vib_dimless(T=T,
                                       wavenumbers=self.vib_wavenumbers,
                                       substitute=self.imaginary_substitute)
        SoR = np.array([
                vib_dimless*np.exp(-vib_dimless)/(1.-np.exp(-vib_dimless))
                - np.log(1. - np.exp(-vib_dimless))])
        return np.sum(SoR)

    def get_FoRT(self, T):
        """Calculates the dimensionless Helmholtz energy

        :math:`\\frac{A^{vib}}{RT}=\\frac{U^{vib}}{RT}-\\frac{S^{vib}}{R}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            FoRT_vib : float
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

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        return {'class': str(self.__class__),
                'vib_wavenumbers': list(self.vib_wavenumbers),
                'imaginary_substitute': self.imaginary_substitute}

    @classmethod
    def from_dict(cls, json_obj):
        """Recreate an object from the JSON representation.

        Parameters
        ----------
            json_obj : dict
                JSON representation
        Returns
        -------
            HarmonicVib : HarmonicVib object
        """
        json_obj = remove_class(json_obj)
        return cls(**json_obj)

    def print_calc_wavenumbers(self):
        """Prints the wavenumbers that will be used in a thermodynamic
        calculation. If ``self.imaginary_substitute`` is a float, then
        imaginary frequencies are replaced with that value. Otherwise,
        imaginary frequencies are ignored."""
        print(_get_valid_vib_wavenumbers(wavenumbers=self.vib_wavenumbers,
                                         substitute=self.imaginary_substitute))


class QRRHOVib:
    """Vibrational modes using the Quasi Rigid Rotor Harmonic Oscillator
    approximation. Equations found in

    1. Li, Y. P.; Gomes, J.; Sharada, S. M.; Bell, A. T.; Head-Gordon, M. J.
    Phys. Chem. C 2015, 119 (4), 1840–1850.

    2. Grimme, S. Chem. - A Eur. J. 2012, 18 (32), 9955–9964.

    Attributes
    ----------
        vib_wavenumber : list of float
            Vibrational wavenumbers in 1/cm
        Bav : float, optional
            Average molecular moment of inertia as a limiting value of small
            wavenumbers. Default is 1.e-44 kg m2
        v0 : float, optional
            Wavenumber to scale vibrations. Default is 100 cm :sup:`-1`
        alpha : int, optional
            Power to raise ratio of wavenumbers. Default is 4
        imaginary_substitute : float, optional
            If this value is set, imaginary frequencies are substituted with
            this value for calculations. Otherwise, imaginary frequencies are
            ignored. Default is None
    """

    def __init__(self, vib_wavenumbers, Bav=1.e-44, v0=100., alpha=4,
                 imaginary_substitute=None):
        self.Bav = Bav
        self.v0 = v0
        self.vib_wavenumbers = vib_wavenumbers
        self.alpha = alpha
        self.imaginary_substitute = imaginary_substitute

    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict

    def _get_scaled_wavenumber(self, wavenumber):
        """Calculates the scaled wavenumber determining mixture of RRHO to
        add.

        :math:`\\omega = \\frac {1}{1 + (\\frac{\\nu_0}{\\nu})^\\alpha}`

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

        :math:`\\mu'=\\frac {\\mu B_{av}} {\\mu + B_{av}}`

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

        :math:`\\frac {C_{v}^{qRRHO}}{R} = \\sum_{i}\\omega_i\\frac{C_{v,i}
        ^{RRHO}}{R} + \\frac{1}{2}(1-\\omega_i)`

        :math:`\\frac{C_{v}^{RRHO}}{R} = \\sum_{i}\\exp \\bigg(-\\frac{
        \\Theta_i}{T}\\bigg) \\bigg(\\frac{\\Theta_i}{T}\\frac{1}{1-\\exp(-
        \\frac{\\Theta_i}{T})}\\bigg)^2`

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

        valid_wavenumbers = _get_valid_vib_wavenumbers(
                wavenumbers=self.vib_wavenumbers,
                substitute=self.imaginary_substitute)
        vib_temperatures = c.wavenumber_to_temp(valid_wavenumbers)
        vib_dimless = vib_temperatures/T
        scaled_wavenumbers = self._get_scaled_wavenumber(valid_wavenumbers)
        for vib_dimless_i, w_i in zip(vib_dimless, scaled_wavenumbers):
            CvoR_RRHO = np.exp(-vib_dimless_i) \
                        * (vib_dimless_i/(1. - np.exp(-vib_dimless_i)))**2
            CvoR.append(w_i*CvoR_RRHO + 0.5*(1.-w_i))
        return np.sum(CvoR)

    def get_CpoR(self, T):
        """Calculates the dimensionless heat capacity at constant pressure

        :math:`\\frac{C_{P}^{qRRHO}} {R} = \\frac{C_{V}^{qRRHO}} {R}`

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

        :math:`ZPE=\\frac{1}{2}k_b\\sum_i \\omega_i\\Theta_{V,i}`

        Returns
        -------
            zpe : float
                Zero point energy in eV
        """
        zpe = []
        valid_wavenumbers = _get_valid_vib_wavenumbers(
                wavenumbers=self.vib_wavenumbers,
                substitute=self.imaginary_substitute)
        vib_temperatures = c.wavenumber_to_temp(valid_wavenumbers)
        scaled_wavenumbers = self._get_scaled_wavenumber(valid_wavenumbers)

        for theta_i, w_i in zip(vib_temperatures, scaled_wavenumbers):
            zpe = 0.5*c.kb('eV/K')*theta_i*w_i
        return np.sum(zpe)

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
        """Calculates the dimensionless internal energy

        :math:`\\frac {U^{qRRHO}}{RT} = \\sum_{i}\\omega_i\\frac{U^{RRHO}}{RT}
        + \\frac{1}{2}(1-\\omega_i)`

        :math:`\\frac {U^{RRHO}_{i}}{RT} = \\frac{\\Theta_i}{T} \\bigg(
        \\frac{1}{2} + \\frac{\\exp(-\\frac{\\Theta_i}{T})}{1-\\exp(-\\frac{
        \\Theta_i}{T})}\\bigg)`

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
        valid_wavenumbers = _get_valid_vib_wavenumbers(
                wavenumbers=self.vib_wavenumbers,
                substitute=self.imaginary_substitute)
        vib_temperatures = c.wavenumber_to_temp(valid_wavenumbers)
        scaled_wavenumbers = self._get_scaled_wavenumber(valid_wavenumbers)

        for theta_i, w_i in zip(vib_temperatures, scaled_wavenumbers):
            UoRT_RRHO = self._get_UoRT_RRHO(T=T, vib_temperature=theta_i)
            UoRT_QRRHO.append(w_i*UoRT_RRHO + (1.-w_i)*0.5)
        return np.sum(UoRT_QRRHO)

    def get_HoRT(self, T):
        """Calculates the dimensionless enthalpy

        :math:`\\frac{H^{qRRHO}} {RT} = \\frac{U^{qRRHO}} {RT}`

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
            - np.log(1-np.exp(-vib_temperature/T))

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
                            / c.h('J s')**2)**0.5)

    def get_SoR(self, T):
        """Calculates the dimensionless entropy

        :math:`\\frac{S^{qRRHO}}{R}=\\sum_i\\omega_i\\frac{S_i^{H}}{R}+(1-
        \\omega_i)\\frac{S_i^{RRHO}}{R}`

        :math:`\\frac {S^{RRHO}_i}{R} = \\frac{1}{2} + \\log \\bigg(\\bigg[
        \\frac{8\\pi^3\\mu'_ik_BT}{h^2}\\bigg]^{\\frac{1}{2}}\\bigg)`

        :math:`\\frac {S^{H}_i}{R}=\\bigg(\\frac{\\Theta_i}{T}\\bigg)\\frac{1}
        {\\exp(\\frac{\\Theta_i}{T})-1}-\\log\\bigg(1-\\exp(\\frac{-\\Theta_i}
        {T})\\bigg)`

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

        valid_wavenumbers = _get_valid_vib_wavenumbers(
                wavenumbers=self.vib_wavenumbers,
                substitute=self.imaginary_substitute)
        vib_temperatures = c.wavenumber_to_temp(valid_wavenumbers)
        scaled_inertia = self._get_scaled_inertia(valid_wavenumbers)
        scaled_wavenumbers = self._get_scaled_wavenumber(valid_wavenumbers)

        for theta_i, mu_i, w_i in zip(vib_temperatures, scaled_inertia,
                                      scaled_wavenumbers):
            SoR_H = self._get_SoR_H(T=T, vib_temperature=theta_i)
            SoR_RRHO = self._get_SoR_RRHO(T=T, vib_inertia=mu_i)
            SoR_QRRHO.append(w_i*SoR_H + (1.-w_i)*SoR_RRHO)
        return np.sum(SoR_QRRHO)

    def get_FoRT(self, T):
        """Calculates the dimensionless Helmholtz energy

        :math:`\\frac{A^{qRRHO}}{RT} = \\frac{U^{qRRHO}}{RT}-
        \\frac{S^{qRRHO}}{R}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            FoRT_vib : float
                Vibrational dimensionless Helmholtz energy
        """
        return self.get_UoRT(T=T) - self.get_SoR(T=T)

    def get_GoRT(self, T):
        """Calculates the dimensionless Gibbs energy

        :math:`\\frac{G^{qRRHO}}{RT} = \\frac{H^{qRRHO}}{RT}-
        \\frac{S^{qRRHO}}{R}`

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

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        return {'class': str(self.__class__),
                'vib_wavenumbers': list(self.vib_wavenumbers),
                'Bav': self.Bav,
                'v0': self.v0,
                'alpha': self.alpha,
                'imaginary_substitute': self.imaginary_substitute}

    @classmethod
    def from_dict(cls, json_obj):
        """Recreate an object from the JSON representation.

        Parameters
        ----------
            json_obj : dict
                JSON representation
        Returns
        -------
            QRRHOVib : QRRHOVib object
        """
        json_obj = remove_class(json_obj)
        return cls(**json_obj)

    def print_calc_wavenumbers(self):
        """Prints the wavenumbers that will be used in a thermodynamic
        calculation. If ``self.imaginary_substitute`` is a float, then
        imaginary frequencies are replaced with that value. Otherwise,
        imaginary frequencies are ignored."""
        print(_get_valid_vib_wavenumbers(wavenumbers=self.vib_wavenumbers,
                                         substitute=self.imaginary_substitute))


class EinsteinVib:
    """Einstein model of a crystal. Equations found in
    Sandler, S. I. An Introduction to Applied Statistical Thermodynamics;
    John Wiley & Sons, 2010.

    Attributes
    ----------
        einstein_temperature : float
            Einstein temperature (:math:`\\Theta_E`) in K
        interaction_energy : float, optional
            Interaction energy (:math:`u`) per atom in eV. Default is 0 eV
    """

    def __init__(self, einstein_temperature, interaction_energy=0.):
        self.einstein_temperature = einstein_temperature
        self.interaction_energy = interaction_energy

    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict

    def get_q(self, T):
        """Calculates the partition function

        :math:`q^{vib}=\\exp\\bigg({\\frac{-u}{k_BT}}\\bigg)\\bigg(\\frac{
        \\exp(-\\frac{\\Theta_E}{2T})}{1-\\exp(\\frac{-\\Theta_E}{T})}\\bigg)`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            q_vib : float
                Vibrational partition function
        """
        u = self.interaction_energy
        theta_E = self.einstein_temperature
        return np.exp(-u/c.kb('eV/K')/T) \
            * (np.exp(-theta_E/2./T)/(1. - np.exp(-theta_E/T)))

    def get_CvoR(self, T):
        """Calculates the dimensionless heat capacity at constant volume

        :math:`\\frac{C_V^{vib}}{R}=3\\bigg(\\frac{\\Theta_E}{T}\\bigg)^2
        \\frac{\\exp(-\\frac{\\Theta_E}{T})}{\\big(1-\\exp(\\frac{-
        \\Theta_E}{T})\\big)^2}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            CvoR_vib : float
                Vibrational dimensionless heat capacity at constant volume
        """
        theta_E = self.einstein_temperature
        return 3.*(theta_E/T)**2*np.exp(-theta_E/T)/(1-np.exp(-theta_E/T))**2

    def get_CpoR(self, T):
        """Calculates the dimensionless heat capacity at constant pressure

        :math:`\\frac{C_P^{vib}}{R}=\\frac{C_V^{vib}}{R}=3\\bigg(\\frac{
        \\Theta_E}{T}\\bigg)^2\\frac{\\exp(-\\frac{\\Theta_E}{T})}{\\big(1-
        \\exp(\\frac{-\\Theta_E}{T})\\big)^2}`

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

        :math:`u^0_E=u+\\frac{3}{2}\\Theta_E k_B`

        Returns
        -------
            ZPE : float
                Zero point energy in eV
        """
        return self.interaction_energy \
            + 1.5*self.einstein_temperature*c.kb('eV/K')

    def get_UoRT(self, T):
        """Calculates the dimensionless internal energy

        :math:`\\frac{U^{vib}}{RT}=\\frac{u^0_E}{k_BT}+3\\frac{\\Theta_E}{T}
        \\bigg(\\frac{\\exp(-\\frac{\\Theta_E}{T})}{1-\\exp(-\\frac{\\Theta_E}
        {T})}\\bigg)`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            UoRT_vib : float
                Vibrational dimensionless internal energy
        """
        theta_E = self.einstein_temperature
        return self.get_ZPE()/c.kb('eV/K')/T \
            + 3.*theta_E/T*np.exp(-theta_E/T)/(1. - np.exp(-theta_E/T))

    def get_HoRT(self, T):
        """Calculates the dimensionless enthalpy

        :math:`\\frac{H^{vib}}{RT}=\\frac{U^{vib}}{RT}=\\frac{N_A u^0_E}{k_BT}
        +3\\frac{\\Theta_E}{T}\\bigg(\\frac{\\exp(-\\frac{\\Theta_E}{T})}{1-
        \\exp(-\\frac{\\Theta_E}{T})}\\bigg)`

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

        :math:`\\frac{S^{vib}}{R}=3\\bigg(\\frac{\\Theta_E}{T}\\frac{\\exp\\big(
        \\frac{-\\Theta_E}{T}\\big)}{1-\\exp\\big(-\\frac{\\Theta_E}{T}\\big)}
        \\bigg)-\\ln\\bigg(1-\\exp\\big(\\frac{-\\Theta_E}{T}\\big)\\bigg)`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            SoR_vib : float
                Vibrational dimensionless entropy
        """
        theta_E = self.einstein_temperature
        exp_term = np.exp(-theta_E/T)
        return 3.*(theta_E/T*exp_term/(1. - exp_term) - np.log(1. - exp_term))

    def get_FoRT(self, T):
        """Calculates the dimensionless Helmholtz energy

        :math:`\\frac{A^{vib}}{RT}=\\frac{U^{vib}}{RT}-\\frac{S^{vib}}{R}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            FoRT_vib : float
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

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        return {'class': str(self.__class__),
                'einstein_temperature': self.einstein_temperature,
                'interaction_energy': self.interaction_energy}

    @classmethod
    def from_dict(cls, json_obj):
        """Recreate an object from the JSON representation.

        Parameters
        ----------
            json_obj : dict
                JSON representation
        Returns
        -------
            HarmonicVib : HarmonicVib object
        """
        json_obj = remove_class(json_obj)
        return cls(**json_obj)


def debye_to_einstein(debye_temperature):
    """Converts Debye temperature to Einstein temperature

    Parameters
    ----------
        debye_temperature : float
            Debye temperature in K
    Returns
    -------
        einstein_temperature : float
            Einstein temperature in K
    """
    return (np.pi/6.)**(1./3.)*debye_temperature


def einstein_to_debye(einstein_temperature):
    """Converts Einstein temperature to Debye temperature

    Parameters
    ----------
        einstein_temperature : float
            Einstein temperature in K
    Returns
    -------
        debye_temperature : float
            Debye temperature in K
    """
    return einstein_temperature/(np.pi/6.)**(1./3.)


def _get_valid_vib_wavenumbers(wavenumbers, substitute=None):
    """Returns wavenumbers to use for vibration calculations. Imaginary
    frequencies are expected to be negative.

    Parameters
    ----------
        wavenumbers : list of float
            Wavenumbers in 1/cm
        substitute : float, optional
            Value to use to replace imaginary frequencies. If not specified,
            imaginary frequencies are ignored. Default is None
    Returns
    -------
        wavenumbers_out : (N,) np.ndarray
            Valid wavenumbers
    """
    wavenumbers_out = []
    for wavenumber in wavenumbers:
        if wavenumber > 0.:
            # Real wavenumbers always added
            wavenumbers_out.append(wavenumber)
        elif substitute is not None:
            # Substitute added if imaginary frequency encountered
            wavenumbers_out.append(substitute)
    return np.array(wavenumbers_out)


def _get_vib_dimless(wavenumbers, T, substitute=None):
    """Calculates dimensionless temperatures for the wavenumbers and
    temperature specified

    Parameters
    ----------
        wavenumbers : (N,) np.ndarray
            Wavenumbers in 1/cm
        T : float
            Temperature in K
        substitute : float, optional
            Value to use to replace imaginary frequencies. If not specified,
            imaginary frequencies are ignored. Default is None
    Returns
    -------
        vib_dimless : (N,) np.ndarray
            Vibrational temperatures normalized by T
    """
    valid_wavenumbers = _get_valid_vib_wavenumbers(wavenumbers=wavenumbers,
                                                   substitute=substitute)
    vib_dimless = c.wavenumber_to_temp(valid_wavenumbers)/T
    return vib_dimless
