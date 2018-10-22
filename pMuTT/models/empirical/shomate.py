# -*- coding: utf-8 -*-
"""
pMuTT.models.empirical.shomate

Operations related to Shomate polynomials
"""

import inspect
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import variation
from pMuTT import constants as c
from pMuTT.io_.jsonio import json_to_pMuTT, remove_class
from pMuTT.models.empirical import BaseThermo

class Shomate(BaseThermo):
    """Stores the information for an individual Shomate specie
    Inherits from pMuTT.models.empirical.BaseThermo

    The thermodynamic properties are calculated using the following form:

    :math:`\\frac{c_P}{R}=\\frac{1}{R}\\bigg(A+Bt+Ct^2+Dt^3+\\frac{E}{t^2}
    \\bigg)`
    :math:`\\frac{H}{RT}=\\frac{1}{RT}\\bigg(At+B\\frac{t^2}{2}+C\\frac{t^3}{3}
    +D\\frac{t^4}{4}-\\frac{E}{t}+F\\bigg)`
    :math:`\\frac{S}{R}=\\frac{1}{R}\\bigg(A\\ln(t)+Bt+C\\frac{t^2}{2}+D
    \\frac{t^3}{3}-\\frac{E}{2t^2}+G\\bigg)`
    
    where :math:`t=\\frac{T}{1000}` in K

    Attributes
    ----------
        T_low : float
            Lower temperature bound (in K)
        T_high : float
            High temperature bound (in K)
        a : (8,) `numpy.ndarray`_
            Shomate polynomial to use between T_low and T_high
    """
    def __init__(self, name, T_low, T_high, a, **kwargs):
        super().__init__(name=name, **kwargs)
        self.T_low = T_low
        self.T_high = T_high
        self.a = a

    def get_CpoR(self, T):
        """Calculate the dimensionless heat capacity

        Parameters
        ----------
            T : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
        Returns
        -------
            CpoR : float or (N,) `numpy.ndarray`_
                Dimensionless heat capacity

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return get_shomate_CpoR(a=self.a, T=T)

    def get_HoRT(self, T):
        """Calculate the dimensionless enthalpy

        Parameters
        ----------
            T : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
        Returns
        -------
            HoRT : float or (N,) `numpy.ndarray`_
                Dimensionless enthalpy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return get_shomate_HoRT(a=self.a, T=T)

    def get_SoR(self, T):
        """Calculate the dimensionless entropy

        Parameters
        ----------
            T : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
        Returns
        -------
            SoR : float or (N,) `numpy.ndarray`_
                Dimensionless entropy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return get_shomate_SoR(a=self.a, T=T)

    def get_GoRT(self, T):
        """Calculate the dimensionless Gibbs free energy

        Parameters
        ----------
            T : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
        Returns
        -------
            GoRT : float or (N,) `numpy.ndarray`_
                Dimensionless Gibbs free energy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return get_shomate_GoRT(a=self.a, T=T)

    @classmethod
    def from_data(cls, name, T, CpoR, T_ref, HoRT_ref, SoR_ref, **kwargs):
        """Calculates the Shomate polynomials using thermodynamic data

        Parameters
        ----------
            name : str
                Name of the species
            T : (N,) `numpy.ndarray`_
                Temperatures in K used for fitting CpoR.
            CpoR : (N,) `numpy.ndarray`_
                Dimensionless heat capacity corresponding to T.
            T_ref : float
                Reference temperature in K used fitting empirical coefficients.
            HoRT_ref : float
                Dimensionless reference enthalpy that corresponds to T_ref.
            SoR_ref : float
                Dimensionless entropy that corresponds to T_ref. 
        Returns
        -------
            shomate : Shomate object
                Shomate object with polynomial terms fitted to data.

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        T_low = min(T)
        T_high = max(T)
        a = _fit_CpoR(T=T, CpoR=CpoR)
        a = _fit_HoRT(T_ref=T_ref, HoRT_ref=HoRT_ref, a=a)
        a = _fit_SoR(T_ref=T_ref, SoR_ref=SoR_ref, a=a)
        return cls(name=name, T_low=T_low, T_high=T_high, a=a, **kwargs)

    @classmethod
    def from_statmech(cls, name, statmech_model, T_low, T_high, references=None,
                      elements=None, **kwargs):
        """Calculates the Shomate polynomial using statistical mechanic models

        Parameters
        ----------
            name : str
                Name of the species
            statmech_model : `pMuTT.models.statmech.StatMech` object or class
                Statistical Mechanics model to generate data
            T_low : float
                Lower limit temerature in K
            T_high : float
                Higher limit temperature in K
            references : `pMuTT.models.empirical.references.References` object
                Reference to adjust enthalpy
            **kwargs : keyword arguments
                Used to initalize ``statmech_model`` or ``BaseThermo`` 
                attributes to be stored.
        Returns
        -------
            shomate : Shomate object
                Shomate object with polynomial terms fitted to data.
        """
        # Initialize the StatMech object
        if inspect.isclass(statmech_model):
            statmech_model = statmech_model(**kwargs)

        # Generate data
        T = np.linspace(T_low, T_high)
        CpoR = np.array([statmech_model.get_CpoR(T=T_i) for T_i in T])
        T_ref = c.T0('K')
        HoRT_ref = statmech_model.get_HoRT(T=T_ref)
        # Add contribution of references
        if references is not None:
            HoRT_ref += references.get_HoRT_offset(elements=elements, T=T_ref)
        SoR_ref = statmech_model.get_SoR(T=T_ref)

        return cls.from_data(name=name, T=T, CpoR=CpoR, T_ref=T_ref, 
                             HoRT_ref=HoRT_ref, SoR_ref=SoR_ref,
                             statmech_model=statmech_model, elements=elements,
                             references=references, **kwargs)

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes
        
        Returns
        -------
            obj_dict : dict
        """
        obj_dict = super().to_dict()
        obj_dict['class'] = str(self.__class__)
        obj_dict['a'] = list(self.a)
        obj_dict['T_low'] = self.T_low
        obj_dict['T_high'] = self.T_high
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
            shomate : Shomate object
        """
        json_obj = remove_class(json_obj)
        # Reconstruct statmech model
        json_obj['statmech_model'] = \
                json_to_pMuTT(json_obj['statmech_model'])
        json_obj['references'] = \
                json_to_pMuTT(json_obj['references'])

        return cls(**json_obj)

def _fit_CpoR(T, CpoR):
    """Fit a[0]-a[4] coefficients given the dimensionless heat capacity data

    Parameters
    ----------
        T : (N,) `numpy.ndarray`_
            Temperatures in K
        CpoR : (N,) `numpy.ndarray`_
            Dimensionless heat capacity
    Returns
    -------
        a : (8,) `numpy.ndarray`_
            Lower coefficients of Shomate polynomial

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
    """
    # If the Cp/R does not vary with temperature (occurs when no
    # vibrational frequencies are listed), return default values
    if (np.isclose(np.mean(CpoR), 0.) and np.isnan(variation(CpoR))) \
        or np.isclose(variation(CpoR), 0.) \
        or any([np.isnan(x) for x in CpoR]):
        return np.zeros(7)
    else:
        [a, pcov] = curve_fit(_shomate_CpoR, T, np.array(CpoR))
        a = np.append(a, [0., 0., 0.])
        return a

def _fit_HoRT(T_ref, HoRT_ref, a):
    """Fit a[5] coefficient in a_low and a_high attributes given the
    dimensionless enthalpy

    Parameters
    ----------
        T_ref : float
            Reference temperature in K
        HoRT_ref : float
            Reference dimensionless enthalpy
        T_mid : float
            Temperature to fit the offset
    Returns
    -------
        a : (8,) `numpy.ndarray`_
            Lower coefficients of Shomate polynomial

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
    """
    a[5] = (HoRT_ref - get_shomate_HoRT(T=T_ref, a=a))*c.R('kJ/mol/K')*T_ref
    a[7] = - get_shomate_HoRT(T=c.T0('K'), a=a)*c.R('kJ/mol/K')*c.T0('K')
    return a

def _fit_SoR(T_ref, SoR_ref, a):
    """Fit a[6] coefficient in a_low and a_high attributes given the
    dimensionless entropy

    Parameters
    ----------
        T_ref : float
            Reference temperature in K
        SoR_ref : float
            Reference dimensionless entropy
    Returns
    -------
        a : (8,) `numpy.ndarray`_
            Lower coefficients of Shomate polynomial

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
    """
    a[6] = (SoR_ref - get_shomate_SoR(T=T_ref, a=a))*c.R('J/mol/K')
    return a

def get_shomate_CpoR(a, T):
    """Calculates the dimensionless heat capacity using Shomate polynomial form

    Parameters
    ----------
        a : (8,) `numpy.ndarray`_
            Coefficients of Shomate polynomial
        T : float or iterable
            Temperature in K
    Returns
    -------
        CpoR: float
            Dimensionless heat capacity

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
    """
    try:
        iter(T)
    except TypeError:
        T = [T]
    T = np.array(T)
    t = T/1000.
    t_arr = np.array([[1., x, x**2, x**3, 1./x**2, 0., 0., 0.] for x in t])
    return np.dot(t_arr, a)/c.R('J/mol/K')


def get_shomate_HoRT(a, T):
    """Calculates the dimensionless enthalpy using Shomate polynomial form

    Parameters
    ----------
        a : (8,) `numpy.ndarray`_
            Coefficients of Shomate polynomial
        T : float or iterable
            Temperature in K
    Returns
    -------
        HoRT : float
            Dimensionless enthalpy

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
    """
    try:
        iter(T)
    except TypeError:
        T = [T]
    T = np.array(T)
    t = T/1000.
    t_arr = np.array([[x, x**2/2., x**3/3., x**4/4., -1./x, 1., 0., 0.] 
                     for x in t])
    return np.dot(t_arr, a)/(c.R('kJ/mol/K')*T)

def get_shomate_SoR(a, T):
    """Calculates the dimensionless entropy using Shomate polynomial form

    Parameters
    ----------
        a : (8,) `numpy.ndarray`_
            Coefficients of Shomate polynomial
        T : float or iterable
            Temperature in K
    Returns
    -------
        SoR : float
            Dimensionless entropy

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
    """
    try:
        iter(T)
    except TypeError:
        T = [T]
    T = np.array(T)
    t = T/1000.
    t_arr = np.array([[np.log(x), x, x**2/2., x**3/3., -1./2./x**2, 0., 1., 0.] 
                     for x in t])
    return np.dot(t_arr, a)/c.R('J/mol/K')


def get_shomate_GoRT(a, T):
    """Calculates the dimensionless Gibbs free energy using Shomate
    polynomial form

    Parameters
    ----------
        a : (8,) `numpy.ndarray`_
            Coefficients of Shomate polynomial
        T : float or iterable
            Temperature in K
    Returns
    -------
        GoRT : float
            Dimensionless Gibbs energy

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
    """
    return get_shomate_HoRT(a=a, T=T) - get_shomate_SoR(a=a, T=T)

def _shomate_CpoR(T, A, B, C, D, E):
    """
    Helper function to fit shomate heat capacity.

    Paramters
    ---------
        T - float
            Temperature in K
        A, B, C, D, E - float
            Shomate parameters
    Returns
    -------
        CpoR - float
            Dimensionless heat capacity
    """
    a = np.array([A, B, C, D, E, 0., 0., 0.])
    return get_shomate_CpoR(a=a, T=T)