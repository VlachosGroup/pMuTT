# -*- coding: utf-8 -*-
"""
pmutt.empirical.shomate

Operations related to Shomate polynomials
"""

import inspect
from warnings import warn

import numpy as np
from scipy.optimize import curve_fit

from pmutt import _get_R_adj, _is_iterable
from pmutt import constants as c
from pmutt.empirical import EmpiricalBase
from pmutt.io.cantera import obj_to_cti
from pmutt.io.json import json_to_pmutt, remove_class
from pmutt.mixture import _get_mix_quantity


class Shomate(EmpiricalBase):
    """Stores the information for an individual Shomate specie
    Inherits from :class:`~pmutt.empirical.EmpiricalBase`

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
        units : str, optional
            Units used to fit the Shomate polynomial. Units should be supported
            by :class:`~pmutt.constants.R` (e.g. J/mol/K, cal/mol/K, eV/K).
            Default is J/mol/K.
    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
    """

    def __init__(self,
                 name,
                 T_low,
                 T_high,
                 a,
                 units='J/mol/K',
                 n_sites=None,
                 **kwargs):
        super().__init__(name=name, **kwargs)
        self.T_low = T_low
        self.T_high = T_high
        self.a = a
        self.units = units
        self.n_sites = n_sites

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, val):
        try:
            c.R(val)
        except KeyError:
            err_msg = ('Units, "{}", inputted into pmutt.empirical.Shomate '
                       'object are invalid. See pmutt.constants.R for '
                       'supported units.'.format(val))
            raise ValueError(err_msg)
        else:
            self._units = val

    def get_CpoR(self, T, raise_error=True, raise_warning=True, **kwargs):
        """Calculate the dimensionless heat capacity

        Parameters
        ----------
            T : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Arguments to calculate mixture model properties, if any
        Returns
        -------
            CpoR : float or (N,) `numpy.ndarray`_
                Dimensionless heat capacity

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        if isinstance(T, np.ndarray):
            T = T.astype(float)
        else:
            T = float(T)
        # Convert T to 1D numpy format
        if not _is_iterable(T):
            T = [T]
        T = np.array(T)
        self._check_T(T)

        # Calculate pure properties
        CpoR = get_shomate_CpoR(a=self.a, T=T, units=self.units)
        # Calculate mixing properties
        for T_i in T:
            CpoR_mix = _get_mix_quantity(misc_models=self.misc_models,
                                         method_name='get_CpoR',
                                         raise_error=raise_error,
                                         raise_warning=raise_warning,
                                         default_value=0.,
                                         T=T_i,
                                         **kwargs)
        # Add mixing quantity in appropriate format
        CpoR = CpoR + CpoR_mix
        if len(T) == 1:
            CpoR = CpoR.item(0)
        return CpoR

    def get_Cp(self, T, units, raise_error=True, raise_warning=True, **kwargs):
        """Calculate the heat capacity

        Parameters
        ----------
            T : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
                units.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Arguments to calculate mixture model properties, if any
        Returns
        -------
            Cp : float or (N,) `numpy.ndarray`_
                Heat capacity

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        R_adj = _get_R_adj(units=units, elements=self.elements)
        return self.get_CpoR(T=T,
                             raise_error=raise_error,
                             raise_warning=raise_warning,
                             **kwargs) * R_adj

    def get_HoRT(self, T, raise_error=True, raise_warning=True, **kwargs):
        """Calculate the dimensionless enthalpy

        Parameters
        ----------
            T : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Arguments to calculate mixture model properties, if any
        Returns
        -------
            HoRT : float or (N,) `numpy.ndarray`_
                Dimensionless enthalpy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        if isinstance(T, np.ndarray):
            T = T.astype(float)
        else:
            T = float(T)
        # Convert T to 1D numpy format
        if not _is_iterable(T):
            T = [T]
        T = np.array(T)
        self._check_T(T)

        # Calculate pure properties
        HoRT = get_shomate_HoRT(a=self.a, T=T, units=self.units)
        # Calculate mixing properties
        for T_i in T:
            HoRT_mix = _get_mix_quantity(misc_models=self.misc_models,
                                         method_name='get_HoRT',
                                         raise_error=raise_error,
                                         raise_warning=raise_warning,
                                         default_value=0.,
                                         T=T_i,
                                         **kwargs)
        # Add mixing quantity in appropriate format
        HoRT = HoRT + HoRT_mix
        if len(T) == 1:
            HoRT = HoRT.item(0)
        return HoRT

    def get_H(self, T, units, raise_error=True, raise_warning=True, **kwargs):
        """Calculate the enthalpy

        Parameters
        ----------
            T : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Arguments to calculate mixture model properties, if any
        Returns
        -------
            H : float or (N,) `numpy.ndarray`_
                Enthalpy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        units = '{}/K'.format(units)
        R_adj = _get_R_adj(units=units, elements=self.elements)
        return self.get_HoRT(T=T,
                             raise_error=raise_error,
                             raise_warning=raise_warning,
                             **kwargs) * T * R_adj

    def get_Selements(self):
        """Calculate the dimensionless entropy of the elements in the molecule

        Parameters
        ----------
        None

        Returns
        -------
        SoR : float
              Entropy
        """
        elements = self.elements
        S_ele = 0
        for element in elements:
            S_ele += c.S_elements[element]*elements[element]
        return S_ele

    def get_SoR(self, T, raise_error=True, raise_warning=True,
                S_elements=None, **kwargs):
        """Calculate the dimensionless entropy

        Parameters
        ----------
            T : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            S_elements : bool, optional
                Includes the entropy of the elements to compute an entropy of
                formation. Defauly is None
            kwargs : key-word arguments
                Arguments to calculate mixture model properties, if any
        Returns
        -------
            SoR : float or (N,) `numpy.ndarray`_
                Dimensionless entropy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        if isinstance(T, np.ndarray):
            T = T.astype(float)
        else:
            T = float(T)
        # Convert T to 1D numpy format
        if not _is_iterable(T):
            T = [T]
        T = np.array(T)
        self._check_T(T)

        # Calculate pure properties
        SoR = get_shomate_SoR(a=self.a, T=T, units=self.units)
        # Calculate mixing properties
        for T_i in T:
            SoR_mix = _get_mix_quantity(misc_models=self.misc_models,
                                        method_name='get_SoR',
                                        raise_error=raise_error,
                                        raise_warning=raise_warning,
                                        default_value=0.,
                                        T=T_i,
                                        **kwargs)
        SoR = SoR + SoR_mix
        # If only one T specified, converted from (1,) numpy array to float
        if not S_elements:
            S_ele = 0
        else:
            S_ele = self.get_Selements()

        if len(T) == 1:
            SoR = SoR.item(0)
        return SoR - S_ele

    def get_S(self, T, units, raise_error=True, raise_warning=True,
              S_elements=None, **kwargs):
        """Calculate the entropy

        Parameters
        ----------
            T : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
                units.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            S_elements : bool, optional
                Includes the entropy of the elements to compute an entropy of
                formation. Defauly is None
            kwargs : key-word arguments
                Arguments to calculate mixture model properties, if any
        Returns
        -------
            S : float or (N,) `numpy.ndarray`_
                Entropy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        R_adj = _get_R_adj(units=units, elements=self.elements)
        return self.get_SoR(T=T, S_elements=S_elements) * R_adj

    def get_GoRT(self, T, raise_error=True, raise_warning=True,
                 S_elements=None, **kwargs):
        """Calculate the dimensionless Gibbs free energy

        Parameters
        ----------
            T : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            S_elements : bool, optional
                Includes the entropy of the elements to compute an entropy of
                formation. Defauly is None
            kwargs : key-word arguments
                Arguments to calculate mixture model properties, if any
        Returns
        -------
            GoRT : float or (N,) `numpy.ndarray`_
                Dimensionless Gibbs free energy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        return self.get_HoRT(T=T, raise_error=raise_error,
                             raise_warning=raise_warning, **kwargs) \
            - self.get_SoR(T=T, raise_error=raise_error,
                           raise_warning=raise_warning,
                           S_elements=S_elements, **kwargs)

    def get_G(self, T, units, raise_error=True, raise_warning=True,
              S_elements=None, **kwargs):
        """Calculate the Gibbs energy

        Parameters
        ----------
            T : float or (N,) `numpy.ndarray`_
                Temperature(s) in K
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            S_elements : bool, optional
                Includes the entropy of the elements to compute an entropy of
                formation. Defauly is None
            kwargs : key-word arguments
                Arguments to calculate mixture model properties, if any
        Returns
        -------
            G : float or (N,) `numpy.ndarray`_
                Gibbs energy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        units = '{}/K'.format(units)
        R_adj = _get_R_adj(units=units, elements=self.elements)
        return self.get_GoRT(T=T,
                             raise_error=raise_error,
                             raise_warning=raise_warning,
                             S_elements=S_elements,
                             **kwargs) * T * R_adj

    @classmethod
    def from_data(cls,
                  name,
                  T,
                  CpoR,
                  T_ref,
                  HoRT_ref,
                  SoR_ref,
                  units='J/mol/K',
                  **kwargs):
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
            units : str, optional
                Units used to fit the Shomate polynomial. Units should be
                supported by :class:`~pmutt.constants.R` (e.g. J/mol/K,
                cal/mol/K, eV/K). Default is J/mol/K.
        Returns
        -------
            shomate : Shomate object
                Shomate object with polynomial terms fitted to data.

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        T_low = min(T)
        T_high = max(T)
        a = _fit_CpoR(T=T, CpoR=CpoR, units=units)
        a = _fit_HoRT(T_ref=T_ref, HoRT_ref=HoRT_ref, a=a, units=units)
        a = _fit_SoR(T_ref=T_ref, SoR_ref=SoR_ref, a=a, units=units)
        return cls(name=name,
                   T_low=T_low,
                   T_high=T_high,
                   a=a,
                   units=units,
                   **kwargs)

    @classmethod
    def from_statmech(cls,
                      name,
                      statmech_model,
                      T_low,
                      T_high,
                      references=None,
                      elements=None,
                      **kwargs):
        """Calculates the Shomate polynomial using statistical mechanic models.
        Deprecated as of Version 1.2.13. Please use ``from_model`` instead.

        Parameters
        ----------
            name : str
                Name of the species
            statmech_model : `pmutt.statmech.StatMech` object or class
                Statistical Mechanics model to generate data
            T_low : float
                Lower limit temerature in K
            T_high : float
                Higher limit temperature in K
            references : `pmutt.empirical.references.References` object
                Reference to adjust enthalpy
            **kwargs : keyword arguments
                Used to initalize ``statmech_model`` or ``EmpiricalBase``
                attributes to be stored.
        Returns
        -------
            shomate : Shomate object
                Shomate object with polynomial terms fitted to data.
        """
        err_msg = ('Shomate.from_statmech is deprecated as of Version 1.2.13. '
                   'Please use the more generic function, Shomate.from_model.')
        raise RuntimeError(err_msg)

    @classmethod
    def from_model(cls,
                   model,
                   name=None,
                   T_low=None,
                   T_high=None,
                   elements=None,
                   n_T=50,
                   units='J/mol/K',
                   **kwargs):
        """Calculates the NASA polynomials using the model passed

        Parameters
        ----------
            model : Model object or class
                Model to generate data. Must contain the methods `get_CpoR`,
                `get_HoRT` and `get_SoR`
            name : str, optional
                Name of the species. If not passed, `model.name` will be used.
            T_low : float, optional
                Lower limit temerature in K. If not passed, `model.T_low` will
                be used.
            T_high : float, optional
                Higher limit temperature in K. If not passed, `model.T_high`
                will be used.
            elements : dict, optional
                Composition of the species. If not passed, `model.elements`
                will be used. Keys of dictionary are elements, values are
                stoichiometric values in a formula unit.
                e.g. CH3OH can be represented as:
                {'C': 1, 'H': 4, 'O': 1,}.
            n_T : int, optional
                Number of data points between `T_low` and `T_high` for fitting
                heat capacity. Default is 50.
            units : str, optional
                Units used to fit the Shomate polynomial. Units should be
                supported by :class:`~pmutt.constants.R` (e.g. J/mol/K,
                cal/mol/K, eV/K). Default is J/mol/K.
            kwargs : keyword arguments
                Used to initalize model if a class is passed.
        Returns
        -------
            Shomate : Shomate object
                Shomate object with polynomial terms fitted to data.
        """
        # Initialize the model object
        if inspect.isclass(model):
            model = model(name=name, elements=elements, **kwargs)

        if name is None:
            try:
                name = model.name
            except AttributeError:
                err_msg = ('Name must either be passed to from_model directly '
                           'or be an attribute of model.')
                raise AttributeError(err_msg)
        if T_low is None:
            try:
                T_low = model.T_low
            except AttributeError:
                err_msg = ('T_low must either be passed to from_model '
                           'directly or be an attribute of model.')
                raise AttributeError(err_msg)
        if T_high is None:
            try:
                T_high = model.T_high
            except AttributeError:
                err_msg = ('T_high must either be passed to from_model '
                           'directly or be an attribute of model.')
                raise AttributeError(err_msg)
        if elements is None:
            try:
                elements = model.elements
            except AttributeError:
                pass
        # Check if inputted T_low and T_high are outside model's T_low and
        # T_high range
        try:
            model.T_low
        except AttributeError:
            pass
        else:
            if T_low < model.T_low:
                warn_msg = ('Inputted T_low ({} K) is lower than model '
                            'T_low ({} K). Fitted empirical object may not be '
                            'valid.'
                            ''.format(T_low, model.T_low))
                warn(warn_msg, UserWarning)

        try:
            if T_high > model.T_high:
                warn_msg = ('Inputted T_high ({} K) is higher than model '
                            'T_high ({} K). Fitted empirical object may not '
                            'be valid.'
                            ''.format(T_high, model.T_high))
                warn(warn_msg, UserWarning)
        except AttributeError:
            pass

        # Generate heat capacity data
        T = np.linspace(T_low, T_high, n_T)
        try:
            CpoR = model.get_CpoR(T=T)
        except ValueError:
            CpoR = np.array([model.get_CpoR(T=T_i) for T_i in T])
        else:
            if not _is_iterable(CpoR) or len(CpoR) != len(T):
                CpoR = np.array([model.get_CpoR(T=T_i) for T_i in T])

        # Generate enthalpy and entropy data
        T_mean = (T_low + T_high) / 2.
        HoRT_ref = model.get_HoRT(T=T_mean)
        SoR_ref = model.get_SoR(T=T_mean)
        return cls.from_data(name=name,
                             T=T,
                             CpoR=CpoR,
                             T_ref=T_mean,
                             HoRT_ref=HoRT_ref,
                             SoR_ref=SoR_ref,
                             model=model,
                             elements=elements,
                             units=units,
                             **kwargs)

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = super().to_dict()
        obj_dict['class'] = str(self.__class__)
        obj_dict['type'] = 'shomate'
        obj_dict['a'] = list(self.a)
        obj_dict['T_low'] = self.T_low
        obj_dict['T_high'] = self.T_high
        obj_dict['units'] = self.units
        return obj_dict

    def to_omkm_yaml(self):
        """Returns a dictionary compatible with Cantera's YAML format

        Returns
        -------
            yaml_dict : dict
                Dictionary compatible with Cantera's YAML format
        """
        yaml_dict = {
            'name': self.name,
            'composition': self.elements,
            'thermo': {'model': 'Shomate',
                       'temperature-ranges': [float(self.T_low),
                                              float(self.T_high)],
                       'data': [self.a.tolist()[:-1]]}
        }
        if self.n_sites is not None:
            yaml_dict['sites'] = self.n_sites,
        return yaml_dict

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
        json_obj['model'] = \
            json_to_pmutt(json_obj['model'])
        json_obj['misc_models'] = json_to_pmutt(json_obj['misc_models'])

        return cls(**json_obj)

    def to_cti(self):
        """Writes the object in Cantera's CTI format.

        Returns
        -------
            CTI_str : str
                Object represented as a CTI string.
        """
        if self.n_sites is None:
            size_str = ''
        else:
            size_str = ' size={},'.format(self.n_sites)
        cti_str = ('species(name="{}", atoms={},{}\n'
                   '        thermo=Shomate([{}, {}],\n'
                   '                       [{: 2.8E}, {: 2.8E}, {: 2.8E},\n'
                   '                        {: 2.8E}, {: 2.8E}, {: 2.8E},\n'
                   '                        {: 2.8E}]))').format(
                       self.name, obj_to_cti(self.elements), size_str,
                       self.T_low, self.T_high, self.a[0], self.a[1],
                       self.a[2], self.a[3], self.a[4], self.a[5], self.a[6])
        return cti_str

    def _check_T(self, T):
        for T_i in T:
            if T_i < self.T_low:
                warn_msg = ('Requested temperature ({} K), below T_low ({} K)'
                            'for Shomate object, {}'
                            ''.format(T, self.T_low, self.name))
                warn(warn_msg, RuntimeWarning)
            elif T_i > self.T_high:
                warn_msg = ('Requested temperature ({} K), above T_high ({} K)'
                            'for Shomate object, {}'
                            ''.format(T, self.T_high, self.name))
                warn(warn_msg, RuntimeWarning)


def _fit_CpoR(T, CpoR, units):
    """Fit a[0]-a[4] coefficients given the dimensionless heat capacity data

    Parameters
    ----------
        T : (N,) `numpy.ndarray`_
            Temperatures in K
        CpoR : (N,) `numpy.ndarray`_
            Dimensionless heat capacity
        units : str
            Units corresponding to Shomate polynomial. Units should be
            supported by :class:`~pmutt.constants.R`.
    Returns
    -------
        a : (8,) `numpy.ndarray`_
            Lower coefficients of Shomate polynomial

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
    """
    # If the Cp/R does not vary with temperature (occurs when no
    # vibrational frequencies are listed), return default values
    if all([np.isclose(x, 0.) for x in CpoR]) \
       or any([np.isnan(x) for x in CpoR]):
        return np.zeros(7)
    else:
        # Pass the unit set
        adj_shomate_CpoR = lambda T, A, B, C, D, E: _shomate_CpoR(
            T=T, A=A, B=B, C=C, D=D, E=E, units=units)
        [a, _] = curve_fit(adj_shomate_CpoR, T, np.array(CpoR))
        a = np.append(a, [0., 0., 0.])
        return a


def _fit_HoRT(T_ref, HoRT_ref, a, units):
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
        units : str
            Units corresponding to Shomate polynomial. Units should be
            supported by :class:`~pmutt.constants.R`.
    Returns
    -------
        a : (8,) `numpy.ndarray`_
            Lower coefficients of Shomate polynomial

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
    """
    a[5] = (HoRT_ref
            - get_shomate_HoRT(T=np.array([T_ref]), a=a, units=units))[0] \
        * c.R(units)*T_ref/c.prefixes['k']
    a[7] = - get_shomate_HoRT(T=np.array([c.T0('K')]), a=a, units=units)[0] \
        * c.R(units)*c.T0('K')/c.prefixes['k']
    return a


def _fit_SoR(T_ref, SoR_ref, a, units):
    """Fit a[6] coefficient in a_low and a_high attributes given the
    dimensionless entropy

    Parameters
    ----------
        T_ref : float
            Reference temperature in K
        SoR_ref : float
            Reference dimensionless entropy
        units : str
            Units corresponding to Shomate polynomial. Units should be
            supported by :class:`~pmutt.constants.R`.
    Returns
    -------
        a : (8,) `numpy.ndarray`_
            Lower coefficients of Shomate polynomial

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
    """
    a[6] = c.R(units) * (
        SoR_ref - get_shomate_SoR(T=np.array([T_ref]), a=a, units=units))[0]
    return a


def get_shomate_CpoR(a, T, units):
    """Calculates the dimensionless heat capacity using Shomate polynomial form

    Parameters
    ----------
        a : (8,) `numpy.ndarray`_
            Coefficients of Shomate polynomial
        T : iterable
            Temperature in K
        units : str
            Units corresponding to Shomate polynomial. Units should be
            supported by :class:`~pmutt.constants.R`.
    Returns
    -------
        CpoR: float
            Dimensionless heat capacity

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
    """
    t = T / 1000.
    t_arr = np.array([[1., x, x**2, x**3, 1. / x**2, 0., 0., 0.] for x in t])
    return np.dot(t_arr, a) / c.R(units)


def get_shomate_HoRT(a, T, units):
    """Calculates the dimensionless enthalpy using Shomate polynomial form

    Parameters
    ----------
        a : (8,) `numpy.ndarray`_
            Coefficients of Shomate polynomial
        T : iterable
            Temperature in K
        units : str
            Units corresponding to Shomate polynomial. Units should be
            supported by :class:`~pmutt.constants.R`.
    Returns
    -------
        HoRT : float
            Dimensionless enthalpy

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
    """
    t = T / 1000.
    t_arr = np.array(
        [[x, x**2 / 2., x**3 / 3., x**4 / 4., -1. / x, 1., 0., 0.] for x in t])
    HoRT = np.dot(t_arr, a) / (T * c.R(units) / c.prefixes['k'])
    return HoRT


def get_shomate_SoR(a, T, units):
    """Calculates the dimensionless entropy using Shomate polynomial form

    Parameters
    ----------
        a : (8,) `numpy.ndarray`_
            Coefficients of Shomate polynomial
        T : iterable
            Temperature in K
        units : str
            Units corresponding to Shomate polynomial. Units should be
            supported by :class:`~pmutt.constants.R`.
    Returns
    -------
        SoR : float
            Dimensionless entropy

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
    """
    t = T / 1000.
    t_arr = np.array(
        [[np.log(x), x, x**2 / 2., x**3 / 3., -1. / 2. / x**2, 0., 1., 0.]
         for x in t])
    SoR = np.dot(t_arr, a) / c.R(units)
    return SoR


def get_shomate_GoRT(a, T, units):
    """Calculates the dimensionless Gibbs free energy using Shomate
    polynomial form

    Parameters
    ----------
        a : (8,) `numpy.ndarray`_
            Coefficients of Shomate polynomial
        T : iterable
            Temperature in K
        units : str
            Units corresponding to Shomate polynomial. Units should be
            supported by :class:`~pmutt.constants.R`.
    Returns
    -------
        GoRT : float
            Dimensionless Gibbs energy

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
    """
    GoRT = get_shomate_HoRT(a=a, T=T, units=units) \
        - get_shomate_SoR(a=a, T=T, units=units)
    return GoRT


def _shomate_CpoR(T, A, B, C, D, E, units):
    """
    Helper function to fit Shomate heat capacity.

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
    if not _is_iterable(T):
        T = [T]
    T = np.array(T)
    return get_shomate_CpoR(a=a, T=T, units=units)
