# -*- coding: utf-8 -*-
"""
PyMuTT.models.empirical

Empirical models.
"""

import inspect
from matplotlib import pyplot as plt
import numpy as np
from PyMuTT import _pass_expected_arguments
from PyMuTT import constants as c
from PyMuTT.io_.jsonio import json_to_PyMuTT

class BaseThermo:
    """The Thermodynamic Parent class.
    Holds properties of a specie, the statistical-mechanical thermodynamic
    model.

    Attributes
    ----------
        name : str
            Name of the specie.
        phase : str
            Phase of the specie.
            G - gas.
            S - surface.
        elements : dict
            Composition of the species.
            Keys of dictionary are elements, values are stoichiometric values
            in a formula unit.
            e.g. CH3OH can be represented as:
            {'C': 1, 'H': 4, 'O': 1,}.
        statmech_model : `PyMuTT.models.statmech` object
            Statistical thermodynamic model.
            Object should have the following methods: `get_CpoR`, `get_HoRT`,
            `get_SoR`, `get_GoRT`.
        T_ref : float
            Temperature (in K) at which `HoRT_dft` was calculated. Only used
            for fitting empirical coefficients.
        HoRT_dft : float
            Dimensionless enthalpy calculated using DFT that corresponds to
            `T_ref`. Only used for fitting empirical coefficients.
        HoRT_ref : float
            Reference dimensionless enthalpy corresponding to `T_ref`.
        references : `PyMuTT.models.empirical.References.references` object
            Contains references to calculate `HoRT_ref`. If not specified then
            HoRT_dft will be used without adjustment.
        notes : str
            Any additional details you would like to include such as
            computational set up.

    """

    def __init__(self, name, phase=None, elements=None, statmech_model=None,
                 T_ref=c.T0('K'), HoRT_dft=None, HoRT_ref=None,
                 references=None, notes=None, **kwargs):
        self.name = name
        self.phase = phase
        self.elements = elements
        self.T_ref = T_ref
        self.references = references
        self.notes = notes

        # Assign self.statmech_model
        if inspect.isclass(statmech_model):
            # If you're passing a class. Note that the required
            # arguments will be guessed.
            self.statmech_model = statmech_model(**kwargs)
        else:
            # If it's an object that has already been initialized
            self.statmech_model = statmech_model

        # Calculate dimensionless DFT energy using thermo model
        if (HoRT_dft is None) and (self.statmech_model is not None):
            self.HoRT_dft = self.statmech_model.get_HoRT(T=self.T_ref)
        else:
            self.HoRT_dft = HoRT_dft

        # Assign self.HoRT_ref
        if HoRT_ref is None:
            if (references is None) or (self.HoRT_dft is None):
                self.HoRT_ref = self.HoRT_dft
            else:
                self.HoRT_ref = self.HoRT_dft +\
                    references.get_HoRT_offset(elements=elements,
                                               Ts=self.T_ref)
        else:
            self.HoRT_ref = HoRT_ref

    def __repr__(self):
        out = ['{} object for Name: {}'.format(self.__class__.__name__,
               self.name)]
        for key, val in self.__dict__.items():
            if key != 'name':
                out.append('\t{}: {}'.format(key, val))
        return '\n'.join(out)

    def plot_empirical(self, T_low=None, T_high=None, Cp_units=None,
                       H_units=None, S_units=None, G_units=None):
        """Plots the thermodynamic profiles between ``T_low`` and ``T_high``
        using empirical relationship

        Parameters
        ----------
            T_low : float
                Lower temperature in K. If not specified,
                ``T_low`` attribute used.
            T_high : float
                Upper temperature in K. If not specified,
                ``T_high`` attribute used.
            Cp_units : str
                Units to plot heat capacity. See ``PyMuTT.constants.R``
                for accepted units. If not specified, dimensionless units used.
            H_units : str
                Units to plot enthalpy. See ``PyMuTT.constants.R`` for accepted
                units but omit the '/K' (e.g. J/mol). If not specified,
                dimensionless units used.
            S_units : str
                Units to plot entropy. See ``PyMuTT.constants.R`` for accepted
                units. If not specified, dimensionless units used.
            G_units : str
                Units to plot Gibbs free energy. See ``PyMuTT.constants.R``
                for accepted units but omit the '/K' (e.g. J/mol). If not
                specified, dimensionless units used.
        Returns
        -------
            figure : `matplotlib.figure.Figure`_
                Figure
            axes : tuple of `matplotlib.axes.Axes.axis`_
                Axes of the plots.
                0. Cp
                1. H
                2. S
                3. G

        .. _`matplotlib.figure.Figure`: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html
        .. _`matplotlib.axes.Axes.axis`: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axis.html
        """
        if T_low is None:
            T_low = self.T_low
        if T_high is None:
            T_high = self.T_high
        Ts = np.linspace(T_low, T_high)

        f, ax = plt.subplots(4, sharex=True)
        '''
        Heat Capacity
        '''
        ax[0].set_title('Specie: {}'.format(self.name))
        Cp_plot = self.get_CpoR(Ts=Ts)
        if Cp_units is None:
            ax[0].set_ylabel('Cp/R')
        else:
            ax[0].set_ylabel('Cp ({})'.format(Cp_units))
            Cp_plot = Cp_plot * c.R(Cp_units)
        ax[0].plot(Ts, Cp_plot, 'r-')

        '''
        Enthalpy
        '''
        H_plot = self.get_HoRT(Ts=Ts)
        if H_units is None:
            ax[1].set_ylabel('H/RT')
        else:
            ax[1].set_ylabel('H ({})'.format(H_units))
            H_plot = H_plot * c.R('{}/K'.format(H_units)) * Ts
        ax[1].plot(Ts, H_plot, 'g-')

        '''
        Entropy
        '''
        S_plot = self.get_SoR(Ts=Ts)
        if S_units is None:
            ax[2].set_ylabel('S/R')
        else:
            ax[2].set_ylabel('S ({})'.format(S_units))
            S_plot = S_plot * c.R(S_units)
        ax[2].plot(Ts, S_plot, 'b-')

        '''
        Gibbs energy
        '''
        ax[3].set_xlabel('Temperature (K)')
        G_plot = self.get_GoRT(Ts=Ts)
        if G_units is None:
            ax[3].set_ylabel('G/RT')
        else:
            ax[3].set_ylabel('G ({})'.format(G_units))
            G_plot = G_plot * c.R('{}/K'.format(G_units)) * Ts
        ax[3].plot(Ts, G_plot, 'k-')

        return f, ax

    def plot_statmech(self, T_low=None, T_high=None, Cp_units=None,
                          H_units=None, S_units=None, G_units=None):
        """Plots the thermodynamic profiles between ``T_low`` and ``T_high``
        using empirical relationship

        Parameters
        ----------
            T_low : float
                Lower temperature in K. If not specified,
                ``T_low`` attribute used
            T_high : float
                Upper temperature in K. If not specified,
                ``T_high`` attribute used
            Cp_units : str
                Units to plot heat capacity. See ``PyMuTT.constants.R``
                for accepted units. If not specified, dimensionless units used.
            H_units : str
                Units to plot enthalpy. See ``PyMuTT.constants.R`` for
                accepted units but omit the '/K' (e.g. J/mol). If not
                specified, dimensionless units used.
            S_units : str
                Units to plot entropy. See ``PyMuTT.constants.R`` for
                accepted units. If not specified, dimensionless units used.
            G_units : str
                Units to plot Gibbs free energy. See ``PyMuTT.constants.R``
                for accepted units but omit the '/K' (e.g. J/mol). If not
                specified, dimensionless units used.
        Returns
        -------
            figure : `matplotlib.figure.Figure`_
                Figure
            axes : tuple of `matplotlib.axes.Axes.axis`_
                Axes of the plots.
                0. Cp
                1. H
                2. S
                3. G

        .. _`matplotlib.figure.Figure`: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html
        .. _`matplotlib.axes.Axes.axis`: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axis.html
        """
        if T_low is None:
            T_low = self.T_low
        if T_high is None:
            T_high = self.T_high
        Ts = np.linspace(T_low, T_high)

        f, ax = plt.subplots(4, sharex=True)
        '''
        Heat Capacity
        '''
        ax[0].set_title('Specie: {}'.format(self.name))
        Cp_plot = [self.statmech_model.get_CpoR(T=T) for T in Ts]
        if Cp_units is None:
            ax[0].set_ylabel('Cp/R')
        else:
            ax[0].set_ylabel('Cp ({})'.format(Cp_units))
            Cp_plot = Cp_plot * c.R(Cp_units)
        ax[0].plot(Ts, Cp_plot, 'r-')

        '''
        Enthalpy
        '''
        H_plot = [self.statmech_model.get_HoRT(T=T) for T in Ts]
        if self.references is not None:
            H_plot += self.references.get_HoRT_offset(elements=self.elements,
                                                      Ts=Ts)

        if H_units is None:
            ax[1].set_ylabel('H/RT')
        else:
            ax[1].set_ylabel('H ({})'.format(H_units))
            H_plot = H_plot * c.R('{}/K'.format(H_units)) * Ts
        ax[1].plot(Ts, H_plot, 'g-')

        '''
        Entropy
        '''
        S_plot = [self.statmech_model.get_SoR(T=T) for T in Ts]
        if S_units is None:
            ax[2].set_ylabel('S/R')
        else:
            ax[2].set_ylabel('S ({})'.format(S_units))
            S_plot = S_plot * c.R(S_units)
        ax[2].plot(Ts, S_plot, 'b-')

        '''
        Gibbs energy
        '''
        ax[3].set_xlabel('Temperature (K)')
        G_plot = [self.statmech_model.get_GoRT(T=T) for T in Ts]
        if self.references is not None:
            G_plot += self.references.get_HoRT_offset(elements=self.elements,
                                                      Ts=Ts)

        if G_units is None:
            ax[3].set_ylabel('G/RT')
        else:
            ax[3].set_ylabel('G ({})'.format(G_units))
            G_plot = G_plot * c.R('{}/K'.format(G_units)) * Ts
        ax[3].plot(Ts, G_plot, 'k-')

        return f, ax

    def plot_statmech_and_empirical(self, T_low=None, T_high=None,
                                        Cp_units=None, H_units=None,
                                        S_units=None, G_units=None):
        """Plots the thermodynamic profiles between ``T_low`` and ``T_high``
        using empirical relationship

        Parameters
        ----------
            T_low : float
                Lower temperature in K. If not specified, ``T_low``
                attribute used
            T_high : float
                Upper temperature in K. If not specified, ``T_high``
                attribute used
            Cp_units : str
                Units to plot heat capacity. See ``PyMuTT.constants.R`` for
                accepted units. If not specified, dimensionless units used.
            H_units : str
                Units to plot enthalpy. See ``PyMuTT.constants.R`` for accepted
                units but omit the '/K' (e.g. J/mol). If not specified,
                dimensionless units used.
            S_units : str
                Units to plot entropy. See ``PyMuTT.constants.R`` for accepted
                units. If not specified, dimensionless units used.
            G_units : str
                Units to plot Gibbs free energy. See ``PyMuTT.constants.R``
                for accepted units but omit the '/K' (e.g. J/mol). If not
                specified, dimensionless units used.
        Returns
        -------
            figure : `matplotlib.figure.Figure`_
                Figure
            axes : tuple of `matplotlib.axes.Axes.axis`_
                Axes of the plots.
                0. Cp
                1. H
                2. S
                3. G

        .. _`matplotlib.figure.Figure`: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html
        .. _`matplotlib.axes.Axes.axis`: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axis.html
        """
        if T_low is None:
            T_low = self.T_low
        if T_high is None:
            T_high = self.T_high
        Ts = np.linspace(T_low, T_high)

        f, ax = plt.subplots(4, sharex=True)
        '''
        Heat Capacity
        '''
        ax[0].set_title('Specie: {}'.format(self.name))
        Ts, Cp_plot_statmech, Cp_plot_empirical = self.compare_CpoR(Ts=Ts)
        if Cp_units is None:
            ax[0].set_ylabel('Cp/R')
        else:
            ax[0].set_ylabel('Cp ({})'.format(Cp_units))
            Cp_plot_statmech = Cp_plot_statmech * c.R(Cp_units)
            Cp_plot_empirical = Cp_plot_empirical * c.R(Cp_units)

        ax[0].plot(Ts, Cp_plot_statmech, 'r-', label='Stat Mech Model')
        ax[0].plot(Ts, Cp_plot_empirical, 'b-', label='Empirical Model')
        ax[0].legend()

        '''
        Enthalpy
        '''
        Ts, H_plot_statmech, H_plot_empirical = self.compare_HoRT(Ts=Ts)

        if H_units is None:
            ax[1].set_ylabel('H/RT')
        else:
            ax[1].set_ylabel('H ({})'.format(H_units))
            H_plot_statmech = H_plot_statmech *\
                c.R('{}/K'.format(H_units)) * Ts
            H_plot_empirical = H_plot_empirical *\
                c.R('{}/K'.format(H_units)) * Ts
        ax[1].plot(Ts, H_plot_statmech, 'r-')
        ax[1].plot(Ts, H_plot_empirical, 'b-')

        '''
        Entropy
        '''
        Ts, S_plot_statmech, S_plot_empirical = self.compare_SoR(Ts=Ts)
        if S_units is None:
            ax[2].set_ylabel('S/R')
        else:
            ax[2].set_ylabel('S ({})'.format(S_units))
            S_plot_statmech = S_plot_statmech * c.R(S_units)
            S_plot_empirical = S_plot_empirical * c.R(S_units)
        ax[2].plot(Ts, S_plot_statmech, 'r-')
        ax[2].plot(Ts, S_plot_empirical, 'b-')

        '''
        Gibbs energy
        '''
        ax[3].set_xlabel('Temperature (K)')
        Ts, G_plot_statmech, G_plot_empirical = self.compare_GoRT(Ts=Ts)
        if G_units is None:
            ax[3].set_ylabel('G/RT')
        else:
            ax[3].set_ylabel('G ({})'.format(G_units))
            G_plot_statmech = G_plot_statmech *\
                c.R('{}/K'.format(G_units)) * Ts
            G_plot_empirical = G_plot_empirical *\
                c.R('{}/K'.format(G_units)) * Ts
        ax[3].plot(Ts, G_plot_statmech, 'r-')
        ax[3].plot(Ts, G_plot_empirical, 'b-')

        return f, ax

    def compare_CpoR(self, Ts=None):
        """Compares the dimensionless heat capacity of the statistical model
        and the empirical model

        Parameters
        ----------
            Ts : (N,) `numpy.ndarray`_ or float, optional
                Temperatures (in K) to calculate CpoR. If None, generates a
                list of temperatures between self.T_low and self.T_high
        Returns
        -------
            Ts : (N,) `numpy.ndarray`_ or float
                Temperatures in K
            CpoR_statmech : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of statistical thermodynamic model
            CpoR_empirical :((N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of empirical model

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """

        if Ts is None:
            Ts = np.linspace(self.T_low, self.T_high)

        try:
            iter(Ts)
        except TypeError:
            CpoR_statmech = [self.statmech_model.get_CpoR(T=T) for T in Ts]
            CpoR_empirical = self.get_CpoR(Ts=Ts)
        else:
            CpoR_statmech = np.zeros_like(Ts)
            CpoR_empirical = np.zeros_like(Ts)
            for i, T in enumerate(Ts):
                CpoR_statmech[i] = self.statmech_model.get_CpoR(T=T)
                CpoR_empirical[i] = self.get_CpoR(Ts=T)
        return (Ts, CpoR_statmech, CpoR_empirical)

    def compare_HoRT(self, Ts=None):
        """Compares the dimensionless enthalpy of the statistical model and
        the empirical model

        Parameters
        ----------
            Ts : (N,) `numpy.ndarray`_ or float, optional
                Temperatures (in K) to calculate CpoR. If None, generates a
                list of temperatures between self.T_low and self.T_high
        Returns
        -------
            Ts : (N,) `numpy.ndarray`_ or float
                Temperatures in K
            CpoR_statmech : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of statistical thermodynamic model
            CpoR_empirical :((N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of empirical model

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        if Ts is None:
            Ts = np.linspace(self.T_low, self.T_high)

        if self.references is not None:
            H_offset = self.references.get_HoRT_offset(elements=self.elements,
                                                       Ts=Ts)

        try:
            iter(Ts)
        except TypeError:
            HoRT_statmech = [self.statmech_model.get_HoRT(T=T) + H_offset \
                             for T in Ts]
            HoRT_empirical = self.get_HoRT(Ts=Ts)
        else:
            HoRT_statmech = np.zeros_like(Ts)
            HoRT_empirical = np.zeros_like(Ts)
            for i, T in enumerate(Ts):
                HoRT_statmech[i] = self.statmech_model.get_HoRT(T=T) \
                                  + H_offset[i]
                HoRT_empirical[i] = self.get_HoRT(Ts=T)
        return (Ts, HoRT_statmech, HoRT_empirical)

    def compare_SoR(self, Ts=None):
        """Compares the dimensionless entropy of the statistical model and the
        empirical model

        Parameters
        ----------
            Ts : (N,) `numpy.ndarray`_ or float, optional
                Temperatures (in K) to calculate CpoR. If None, generates a
                list of temperatures between self.T_low and self.T_high
        Returns
        -------
            Ts : (N,) `numpy.ndarray`_ or float
                Temperatures in K
            CpoR_statmech : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of statistical thermodynamic model
            CpoR_empirical :((N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of empirical model

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        if Ts is None:
            Ts = np.linspace(self.T_low, self.T_high)

        try:
            iter(Ts)
        except TypeError:
            SoR_statmech = [self.statmech_model.get_SoR(T=T) for T in Ts]
            SoR_empirical = self.get_SoR(Ts=Ts)
        else:
            SoR_statmech = np.zeros_like(Ts)
            SoR_empirical = np.zeros_like(Ts)
            for i, T in enumerate(Ts):
                SoR_statmech[i] = self.statmech_model.get_SoR(T=T)
                SoR_empirical[i] = self.get_SoR(Ts=T)
        return (Ts, SoR_statmech, SoR_empirical)

    def compare_GoRT(self, Ts=None):
        """Compares the dimensionless Gibbs energy of the statistical model and
        the empirical model

        Parameters
        ----------
            Ts : (N,) `numpy.ndarray`_ or float, optional
                Temperatures (in K) to calculate CpoR. If None, generates a
                list of temperatures between self.T_low and self.T_high
        Returns
        -------
            Ts : (N,) `numpy.ndarray`_ or float
                Temperatures in K
            CpoR_statmech : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of statistical thermodynamic model
            CpoR_empirical : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of empirical model

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        if Ts is None:
            Ts = np.linspace(self.T_low, self.T_high)

        if self.references is not None:
            G_offset = self.references.get_HoRT_offset(elements=self.elements,
                                                       Ts=Ts)

        try:
            iter(Ts)
        except TypeError:
            GoRT_statmech = [self.statmech_model.get_GoRT(T=T) + G_offset for T in Ts]
            GoRT_empirical = self.get_GoRT(Ts=Ts)
        else:
            GoRT_statmech = np.zeros_like(Ts)
            GoRT_empirical = np.zeros_like(Ts)
            for i, T in enumerate(Ts):
                GoRT_statmech[i] = self.statmech_model.get_GoRT(T=T) +\
                    G_offset[i]
                GoRT_empirical[i] = self.get_GoRT(Ts=T)
        return (Ts, GoRT_statmech, GoRT_empirical)

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes
        
        Returns
        -------
            obj_dict : dict
        """
        obj_dict = {'class': str(self.__class__),
                    'name': self.name,
                    'phase': self.phase,
                    'elements': self.elements,
                    'T_ref': self.T_ref,
                    'notes': self.notes,
                    'HoRT_dft': self.HoRT_dft,
                    'HoRT_ref': self.HoRT_ref}
        try:
            obj_dict['references'] = self.references.to_dict()
        except AttributeError:
            obj_dict['references'] = self.references

        try:
            obj_dict['statmech_model'] = self.statmech_model.to_dict()
        except AttributeError:
            obj_dict['statmech_model'] = self.statmech_model

        return obj_dict

    @classmethod
    def from_dict(cls, json_obj):
        try:
            del json_obj['class']
        except KeyError:
            pass

        # Reconstruct statmech model
        json_obj['statmech_model'] = \
                json_to_PyMuTT(json_obj['statmech_model'])
        json_obj['references'] = \
                json_to_PyMuTT(json_obj['references'])

        return cls(**json_obj)
