# -*- coding: utf-8 -*-
"""
pMuTT.empirical

Empirical models.
"""

import inspect
from matplotlib import pyplot as plt
import numpy as np
from pMuTT import _is_iterable
from pMuTT import constants as c
from pMuTT.io_.jsonio import json_to_pMuTT, remove_class


class EmpiricalBase:
    """The empirical parent class.
    Holds properties of a specie, the statistical-mechanical thermodynamic
    model.

    Attributes
    ----------
        name : str, optional
            Name of the specie. Default is None
        phase : str, optional
            Phase of the specie. Default is None
            G - gas.
            S - surface.
        elements : dict, optional
            Composition of the species. Default is None.
            Keys of dictionary are elements, values are stoichiometric values
            in a formula unit.
            e.g. CH3OH can be represented as:
            {'C': 1, 'H': 4, 'O': 1,}.
        statmech_model : ``pMuTT.statmech`` object, optional
            Statistical thermodynamic model. Default is None.
            Object should have the following methods: ``get_CpoR``,
            ``get_HoRT``, ``get_SoR``, ``get_GoRT``.
        references : ``pMuTT.empirical.References.references`` object, optional
            Contains references to calculate ``HoRT_ref``. If not specified
            then HoRT_dft will be used without adjustment. Default is None
        mix_models : list of ``pMuTT.mixture`` objects, optional
            Mixture models that calculate excess properties.
        smiles : str, optional
            Smiles representation of species. Default is None
        notes : str, optional
            Any additional details you would like to include such as
            computational set up. Default is None
    """

    def __init__(self, name=None, phase=None, elements=None,
                 statmech_model=None, references=None, mix_models=None,
                 smiles=None, notes=None, **kwargs):
        self.name = name
        self.phase = phase
        self.elements = elements
        self.references = references
        self.smiles = smiles
        self.notes = notes

        # Assign self.statmech_model
        if inspect.isclass(statmech_model):
            # If you're passing a class. Note that the required
            # arguments will be guessed.
            self.statmech_model = statmech_model(**kwargs)
        else:
            # If it's an object that has already been initialized
            self.statmech_model = statmech_model

        # Assign mixing models
        # TODO Mixing models can not be initialized by passing the class
        # because all the models will have the same attributes. Figure out a
        # way to pass them. Perhaps have a dictionary that contains the
        # attributes separated by species
        if not _is_iterable(mix_models) and mix_models is not None:
            mix_models = [mix_models]
        self.mix_models = mix_models

    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict

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
                Units to plot heat capacity. See :func:`~pMuTT.constants.R`
                for accepted units. If not specified, dimensionless units used.
            H_units : str
                Units to plot enthalpy. See :func:`~pMuTT.constants.R` for
                accepted units but omit the '/K' (e.g. J/mol). If not
                specified, dimensionless units used.
            S_units : str
                Units to plot entropy. See :func:`~pMuTT.constants.R` for
                accepted units. If not specified, dimensionless units used.
            G_units : str
                Units to plot Gibbs free energy. See :func:`~pMuTT.constants.R`
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
        T = np.linspace(T_low, T_high)

        f, ax = plt.subplots(4, sharex=True)
        '''
        Heat Capacity
        '''
        ax[0].set_title('Specie: {}'.format(self.name))
        Cp_plot = self.get_CpoR(T=T)
        if Cp_units is None:
            ax[0].set_ylabel('Cp/R')
        else:
            ax[0].set_ylabel('Cp ({})'.format(Cp_units))
            Cp_plot = Cp_plot*c.R(Cp_units)
        ax[0].plot(T, Cp_plot, 'r-')

        '''
        Enthalpy
        '''
        H_plot = self.get_HoRT(T=T)
        if H_units is None:
            ax[1].set_ylabel('H/RT')
        else:
            ax[1].set_ylabel('H ({})'.format(H_units))
            H_plot = H_plot*c.R('{}/K'.format(H_units))*T
        ax[1].plot(T, H_plot, 'g-')

        '''
        Entropy
        '''
        S_plot = self.get_SoR(T=T)
        if S_units is None:
            ax[2].set_ylabel('S/R')
        else:
            ax[2].set_ylabel('S ({})'.format(S_units))
            S_plot = S_plot*c.R(S_units)
        ax[2].plot(T, S_plot, 'b-')

        '''
        Gibbs energy
        '''
        ax[3].set_xlabel('Temperature (K)')
        G_plot = self.get_GoRT(T=T)
        if G_units is None:
            ax[3].set_ylabel('G/RT')
        else:
            ax[3].set_ylabel('G ({})'.format(G_units))
            G_plot = G_plot*c.R('{}/K'.format(G_units))*T
        ax[3].plot(T, G_plot, 'k-')

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
                Units to plot heat capacity. See :func:`~pMuTT.constants.R`
                for accepted units. If not specified, dimensionless units used.
            H_units : str
                Units to plot enthalpy. See :func:`~pMuTT.constants.R` for
                accepted units but omit the '/K' (e.g. J/mol). If not
                specified, dimensionless units used.
            S_units : str
                Units to plot entropy. See :func:`~pMuTT.constants.R` for
                accepted units. If not specified, dimensionless units used.
            G_units : str
                Units to plot Gibbs free energy. See :func:`~pMuTT.constants.R`
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
        T = np.linspace(T_low, T_high)

        f, ax = plt.subplots(4, sharex=True)
        '''
        Heat Capacity
        '''
        ax[0].set_title('Specie: {}'.format(self.name))
        Cp_plot = np.array([self.statmech_model.get_CpoR(T=T_i) for T_i in T])
        if Cp_units is None:
            ax[0].set_ylabel('Cp/R')
        else:
            ax[0].set_ylabel('Cp ({})'.format(Cp_units))
            Cp_plot = Cp_plot*c.R(Cp_units)
        ax[0].plot(T, Cp_plot, 'r-')

        '''
        Enthalpy
        '''
        H_plot = self.statmech_model.get_HoRT(T=T)
        if self.references is not None:
            H_plot += self.references.get_HoRT_offset(
                    descriptors=getattr(self, self.references.descriptor), T=T)

        if H_units is None:
            ax[1].set_ylabel('H/RT')
        else:
            ax[1].set_ylabel('H ({})'.format(H_units))
            H_plot = H_plot*c.R('{}/K'.format(H_units))*T
        ax[1].plot(T, H_plot, 'g-')

        '''
        Entropy
        '''
        S_plot = self.statmech_model.get_SoR(T=T)
        if S_units is None:
            ax[2].set_ylabel('S/R')
        else:
            ax[2].set_ylabel('S ({})'.format(S_units))
            S_plot = S_plot*c.R(S_units)
        ax[2].plot(T, S_plot, 'b-')

        '''
        Gibbs energy
        '''
        ax[3].set_xlabel('Temperature (K)')
        G_plot = self.statmech_model.get_GoRT(T=T)
        if self.references is not None:
            G_plot += self.references.get_HoRT_offset(
                    descriptors=getattr(self, self.references.descriptor), T=T)

        if G_units is None:
            ax[3].set_ylabel('G/RT')
        else:
            ax[3].set_ylabel('G ({})'.format(G_units))
            G_plot = G_plot*c.R('{}/K'.format(G_units))*T
        ax[3].plot(T, G_plot, 'k-')

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
                Units to plot heat capacity. See :func:`~pMuTT.constants.R` for
                accepted units. If not specified, dimensionless units used.
            H_units : str
                Units to plot enthalpy. See :func:`~pMuTT.constants.R` for
                accepted units but omit the '/K' (e.g. J/mol). If not
                specified, dimensionless units used.
            S_units : str
                Units to plot entropy. See :func:`~pMuTT.constants.R` for
                accepted units. If not specified, dimensionless units used.
            G_units : str
                Units to plot Gibbs free energy. See :func:`~pMuTT.constants.R`
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
        T = np.linspace(T_low, T_high)

        f, ax = plt.subplots(4, sharex=True)
        '''
        Heat Capacity
        '''
        ax[0].set_title('Specie: {}'.format(self.name))
        T, Cp_plot_statmech, Cp_plot_empirical = self.compare_CpoR(T=T)
        if Cp_units is None:
            ax[0].set_ylabel('Cp/R')
        else:
            ax[0].set_ylabel('Cp ({})'.format(Cp_units))
            Cp_plot_statmech = Cp_plot_statmech*c.R(Cp_units)
            Cp_plot_empirical = Cp_plot_empirical*c.R(Cp_units)

        ax[0].plot(T, Cp_plot_statmech, 'r-', label='Stat Mech Model')
        ax[0].plot(T, Cp_plot_empirical, 'b-', label='Empirical Model')
        ax[0].legend()

        '''
        Enthalpy
        '''
        T, H_plot_statmech, H_plot_empirical = self.compare_HoRT(T=T)

        if H_units is None:
            ax[1].set_ylabel('H/RT')
        else:
            ax[1].set_ylabel('H ({})'.format(H_units))
            H_plot_statmech = H_plot_statmech *\
                c.R('{}/K'.format(H_units))*T
            H_plot_empirical = H_plot_empirical *\
                c.R('{}/K'.format(H_units))*T
        ax[1].plot(T, H_plot_statmech, 'r-')
        ax[1].plot(T, H_plot_empirical, 'b-')

        '''
        Entropy
        '''
        T, S_plot_statmech, S_plot_empirical = self.compare_SoR(T=T)
        if S_units is None:
            ax[2].set_ylabel('S/R')
        else:
            ax[2].set_ylabel('S ({})'.format(S_units))
            S_plot_statmech = S_plot_statmech*c.R(S_units)
            S_plot_empirical = S_plot_empirical*c.R(S_units)
        ax[2].plot(T, S_plot_statmech, 'r-')
        ax[2].plot(T, S_plot_empirical, 'b-')

        '''
        Gibbs energy
        '''
        ax[3].set_xlabel('Temperature (K)')
        T, G_plot_statmech, G_plot_empirical = self.compare_GoRT(T=T)
        if G_units is None:
            ax[3].set_ylabel('G/RT')
        else:
            ax[3].set_ylabel('G ({})'.format(G_units))
            G_plot_statmech = G_plot_statmech *\
                c.R('{}/K'.format(G_units))*T
            G_plot_empirical = G_plot_empirical *\
                c.R('{}/K'.format(G_units))*T
        ax[3].plot(T, G_plot_statmech, 'r-')
        ax[3].plot(T, G_plot_empirical, 'b-')

        return f, ax

    def compare_CpoR(self, T=None):
        """Compares the dimensionless heat capacity of the statistical model
        and the empirical model

        Parameters
        ----------
            T : (N,) `numpy.ndarray`_ or float, optional
                Temperatures (in K) to calculate CpoR. If None, generates a
                list of temperatures between self.T_low and self.T_high
        Returns
        -------
            T : (N,) `numpy.ndarray`_ or float
                Temperatures in K
            CpoR_statmech : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of statistical thermodynamic model
            CpoR_empirical :((N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of empirical model

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """

        if T is None:
            T = np.linspace(self.T_low, self.T_high)

        try:
            iter(T)
        except TypeError:
            CpoR_statmech = self.statmech_model.get_CpoR(T=T)
            CpoR_empirical = self.get_CpoR(T=T)
        else:
            CpoR_statmech = np.zeros_like(T)
            CpoR_empirical = np.zeros_like(T)
            for i, T_i in enumerate(T):
                CpoR_statmech[i] = self.statmech_model.get_CpoR(T=T_i)
                CpoR_empirical[i] = self.get_CpoR(T=T_i)
        return (T, CpoR_statmech, CpoR_empirical)

    def compare_HoRT(self, T=None):
        """Compares the dimensionless enthalpy of the statistical model and
        the empirical model

        Parameters
        ----------
            T : (N,) `numpy.ndarray`_ or float, optional
                Temperatures (in K) to calculate CpoR. If None, generates a
                list of temperatures between self.T_low and self.T_high
        Returns
        -------
            T : (N,) `numpy.ndarray`_ or float
                Temperatures in K
            CpoR_statmech : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of statistical thermodynamic model
            CpoR_empirical :((N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of empirical model

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        if T is None:
            T = np.linspace(self.T_low, self.T_high)

        if self.references is not None:
            H_offset = self.references.get_HoRT_offset(
                    descriptors=getattr(self, self.references.descriptor), T=T)
        else:
            H_offset = np.zeros_like(T)

        try:
            iter(T)
        except TypeError:
            HoRT_statmech = np.array([self.statmech_model.get_HoRT(T=T_i)
                                      + H_offset for T_i in T])
            HoRT_empirical = self.get_HoRT(T=T)
        else:
            HoRT_statmech = np.zeros_like(T)
            HoRT_empirical = np.zeros_like(T)
            HoRT_empirical = self.get_HoRT(T=T)
            for i, T_i in enumerate(T):
                HoRT_statmech[i] = self.statmech_model.get_HoRT(T=T_i) \
                                   + H_offset[i]
        return (T, HoRT_statmech, HoRT_empirical)

    def compare_SoR(self, T=None):
        """Compares the dimensionless entropy of the statistical model and the
        empirical model

        Parameters
        ----------
            T : (N,) `numpy.ndarray`_ or float, optional
                Temperatures (in K) to calculate CpoR. If None, generates a
                list of temperatures between self.T_low and self.T_high
        Returns
        -------
            T : (N,) `numpy.ndarray`_ or float
                Temperatures in K
            CpoR_statmech : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of statistical thermodynamic model
            CpoR_empirical :((N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of empirical model

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        if T is None:
            T = np.linspace(self.T_low, self.T_high)

        try:
            iter(T)
        except TypeError:
            SoR_statmech = np.array([self.statmech_model.get_SoR(T=T_i)
                                     for T_i in T])
            SoR_empirical = self.get_SoR(T=T)
        else:
            SoR_statmech = np.zeros_like(T)
            SoR_empirical = np.zeros_like(T)
            for i, T_i in enumerate(T):
                SoR_statmech[i] = self.statmech_model.get_SoR(T=T_i)
                SoR_empirical[i] = self.get_SoR(T=T_i)
        return (T, SoR_statmech, SoR_empirical)

    def compare_GoRT(self, T=None):
        """Compares the dimensionless Gibbs energy of the statistical model and
        the empirical model

        Parameters
        ----------
            T : (N,) `numpy.ndarray`_ or float, optional
                Temperatures (in K) to calculate CpoR. If None, generates a
                list of temperatures between self.T_low and self.T_high
        Returns
        -------
            T : (N,) `numpy.ndarray`_ or float
                Temperatures in K
            CpoR_statmech : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of statistical thermodynamic model
            CpoR_empirical : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of empirical model

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        if T is None:
            T = np.linspace(self.T_low, self.T_high)

        if self.references is not None:
            G_offset = self.references.get_HoRT_offset(
                    descriptors=getattr(self, self.references.descriptor), T=T)
        else:
            G_offset = np.zeros_like(T)

        try:
            iter(T)
        except TypeError:
            GoRT_statmech = self.statmech_model.get_GoRT(T=T) + G_offset
            GoRT_empirical = self.get_GoRT(T=T)
        else:
            GoRT_statmech = np.zeros_like(T)
            GoRT_empirical = np.zeros_like(T)
            for i, T_i in enumerate(T):
                GoRT_statmech[i] = self.statmech_model.get_GoRT(T=T_i) +\
                    G_offset[i]
                GoRT_empirical[i] = self.get_GoRT(T=T_i)
        return (T, GoRT_statmech, GoRT_empirical)

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = {'class': str(self.__class__),
                    'type': 'empiricalbase',
                    'name': self.name,
                    'phase': self.phase,
                    'elements': self.elements,
                    'notes': self.notes,
                    'smiles': self.smiles, }
        try:
            obj_dict['references'] = self.references.to_dict()
        except AttributeError:
            obj_dict['references'] = self.references

        try:
            obj_dict['statmech_model'] = self.statmech_model.to_dict()
        except AttributeError:
            obj_dict['statmech_model'] = self.statmech_model

        if _is_iterable(self.mix_models):
            obj_dict['mix_models'] = \
                    [mix_model.to_dict() for mix_model in self.mix_models]
        else:
            obj_dict['mix_models'] = self.mix_models
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
        # Reconstruct statmech model
        json_obj['statmech_model'] = \
            json_to_pMuTT(json_obj['statmech_model'])
        json_obj['references'] = \
            json_to_pMuTT(json_obj['references'])
        json_obj['mix_models'] = \
            [json_to_pMuTT(mix_model) for mix_model in json_obj['mix_models']]

        return cls(**json_obj)
