# -*- coding: utf-8 -*-
"""
pmutt.empirical

Empirical models.
"""

import inspect
from matplotlib import pyplot as plt
import numpy as np
from pmutt import _is_iterable, _pmuttBase, plot_1D
from pmutt import constants as c
from pmutt.io.json import json_to_pmutt, remove_class


class EmpiricalBase(_pmuttBase):
    """The empirical parent class.
    Holds properties of a species, the statistical-mechanical thermodynamic
    model.

    Attributes
    ----------
        name : str, optional
            Name of the species. Default is None
        phase : str, optional
            Phase of the species. Default is None
            G - gas.
            S - surface.
        elements : dict, optional
            Composition of the species. Default is None.
            Keys of dictionary are elements, values are stoichiometric values
            in a formula unit.
            e.g. CH3OH can be represented as:
            {'C': 1, 'H': 4, 'O': 1,}.
        statmech_model : ``pmutt.statmech`` object, optional
            Statistical thermodynamic model. Default is None.
            Object should have the following methods: ``get_CpoR``,
            ``get_HoRT``, ``get_SoR``, ``get_GoRT``.
        misc_models : list of pmutt model objects, optional
            Extra models to add extra functionality. Commonly used models would
            be ``pmutt.mixture`` models
        smiles : str, optional
            Smiles representation of species. Default is None
        notes : str, optional
            Any additional details you would like to include such as
            computational set up. Default is None
    """

    def __init__(self, name=None, phase=None, elements=None,
                 statmech_model=None, references=None, misc_models=None,
                 smiles=None, notes=None, **kwargs):
        self.name = name
        self.phase = phase
        self.elements = elements
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
        if not _is_iterable(misc_models) and misc_models is not None:
            misc_models = [misc_models]
        self.misc_models = misc_models

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
                Units to plot heat capacity. See :func:`~pmutt.constants.R`
                for accepted units. If not specified, dimensionless units used.
            H_units : str
                Units to plot enthalpy. See :func:`~pmutt.constants.R` for
                accepted units but omit the '/K' (e.g. J/mol). If not
                specified, dimensionless units used.
            S_units : str
                Units to plot entropy. See :func:`~pmutt.constants.R` for
                accepted units. If not specified, dimensionless units used.
            G_units : str
                Units to plot Gibbs free energy. See :func:`~pmutt.constants.R`
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
        # Process temperatures
        if T_low is None:
            T_low = self.T_low
        if T_high is None:
            T_high = self.T_high
        T = np.linspace(T_low, T_high)

        # Process the method names, units, and y labels
        units = (Cp_units, H_units, S_units, G_units)
        methods = ['get_Cp', 'get_H', 'get_S', 'get_G']
        y_labels = []
        kwargs = {}
        for i, (units, method) in enumerate(zip(units, methods)):
            if units is None:
                if method == 'get_Cp' or method == 'get_S':
                    methods[i] = '{}oR'.format(method)
                else:
                    methods[i] = '{}oRT'.format(method)
                y_labels.append(method.replace('o', '/').replace('get_', ''))
            else:
                kwargs['{}_kwargs'.format(method)] = {'units': units}
                y_labels.append('{} ({})'.format(method.replace('get_', ''),
                                                 units))

        fig, ax = plot_1D(self, x_name='T', x_values=T, methods=methods,
                          **kwargs)
        
        # Add titles and labels
        ax[0].set_title('Species: {}'.format(self.name))
        for i, y_label in enumerate(y_labels):
            ax[i].set_xlabel('Temperature (K)')            
            ax[i].set_ylabel(y_label)
        return fig, ax

    def plot_statmech(self, T_low=None, T_high=None, Cp_units=None,
                      H_units=None, S_units=None, G_units=None,
                      use_references=True):
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
                Units to plot heat capacity. See :func:`~pmutt.constants.R`
                for accepted units. If not specified, dimensionless units used.
            H_units : str
                Units to plot enthalpy. See :func:`~pmutt.constants.R` for
                accepted units but omit the '/K' (e.g. J/mol). If not
                specified, dimensionless units used.
            S_units : str
                Units to plot entropy. See :func:`~pmutt.constants.R` for
                accepted units. If not specified, dimensionless units used.
            G_units : str
                Units to plot Gibbs free energy. See :func:`~pmutt.constants.R`
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
        # Process temperatures
        if T_low is None:
            T_low = self.T_low
        if T_high is None:
            T_high = self.T_high
        T = np.linspace(T_low, T_high)

        # Process the method names, units, and y labels
        units = (Cp_units, H_units, S_units, G_units)
        methods = ['get_Cp', 'get_H', 'get_S', 'get_G']
        y_labels = []
        kwargs = {}
        for i, (units, method) in enumerate(zip(units, methods)):
            if units is None:
                if method == 'get_Cp' or method == 'get_S':
                    methods[i] = '{}oR'.format(method)
                else:
                    methods[i] = '{}oRT'.format(method)
                y_labels.append(method.replace('o', '/').replace('get_', ''))
            else:
                kwargs['{}_kwargs'.format(method)] = {'units': units}
                y_labels.append('{} ({})'.format(method.replace('get_', ''),
                                                 units))

        fig, ax = plot_1D(self.statmech_model, x_name='T', x_values=T,
                          methods=methods, use_references=use_references,
                          **kwargs)
        
        # Add titles and labels
        ax[0].set_title('Species: {}'.format(self.name))
        for i, y_label in enumerate(y_labels):
            ax[i].set_xlabel('Temperature (K)')            
            ax[i].set_ylabel(y_label)
        return fig, ax

    def plot_statmech_and_empirical(self, T_low=None, T_high=None,
                                    Cp_units=None, H_units=None,
                                    S_units=None, G_units=None,
                                    use_references=True):
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
                Units to plot heat capacity. See :func:`~pmutt.constants.R` for
                accepted units. If not specified, dimensionless units used.
            H_units : str
                Units to plot enthalpy. See :func:`~pmutt.constants.R` for
                accepted units but omit the '/K' (e.g. J/mol). If not
                specified, dimensionless units used.
            S_units : str
                Units to plot entropy. See :func:`~pmutt.constants.R` for
                accepted units. If not specified, dimensionless units used.
            G_units : str
                Units to plot Gibbs free energy. See :func:`~pmutt.constants.R`
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
        # Process temperatures
        if T_low is None:
            T_low = self.T_low
        if T_high is None:
            T_high = self.T_high
        T = np.linspace(T_low, T_high)

        # Process the method names, units, and y labels
        units = (Cp_units, H_units, S_units, G_units)
        methods = ['get_Cp', 'get_H', 'get_S', 'get_G']
        y_labels = []
        kwargs = {}
        for i, (units, method) in enumerate(zip(units, methods)):
            if units is None:
                if method == 'get_Cp' or method == 'get_S':
                    methods[i] = '{}oR'.format(method)
                else:
                    methods[i] = '{}oRT'.format(method)
                y_labels.append(method.replace('o', '/').replace('get_', ''))
            else:
                kwargs['{}_kwargs'.format(method)] = {'units': units}
                y_labels.append('{} ({})'.format(method.replace('get_', ''),
                                                 units))

        fig, ax = plot_1D(self, x_name='T', x_values=T, methods=methods,
                          **kwargs)
        fig, ax = plot_1D(self.statmech_model, x_name='T', x_values=T,
                          methods=methods, figure=fig, ax=ax,
                          use_references=use_references, **kwargs)
        
        # Add titles and labels
        ax[0].set_title('Species: {}'.format(self.name))
        ax[0].legend(['Empirical', 'StatMech'])
        for i, y_label in enumerate(y_labels):
            ax[i].set_xlabel('Temperature (K)')            
            ax[i].set_ylabel(y_label)
        return fig, ax

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

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
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

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
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

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
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

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
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
            obj_dict['statmech_model'] = self.statmech_model.to_dict()
        except AttributeError:
            obj_dict['statmech_model'] = self.statmech_model

        if _is_iterable(self.misc_models):
            obj_dict['misc_models'] = \
                    [mix_model.to_dict() for mix_model in self.misc_models]
        else:
            obj_dict['misc_models'] = self.misc_models
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
        json_obj['statmech_model'] = json_to_pmutt(json_obj['statmech_model'])
        json_obj['misc_models'] = \
            [json_to_pmutt(mix_model) for mix_model in json_obj['misc_models']]

        return cls(**json_obj)
