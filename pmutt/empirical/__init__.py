# -*- coding: utf-8 -*-
"""
pmutt.empirical

Empirical models.
"""

import inspect

import numpy as np
from matplotlib import pyplot as plt

from pmutt import plot_1D
from pmutt import constants as c
from pmutt.io.json import json_to_pmutt, remove_class
from pmutt import _is_iterable, _ModelBase, _pmuttBase


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
        model : model object, optional
            Object which was used to fit the Nasa polynomial. Default is None.
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
        add_gas_P_adj : bool, optional
            If True, gas species (i.e. `phase`='gas' or `phase`='g') will
            automatically be assigned :class:`~pmutt.empirical.GasPressureAdj`.
            Default is True.
    """

    def __init__(self, name=None, phase=None, elements=None,
                 model=None, misc_models=None, smiles=None, notes=None,
                 add_gas_P_adj=True, **kwargs):
        self.name = name
        self.phase = phase
        self.elements = elements
        self.smiles = smiles
        self.notes = notes

        # Assign self.model
        if inspect.isclass(model):
            # If you're passing a class. Note that the required
            # arguments will be guessed.
            self.model = model(**kwargs)
        else:
            # If it's an object that has already been initialized
            self.model = model

        # Assign mixing models
        # TODO Mixing models can not be initialized by passing the class
        # because all the models will have the same attributes. Figure out a
        # way to pass them. Perhaps have a dictionary that contains the
        # attributes separated by species
        
        # Single misc model is passed
        if not _is_iterable(misc_models) and misc_models is not None:
            misc_models = [misc_models]

        # Assign pressure adjustment if necessary
        dict_entry = {'class': "<class 'pmutt.empirical.GasPressureAdj'>"}
        if self.phase is not None:
            # If the species is a gas
            if (self.phase.lower() == 'g' or self.phase.lower() == 'gas'):
                # If the species has not already been assigned GasPressureAdj
                if misc_models is None:
                    misc_models = [GasPressureAdj()]
                else:
                    for i, model in enumerate(misc_models):
                        if model == dict_entry:
                            misc_models[i] = GasPressureAdj()
                            break
                        elif isinstance(model, GasPressureAdj):
                            break
                    else:
                        misc_models.append(GasPressureAdj())
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

        fig, ax = plot_1D(self.model, x_name='T', x_values=T,
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
        fig, ax = plot_1D(self.model, x_name='T', x_values=T,
                          methods=methods, figure=fig, ax=ax,
                          use_references=use_references, **kwargs)

        # Get class name of current model object
        model_class1 = str(self.__class__)
        model_class1 = model_class1.split('.')[-1].replace("'>", "")

        # Get class name of derived model object
        model_class2 = str(self.model.__class__)
        model_class2 = model_class2.split('.')[-1].replace("'>", "")

        # Add titles and labels
        ax[0].set_title('Species: {}'.format(self.name))
        ax[0].legend([model_class1, model_class2])
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
            CpoR_model : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of original model
            CpoR_empirical :((N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of empirical model

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """

        if T is None:
            T = np.linspace(self.T_low, self.T_high)

        CpoR_empirical = self.get_CpoR(T=T)
        try:
            CpoR_model = self.model.get_CpoR(T=T)
        except ValueError:
            CpoR_model = np.array([self.model.get_CpoR(T=T_i) for T_i in T])
        return (T, CpoR_model, CpoR_empirical)

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
            CpoR_model : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of original model
            CpoR_empirical :((N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of empirical model

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        if T is None:
            T = np.linspace(self.T_low, self.T_high)

        HoRT_empirical = self.get_HoRT(T=T)
        try:
            HoRT_model = self.model.get_HoRT(T=T)
        except ValueError:
            HoRT_model = np.array([self.model.get_HoRT(T=T_i) for T_i in T])
        return (T, HoRT_model, HoRT_empirical)

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
            CpoR_model : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of original model
            CpoR_empirical :((N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of empirical model

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        if T is None:
            T = np.linspace(self.T_low, self.T_high)

        SoR_empirical = self.get_SoR(T=T)
        try:
            SoR_model = self.model.get_SoR(T=T)
        except ValueError:
            SoR_model = np.array([self.model.get_SoR(T=T_i) for T_i in T])
        return (T, SoR_model, SoR_empirical)

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
            CpoR_model : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of original model
            CpoR_empirical : (N,) `numpy.ndarray`_ or float
                Dimensionless heat capacity of empirical model

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        if T is None:
            T = np.linspace(self.T_low, self.T_high)

        GoRT_empirical = self.get_GoRT(T=T)
        try:
            GoRT_model = self.model.get_GoRT(T=T)
        except ValueError:
            GoRT_model = np.array([self.model.get_GoRT(T=T_i) for T_i in T])
        return (T, GoRT_model, GoRT_empirical)

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
            obj_dict['model'] = self.model.to_dict()
        except AttributeError:
            obj_dict['model'] = self.model

        if _is_iterable(self.misc_models):
            obj_dict['misc_models'] = \
                    [misc_model.to_dict() for misc_model in self.misc_models]
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
        # Reconstruct model
        json_obj['model'] = json_to_pmutt(json_obj['model'])
        if json_obj['misc_models'] is not None:
            json_obj['misc_models'] = [json_to_pmutt(misc_model) \
                for misc_model in json_obj['misc_models']]

        return cls(**json_obj)

class GasPressureAdj(_ModelBase):
    """Includes pressure's effect on entropy for gas molecules"""

    def __init__(self):
        pass

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = {'class': str(self.__class__)}
        return obj_dict

    def get_CvoR(self):
        return 0.

    def get_CpoR(self):
        return 0.

    def get_UoRT(self):
        return 0.

    def get_HoRT(self):
        return 0.

    def get_SoR(self, P=c.P0('bar')):
        """Calculates dimesionless entropy

        :math:`\\frac{S}{R} = -\\ln\\bigg(\\frac{P}{P_0}\\bigg)`

        Parameters
        ----------
            P : float or `numpy.ndarray`_, optional
                Pressure in bar. Default is P0 (1 bar)
        Returns
        -------
            SoR : float or `numpy.ndarray`_
                Dimensionless adjustment to entropy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        return -np.log(P)