# -*- coding: utf-8 -*-
from collections import Counter
import re
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from pMuTT import _force_pass_arguments, _is_iterable
from pMuTT import constants as c
from pMuTT.io_.jsonio import json_to_pMuTT, remove_class


class Reaction:
    """Represents a chemical reaction

    Attributes
    ----------
        reactants : list of pMuTT model objects
            Reactants
        reactants_stoich : list of float
            Stoichiometric quantities of reactants
        products : list of pMuTT model objects
            Products
        products_stoich : list of float
            Stoichiometric quantities of products
        transition_state : pMuTT model object, optional
            Transition state specie. Default is None
        transition_state_stoich : list of float, optional
            Stoichiometric quantities of transition state species.
            Default is None
    """

    def __init__(self, reactants, reactants_stoich, products, products_stoich,
                 transition_state=None, transition_state_stoich=None):
        # If any of the entries were not iterable, assign them to a list
        if not _is_iterable(reactants):
            reactants = [reactants]
        if not _is_iterable(reactants_stoich):
            reactants_stoich = [reactants_stoich]
        if not _is_iterable(products):
            products = [products]
        if not _is_iterable(products_stoich):
            products_stoich = [products_stoich]
        if not _is_iterable(transition_state):
            transition_state = [transition_state]
        if not _is_iterable(transition_state_stoich):
            transition_state_stoich = [transition_state_stoich]
        self.reactants = reactants
        self.reactants_stoich = reactants_stoich
        self.products = products
        self.products_stoich = products_stoich
        self.transition_state = transition_state
        self.transition_state_stoich = transition_state_stoich

    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict

    def __str__(self):
        return self.to_str()

    def check_element_balance(self):
        """Checks the reactants, products and transition state elemental
        composition

        Raises
        ------
            ValueError
                Raised if the reactants, products and/or transition state
                element composition does not agree.
        """
        reactant_elements = _count_elements(self.reactants,
                                            self.reactants_stoich)
        product_elements = _count_elements(self.products,
                                           self.products_stoich)
        if reactant_elements != product_elements:
            raise ValueError('Number of elements in reactants and products do '
                             'not agree.\nReactant count: {}\n'
                             'Product count: {}'.format(reactant_elements,
                                                        product_elements))

        if self.transition_state != [None]:
            TS_elements = _count_elements(self.transition_state,
                                          self.transition_state_stoich)
            if reactant_elements != TS_elements:
                raise ValueError('Number of elements in reactants and '
                                 'transition state do not agree.\n'
                                 'Reactant count: {}\n'
                                 'Product count: {}'.format(reactant_elements,
                                                            TS_elements))

    def get_delta_q(self, rev=False, **kwargs):
        """Gets change in partition function between reactants and products

        Parameters
        ----------
            kwargs : keyword arguments
                Parameters required to calculate partition function
        Returns
        -------
            delta_q : float
                Change in partition function between reactants and products
        """
        delta_q = _get_q_rxn(initial_state=self.reactants,
                             initial_state_stoich=self.reactants_stoich,
                             final_state=self.products,
                             final_state_stoich=self.products_stoich,
                             **kwargs)
        if rev:
            return 1./delta_q
        else:
            return delta_q

    def get_delta_CvoR(self, rev=False, **kwargs):
        """Gets change in dimensionless heat capacity between reactants and
        products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate heat capacity
        Returns
        -------
            delta_CvoR : float
                Change in heat capacity between reactants and products
        """
        delta_CvoR = _get_CvoR_rxn(initial_state=self.reactants,
                                   initial_state_stoich=self.reactants_stoich,
                                   final_state=self.products,
                                   final_state_stoich=self.products_stoich,
                                   **kwargs)
        if rev:
            return -delta_CvoR
        else:
            return delta_CvoR

    def get_delta_CpoR(self, rev=False, **kwargs):
        """Gets change in dimensionless heat capacity between reactants and
        products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate heat capacity
        Returns
        -------
            delta_CvoR : float
                Change in heat capacity between reactants and products
        """
        delta_CpoR = _get_CpoR_rxn(initial_state=self.reactants,
                                   initial_state_stoich=self.reactants_stoich,
                                   final_state=self.products,
                                   final_state_stoich=self.products_stoich,
                                   **kwargs)
        if rev:
            return -delta_CpoR
        else:
            return delta_CpoR

    def get_delta_UoRT(self, rev=False, **kwargs):
        """Gets change in dimensionless internal energy between reactants and
        products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate internal energy
        Returns
        -------
            delta_UoRT : float
                Change in internal energy between reactants and products
        """
        delta_UoRT = _get_UoRT_rxn(initial_state=self.reactants,
                                   initial_state_stoich=self.reactants_stoich,
                                   final_state=self.products,
                                   final_state_stoich=self.products_stoich,
                                   **kwargs)
        if rev:
            return -delta_UoRT
        else:
            return delta_UoRT

    def get_delta_HoRT(self, rev=False, **kwargs):
        """Gets change in dimensionless enthalpy between reactants and products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate enthalpy
        Returns
        -------
            delta_HoRT : float
                Change in enthalpy between reactants and products
        """
        delta_HoRT = _get_HoRT_rxn(initial_state=self.reactants,
                                   initial_state_stoich=self.reactants_stoich,
                                   final_state=self.products,
                                   final_state_stoich=self.products_stoich,
                                   **kwargs)
        if rev:
            return -delta_HoRT
        else:
            return delta_HoRT

    def get_delta_SoR(self, rev=False, **kwargs):
        """Gets change in dimensionless entropy between reactants and products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate entropy
        Returns
        -------
            delta_SoR : float
                Change in entropy between reactants and products
        """
        delta_SoR = _get_SoR_rxn(initial_state=self.reactants,
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=self.products,
                                 final_state_stoich=self.products_stoich,
                                 **kwargs)
        if rev:
            return -delta_SoR
        else:
            return delta_SoR

    def get_delta_AoRT(self, rev=False, **kwargs):
        """Gets change in dimensionless Helmholtz energy between reactants and
        products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate Helmholtz energy
        Returns
        -------
            delta_AoRT : float
                Change in Helmholtz energy between reactants and products
        """
        delta_AoRT = _get_AoRT_rxn(initial_state=self.reactants,
                                   initial_state_stoich=self.reactants_stoich,
                                   final_state=self.products,
                                   final_state_stoich=self.products_stoich,
                                   **kwargs)
        if rev:
            return -delta_AoRT
        else:
            return delta_AoRT

    def get_delta_GoRT(self, rev=False, **kwargs):
        """Gets change in dimensionless Gibbs energy between reactants and
        products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate Gibbs energy
        Returns
        -------
            delta_GoRT : float
                Change in Gibbs energy between reactants and products
        """
        delta_GoRT = _get_GoRT_rxn(initial_state=self.reactants,
                                   initial_state_stoich=self.reactants_stoich,
                                   final_state=self.products,
                                   final_state_stoich=self.products_stoich,
                                   **kwargs)
        if rev:
            return -delta_GoRT
        else:
            return delta_GoRT

    def get_Keq(self, rev=False, **kwargs):
        """Gets equilibrium constant between reactants and products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate equilibrium constant
        Returns
        -------
            delta_GoRT : float
                Change in equilibrium constant between reactants and products
        """
        return np.exp(-self.get_delta_GoRT(rev=rev, **kwargs))

    def get_q_act(self, rev=False, **kwargs):
        """Gets change in partition function between reactants (or products)
        and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate partition function
        Returns
        -------
            q_act : float
                Change in partition function between reactants and transition
                state
        """
        if rev:
            return _get_q_rxn(initial_state=self.products,
                              initial_state_stoich=self.products_stoich,
                              final_state=self.transition_state,
                              final_state_stoich=self.transition_state_stoich,
                              **kwargs)
        else:
            return _get_q_rxn(initial_state=self.reactants,
                              initial_state_stoich=self.reactants_stoich,
                              final_state=self.transition_state,
                              final_state_stoich=self.transition_state_stoich,
                              **kwargs)

    def get_CvoR_act(self, rev=False, **kwargs):
        """Gets change in dimensionless heat capacity between reactants (or
        products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless heat capacity
        Returns
        -------
            CvoR_act : float
                Change in dimensionless heat capacity between reactants and
                transition state
        """
        if rev:
            return _get_CvoR_rxn(initial_state=self.products,
                                 initial_state_stoich=self.products_stoich,
                                 final_state=self.transition_state,
                                 final_state_stoich=self.transition_state_stoich,
                                 **kwargs)
        else:
            return _get_CvoR_rxn(initial_state=self.reactants,
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=self.transition_state,
                                 final_state_stoich=self.transition_state_stoich,
                                 **kwargs)

    def get_CpoR_act(self, rev=False, **kwargs):
        """Gets change in dimensionless heat capacity between reactants (or
        products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless heat capacity
        Returns
        -------
            CpoR_act : float
                Change in dimensionless heat capacity between reactants and
                transition state
        """
        if rev:
            return _get_CpoR_rxn(initial_state=self.products,
                                 initial_state_stoich=self.products_stoich,
                                 final_state=self.transition_state,
                                 final_state_stoich=self.transition_state_stoich,
                                 **kwargs)
        else:
            return _get_CpoR_rxn(initial_state=self.reactants,
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=self.transition_state,
                                 final_state_stoich=self.transition_state_stoich,
                                 **kwargs)

    def get_UoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless internal energy between reactants (or
        products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless internal energy
        Returns
        -------
            UoRT_act : float
                Change in dimensionless internal energy between reactants and
                transition state
        """
        if rev:
            return _get_UoRT_rxn(initial_state=self.products,
                                 initial_state_stoich=self.products_stoich,
                                 final_state=self.transition_state,
                                 final_state_stoich=self.transition_state_stoich,
                                 **kwargs)
        else:
            return _get_UoRT_rxn(initial_state=self.reactants,
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=self.transition_state,
                                 final_state_stoich=self.transition_state_stoich,
                                 **kwargs)

    def get_HoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless enthalpy between reactants
        (or products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless enthalpy
        Returns
        -------
            HoRT_act : float
                Change in dimensionless enthalpy between reactants and
                transition state
        """
        if rev:
            return _get_HoRT_rxn(initial_state=self.products,
                                 initial_state_stoich=self.products_stoich,
                                 final_state=self.transition_state,
                                 final_state_stoich=self.transition_state_stoich,
                                 **kwargs)
        else:
            return _get_HoRT_rxn(initial_state=self.reactants,
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=self.transition_state,
                                 final_state_stoich=self.transition_state_stoich,
                                 **kwargs)

    def get_SoR_act(self, rev=False, **kwargs):
        """Gets change in dimensionless entropy between reactants (or products)
        and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless entropy
        Returns
        -------
            SoR_act : float
                Change in dimensionless entropy between reactants and
                transition state
        """
        if rev:
            return _get_SoR_rxn(initial_state=self.products,
                                initial_state_stoich=self.products_stoich,
                                final_state=self.transition_state,
                                final_state_stoich=self.transition_state_stoich,
                                **kwargs)
        else:
            return _get_SoR_rxn(initial_state=self.reactants,
                                initial_state_stoich=self.reactants_stoich,
                                final_state=self.transition_state,
                                final_state_stoich=self.transition_state_stoich,
                                **kwargs)

    def get_AoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless Helmholtz energy between reactants (or
        products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless Helmholtz energy
        Returns
        -------
            AoRT_act : float
                Change in dimensionless Helmholtz energy between reactants and
                transition state
        """
        if rev:
            return _get_AoRT_rxn(initial_state=self.products,
                                 initial_state_stoich=self.products_stoich,
                                 final_state=self.transition_state,
                                 final_state_stoich=self.transition_state_stoich,
                                 **kwargs)
        else:
            return _get_AoRT_rxn(initial_state=self.reactants,
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=self.transition_state,
                                 final_state_stoich=self.transition_state_stoich,
                                 **kwargs)

    def get_GoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless Gibbs energy between reactants (or
        products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless Gibbs energy
        Returns
        -------
            GoRT_act : float
                Change in dimensionless Gibbs energy between reactants and
                transition state
        """
        if rev:
            return _get_GoRT_rxn(initial_state=self.products,
                                 initial_state_stoich=self.products_stoich,
                                 final_state=self.transition_state,
                                 final_state_stoich=self.transition_state_stoich,
                                 **kwargs)
        else:
            return _get_GoRT_rxn(initial_state=self.reactants,
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=self.transition_state,
                                 final_state_stoich=self.transition_state_stoich,
                                 **kwargs)

    def get_A(self, T=c.T0('K'), rev=False, **kwargs):
        """Gets pre-exponential factor between reactants (or products) and
        transition state in 1/s

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            T : float, optional
                Temperature in K. Default is standard temperature.
            kwargs : keyword arguments
                Parameters required to calculate pre-exponential factor
        Returns
        -------
            A : float
                Pre-exponential factor
        """
        return c.kb('J/K')*T/c.h('J s')\
            * np.exp(self.get_SoR_act(rev=rev, T=c.T0('K'), **kwargs))

    @classmethod
    def from_string(cls, reaction_str, species, species_delimiter='+',
                    reaction_delimiter='='):
        """Create a reaction object using the reaction string

        Parameters
        ----------
            reaction_str : str
                Reaction string.
            species : dict
                Dictionary using the names as keys. If you have a list of
                species, use pMuTT.models.pMuTT_list_to_dict to make a dict.
            species_delimiter : str, optional
                Delimiter that separate species. Leading and trailing spaces
                will be trimmed. Default is '+'
            reaction_delimiter : str, optional
                Delimiter that separate sides of the reaction. Leading and
                trailing spaces will be trimmed. Default is '='
        Returns
        -------
            Reaction : Reaction object
        """
        (react_names, react_stoich, prod_names, prod_stoich,
         ts_names, ts_stoich) = _parse_reaction(
                reaction_str=reaction_str,
                species_delimiter=species_delimiter,
                reaction_delimiter=reaction_delimiter)
        reactants = [species[name] for name in react_names]
        products = [species[name] for name in prod_names]
        if ts_names is None:
            ts = None
        else:
            ts = [species[name] for name in ts_names]
        return cls(reactants=reactants, reactants_stoich=react_stoich,
                   products=products, products_stoich=prod_stoich,
                   transition_state=ts, transition_state_stoich=ts_stoich)

    def to_str(self, species_delimiter='+', reaction_delimiter='='):
        """Writes the Reaction object as a stoichiometric reaction

        Parameters
        ----------
            species_delimiter : str, optional
                Separates species. Default is '+'
            reaction_delimiter : str, optional
                Separates reaction sides. Default is '='
        Returns
        -------
            reaction_str : str
                Reaction string
        """
        # Write reactants
        reaction_str = _write_reaction_side(species=self.reactants,
                                            stoich=self.reactants_stoich,
                                            species_delimiter=species_delimiter)
        reaction_str += reaction_delimiter

        # Write transition state if any
        if self.transition_state is not None:
            reaction_str += _write_reaction_side(
                    species=self.transition_state,
                    stoich=self.transition_state_stoich,
                    species_delimiter=species_delimiter)
            reaction_str += reaction_delimiter

        # Write products
        reaction_str += _write_reaction_side(species=self.products,
                                             stoich=self.products_stoich,
                                             species_delimiter=species_delimiter)
        return reaction_str

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = {
            'class': str(self.__class__),
            'reactants': [reactant.to_dict() for reactant in self.reactants],
            'reactants_stoich': list(self.reactants_stoich),
            'products': [product.to_dict() for product in self.products],
            'products_stoich': list(self.products_stoich)}
        try:
            obj_dict['transition_state'] = self.transition_state.to_dict()
        except AttributeError:
            obj_dict['transition_state'] = self.transition_state

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
            Reaction : Reaction object
        """
        json_obj = remove_class(json_obj)
        json_obj['reactants'] = [json_to_pMuTT(reactant)
                                 for reactant in json_obj['reactants']]
        json_obj['products'] = [json_to_pMuTT(product)
                                for product in json_obj['products']]
        json_obj['transition_state'] = json_to_pMuTT(
                json_obj['transition_state'])
        return cls(**json_obj)


class PhaseDiagram:
    """Generate phase diagrams based on reactions specified.

    Attributes
    ----------
        reactions : list of ``pMuTT.models.reaction.Reaction`` objects
            Formation reactions for each phase. Reactions should be written
            with consistent reference species to obtain meaningful data.
        norm_factors : (N,) `numpy.ndarray`_ of float, optional
            Used for normalizing Gibbs energies. These factors could be
            surface areas when calculating surface energies or if the
            reactions stoichiometry is not consistent. Default is an array of
            1. It should have the same length as reactions.

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
    """

    def __init__(self, reactions, norm_factors=None):
        self.reactions = reactions
        if norm_factors is None:
            self.norm_factors = np.ones(len(reactions))
        else:
            self.norm_factors = norm_factors

    def get_GoRT_1D(self, x_name, x_values, G_units=None, **kwargs):
        """Calculates the Gibbs free energy for all the reactions for 1 varying
        parameter

        Parameters
        ----------
            x_name : str
                Name of variable to vary
            x_values : iterable object
                x values to use
            G_units : str, optional
                Units for G. If None, uses GoRT. Default is None
            kwargs : keyword arguments
                Other variables to use in the calculation
        Returns
        -------
            GoRT : (M, N) `numpy.ndarray`_ of float
                GoRT values. The first index corresponds to the number of
                reactions. The second index corresponds to the conditions
                specified by x_values.
            stable_phases : (N,) `numpy.ndarray`_ of int
                Each element of the array corresponds to the index of the most
                stable phase at the x_values.
        """
        GoRT = np.zeros(shape=(len(self.reactions), len(x_values)))
        for i, (reaction, norm_factor) in enumerate(zip(self.reactions,
                                                        self.norm_factors)):
            for j, x in enumerate(x_values):
                kwargs[x_name] = x
                GoRT[i, j] = reaction.get_delta_GoRT(**kwargs)/norm_factor

                # Add unit corrections
                if G_units is not None:
                    GoRT[i, j] *= c.R('{}/K'.format(G_units))*kwargs['T']
        stable_phases = np.nanargmin(GoRT, axis=1)
        return (GoRT, stable_phases)

    def plot_1D(self, x_name, x_values, G_units=None, **kwargs):
        """Make a 1D phase diagram.

        Parameters
        ----------
            x_name : str
                Name of variable to vary
            x_values : iterable object
                x values to use
            G_units : str, optional
                Units for G. If None, uses GoRT. Default is None
            kwargs : keyword arguments
                Other variables to use in the calculation
        Returns
        -------
            figure : `matplotlib.figure.Figure`_
                Figure
            ax : `matplotlib.axes.Axes.axis`_
                Axes of the plots.

        .. _`matplotlib.figure.Figure`: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html
        .. _`matplotlib.axes.Axes.axis`: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axis.html
        """
        fig, ax = plt.subplots()
        GoRT, stable_phases = self.get_GoRT_1D(x_name=x_name,
                                               x_values=x_values,
                                               G_units=G_units, **kwargs)
        for GoRT_rxn, rxn in zip(GoRT, self.reactions):
            plt.plot(x_values, GoRT_rxn, label=rxn.to_str())
        ax.legend()
        ax.set_xlabel(x_name)
        if G_units is None:
            ax.set_ylabel('G/RT')
        else:
            ax.set_ylabel('G ({})'.format(G_units))
        return (fig, ax)

    def get_GoRT_2D(self, x1_name, x1_values, x2_name, x2_values,
                    G_units=None, **kwargs):
        """Calculates the Gibbs free energy for all the reactions for two
        varying parameters

        Parameters
        ----------
            x1_name : str
                Name of first variable to vary
            x1_values : iterable object
                x1 values to use
            x2_name : str
                Name of first variable to vary
            x2_values : iterable object
                x2 values to use
            G_units : str, optional
                Units for G. If None, uses GoRT. Default is None
            kwargs : keyword arguments
                Other variables to use in the calculation
        Returns
        -------
            GoRT : (M, N, O) `numpy.ndarray`_ of float
                GoRT values. The first index corresponds to the number of
                reactions. The second index corresponds to the conditions
                specified by x_values.
            stable_phases : (N, O) `numpy.ndarray`_ of int
                Each element of the array corresponds to the index of the most
                stable phase at the x_values.
        """
        GoRT = np.zeros(
                shape=(len(self.reactions), len(x1_values), len(x2_values)))
        for i, (reaction, norm_factor) in enumerate(zip(self.reactions,
                                                        self.norm_factors)):
            for j, x1 in enumerate(x1_values):
                kwargs[x1_name] = x1
                for k, x2 in enumerate(x2_values):
                    kwargs[x2_name] = x2
                    GoRT[i, j, k] = \
                        reaction.get_delta_GoRT(**kwargs)/norm_factor
                    # Add unit corrections
                    if G_units is not None:
                        GoRT[i, j, k] *=\
                            c.R('{}/K'.format(G_units))*kwargs['T']
        # Take a transpose
        GoRT_T = GoRT.transpose((1, 2, 0))
        stable_phases = np.zeros((len(x1_values), len(x2_values)))
        for i, GoRT_row in enumerate(GoRT_T):
            stable_phases[i, :] = np.nanargmin(GoRT_row, axis=1)

        return GoRT, stable_phases

    def plot_2D(self, x1_name, x1_values, x2_name, x2_values, G_units=None,
                **kwargs):
        """Make a 2D phase diagram.

        Parameters
        ----------
            x1_name : str
                Name of first variable to vary
            x1_values : iterable object
                x1 values to use
            x2_name : str
                Name of first variable to vary
            x2_values : iterable object
                x2 values to use
            G_units : str, optional
                Units for G. If None, uses GoRT. Default is None
            kwargs : keyword arguments
                Other variables to use in the calculation
        Returns
        -------
            figure : `matplotlib.figure.Figure`_
                Figure
            ax : `matplotlib.axes.Axes.axis`_
                Axes of the plots.
            c : `matplotlib.collections.QuadMesh`_
                Heatmap plot
            cbar : `matplotlib.colorbar.Colorbar`_
                Colorbar for plot

        .. _`matplotlib.figure.Figure`: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html
        .. _`matplotlib.axes.Axes.axis`: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axis.html
        .. _`matplotlib.collections.QuadMesh`: https://matplotlib.org/api/collections_api.html#matplotlib.collections.QuadMesh
        .. _`matplotlib.colorbar.Colorbar`: https://matplotlib.org/api/_as_gen/matplotlib.pyplot.colorbar.html
        """
        # Process input data
        x2_mesh, x1_mesh = np.meshgrid(x2_values, x1_values)
        GoRT, stable_phases = self.get_GoRT_2D(x1_name=x1_name,
                                               x1_values=x1_values,
                                               x2_name=x2_name,
                                               x2_values=x2_values,
                                               G_units=G_units, **kwargs)

        fig, ax = plt.subplots()
        # Choosing color palette
        cmap = plt.get_cmap('viridis')
        norm = matplotlib.colors.BoundaryNorm(np.arange(len(self.reactions)+1),
                                              cmap.N)
        # Create colormap
        c = plt.pcolormesh(x1_mesh, x2_mesh, stable_phases, cmap=cmap,
                           norm=norm, vmin=0, vmax=len(self.reactions))
        # Set colorbar
        cbar = fig.colorbar(c, ticks=np.arange(len(self.reactions))+0.5)
        cbar.ax.set_yticklabels(
                [reaction.to_str() for reaction in self.reactions])
        # Set axis labels
        ax.set_xlabel(x1_name)
        ax.set_ylabel(x2_name)
        return (fig, ax, c, cbar)

def _get_q_rxn(initial_state, initial_state_stoich, final_state,
               final_state_stoich, **kwargs):
    """Helper function to calculate partition function

    Parameters
    ----------
        initial_state : list of pMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of pMuTT model objects
            Final state
        final_state_stoich : list of float
            Stoichiometry of final state
        kwargs : keyword arguments
            Parameters to calculate q
    Returns
    -------
        q : float
            Partition function between initial state and final state
    """
    q = 1.
    for specie, stoich in zip(final_state, final_state_stoich):
        q *= _force_pass_arguments(specie.get_q, **kwargs)**stoich
    for specie, stoich in zip(initial_state, initial_state_stoich):
        q /= _force_pass_arguments(specie.get_q, **kwargs)**stoich
    return q


def _get_CvoR_rxn(initial_state, initial_state_stoich, final_state,
                  final_state_stoich, **kwargs):
    """Helper function to calculate dimensionless heat capacity

    Parameters
    ----------
        initial_state : list of pMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of pMuTT model objects
            Final state
        final_state_stoich : list of float
            Stoichiometry of final state
        kwargs : keyword arguments
            Parameters to calculate CvoR
    Returns
    -------
        CvoR : float
            Dimensionless heat capacity between initial state and final state
    """
    CvoR = 0.
    for specie, stoich in zip(final_state, final_state_stoich):
        CvoR += _force_pass_arguments(specie.get_CvoR, **kwargs)*stoich
    for specie, stoich in zip(initial_state, initial_state_stoich):
        CvoR -= _force_pass_arguments(specie.get_CvoR, **kwargs)*stoich
    return CvoR


def _get_CpoR_rxn(initial_state, initial_state_stoich, final_state,
                  final_state_stoich, **kwargs):
    """Helper function to calculate dimensionless heat capacity

    Parameters
    ----------
        initial_state : list of pMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of pMuTT model objects
            Final state
        final_state_stoich : list of float
            Stoichiometry of final state
        kwargs : keyword arguments
            Parameters to calculate CpoR
    Returns
    -------
        CpoR : float
            Dimensionless heat capacity between initial state and final state
    """
    CpoR = 0.
    for specie, stoich in zip(final_state, final_state_stoich):
        CpoR += _force_pass_arguments(specie.get_CpoR, **kwargs)*stoich
    for specie, stoich in zip(initial_state, initial_state_stoich):
        CpoR -= _force_pass_arguments(specie.get_CpoR, **kwargs)*stoich
    return CpoR


def _get_UoRT_rxn(initial_state, initial_state_stoich, final_state,
                  final_state_stoich, **kwargs):
    """Helper function to calculate dimensionless internal energy

    Parameters
    ----------
        initial_state : list of pMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of pMuTT model objects
            Final state
        final_state_stoich : list of float
            Stoichiometry of final state
        kwargs : keyword arguments
            Parameters to calculate UoRT
    Returns
    -------
        UoRT : float
            Dimensionless internal energy between initial state and final state
    """
    UoRT = 0.
    for specie, stoich in zip(final_state, final_state_stoich):
        UoRT += _force_pass_arguments(specie.get_UoRT, **kwargs)*stoich
    for specie, stoich in zip(initial_state, initial_state_stoich):
        UoRT -= _force_pass_arguments(specie.get_UoRT, **kwargs)*stoich
    return UoRT


def _get_HoRT_rxn(initial_state, initial_state_stoich, final_state,
                  final_state_stoich, **kwargs):
    """Helper function to calculate dimensionless enthalpy

    Parameters
    ----------
        initial_state : list of pMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of pMuTT model objects
            Final state
        final_state_stoich : list of float
            Stoichiometry of final state
        kwargs : keyword arguments
            Parameters to calculate HoRT
    Returns
    -------
        HoRT : float
            Dimensionless enthalpy between initial state and final state
    """
    HoRT = 0.
    for specie, stoich in zip(final_state, final_state_stoich):
        HoRT += _force_pass_arguments(specie.get_HoRT, **kwargs)*stoich
    for specie, stoich in zip(initial_state, initial_state_stoich):
        HoRT -= _force_pass_arguments(specie.get_HoRT, **kwargs)*stoich
    return HoRT


def _get_SoR_rxn(initial_state, initial_state_stoich, final_state,
                 final_state_stoich, **kwargs):
    """Helper function to calculate dimensionless entropy

    Parameters
    ----------
        initial_state : list of pMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of pMuTT model objects
            Final state
        final_state_stoich : list of float
            Stoichiometry of final state
        kwargs : keyword arguments
            Parameters to calculate SoR
    Returns
    -------
        SoR : float
            Dimensionless entropy between initial state and final state
    """
    SoR = 0.
    for specie, stoich in zip(final_state, final_state_stoich):
        SoR += _force_pass_arguments(specie.get_SoR, **kwargs)*stoich
    for specie, stoich in zip(initial_state, initial_state_stoich):
        SoR -= _force_pass_arguments(specie.get_SoR, **kwargs)*stoich
    return SoR


def _get_AoRT_rxn(initial_state, initial_state_stoich, final_state,
                  final_state_stoich, **kwargs):
    """Helper function to calculate dimensionless Helholtz energy

    Parameters
    ----------
        initial_state : list of pMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of pMuTT model objects
            Final state
        final_state_stoich : list of float
            Stoichiometry of final state
        kwargs : keyword arguments
            Parameters to calculate AoRT
    Returns
    -------
        AoRT : float
            Dimensionless Helmholtz energy between initial state and
            final state
    """
    AoRT = 0.
    for specie, stoich in zip(final_state, final_state_stoich):
        AoRT += _force_pass_arguments(specie.get_AoRT, **kwargs)*stoich
    for specie, stoich in zip(initial_state, initial_state_stoich):
        AoRT -= _force_pass_arguments(specie.get_AoRT, **kwargs)*stoich
    return AoRT


def _get_GoRT_rxn(initial_state, initial_state_stoich, final_state,
                  final_state_stoich, **kwargs):
    """Helper function to calculate dimensionless Gibbs energy

    Parameters
    ----------
        initial_state : list of pMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of pMuTT model objects
            Final state
        final_state_stoich : list of float
            Stoichiometry of final state
        kwargs : keyword arguments
            Parameters to calculate GoRT
    Returns
    -------
        GoRT : float
            Dimensionless Gibbs energy between initial state and final state
    """
    GoRT = 0.
    for specie, stoich in zip(final_state, final_state_stoich):
        GoRT += _force_pass_arguments(specie.get_GoRT, **kwargs)*stoich
    for specie, stoich in zip(initial_state, initial_state_stoich):
        GoRT -= _force_pass_arguments(specie.get_GoRT, **kwargs)*stoich
    return GoRT


def _parse_reaction_side(reaction_str, species_delimiter='+'):
    """Takes the reactants/products side of a reaction string and parse it
    into species and stoichiometric amounts

    Parameters
    ----------
        reaction_str : str
            Reactant or product side of reaction
        species_delimiters : str
            Delimiter that separate species. Leading and trailing spaces will
            be trimmed. Default is '+'

    Returns
    -------
        species : list of str
            Names of the species
        stoichiometry : list of int
            Stoichiometric coefficients
    """
    species_str = reaction_str.split(species_delimiter)
    species = []
    stoichiometry = []
    for specie in species_str:
        # Strip spaces for easier searching
        specie = specie.strip()
        # Search for int and float at the start of a string.
        # If there is no numbers, returns None.
        stoich_search = re.search(r'^\d+\.?\d*', specie)
        if stoich_search is None:
            # No stoichiometric coefficient so assign 1.
            species.append(specie.strip())
            stoichiometry.append(1.)
        else:
            # Stoichiometric coefficient present
            specie_stoich = stoich_search.group()
            trim_len = len(specie_stoich)
            species.append(specie[trim_len:].strip())
            stoichiometry.append(float(specie_stoich))
    return (species, stoichiometry)


def _parse_reaction(reaction_str, species_delimiter='+',
                    reaction_delimiter='='):
    """Takes a reaction string and parses it into reactants and products.

    Parameters
    ----------
        reaction_str : str
            Reaction string. A transition state can be specified by using two
            reaction delimiters.
            e.g. H2 + 0.5O2 = H2O_TS = H2O
        species_delimiter : str, optional
            Delimiter that separate species. Leading and trailing spaces will
            be trimmed. Default is '+'
        reaction_delimiter : str, optional
            Delimiter that separate sides of the reaction. Leading and trailing
            spaces will be trimmed. Default is '='
    Returns
    -------
        reactants : list of str
            Reactant names
        reactants_stoich : list of float
            Stoichiometry of reactants
        products : list of str
            Product names
        products_stoich : list of float
            Stoichiometry of products
        transition_state : list of str
            Transition state names. Returns None if reaction does not have a
            transition state.
        transition_state_stoich : list of float
            Stoichiometry of transition state. Returns None if the reaction
            does not have a transition state.
    """
    # Separate sides of reaction
    reaction_sides = reaction_str.split(reaction_delimiter)

    reactants_side = reaction_sides[0]
    products_side = reaction_sides[-1]

    # Separate the species in each side
    reactants, reactants_stoich = _parse_reaction_side(
            reaction_str=reactants_side, species_delimiter=species_delimiter)
    products, products_stoich = _parse_reaction_side(
            reaction_str=products_side, species_delimiter=species_delimiter)

    # Check for transition state
    if len(reaction_sides) > 2:
        transition_state_side = reaction_sides[1]
        transition_state, transition_state_stoich = _parse_reaction_side(
                reaction_str=transition_state_side,
                species_delimiter=species_delimiter)
    else:
        transition_state = None
        transition_state_stoich = None

    return (reactants, reactants_stoich, products, products_stoich,
            transition_state, transition_state_stoich)


def _write_reaction_side(species, stoich, species_delimiter='+'):
    """Writes one section of the reaction string

    Parameters
    ----------
        species : list of ``pMuTT.model`` objects
            Species to write
        stoich : list of float
            Stoichiometry corresponding to species
        species_delimiter : str, optional
            Delimiter that separates species. Default is '+'
    Returns
    -------
        reaction_str : str
            One section of the reaction string
    """
    for i, (specie, stoich_val) in enumerate(zip(species, stoich)):
        if i == 0:
            reaction_str = '{}{}'.format(stoich_val, specie.name)
        else:
            reaction_str += '{}{}{}'.format(species_delimiter, stoich_val,
                                            specie.name)
    return reaction_str


def _count_elements(species, stoich):
    """Total the number of elements

    Parameters
    ----------
        species : list of pMuTT model objects
            Species to count the elements
        stoich : list of float
            Stoichiometric coefficients of each element
    Returns
    -------
        element_count : collections.Counter object
            Sum of elements
    """
    element_count = Counter()
    for specie, stoich_specie in zip(species, stoich):
        for element, coeff in specie.elements.items():
            element_count += Counter({element: coeff*stoich_specie})
    return element_count
