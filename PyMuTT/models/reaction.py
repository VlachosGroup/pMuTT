# -*- coding: utf-8 -*-
from collections import Counter
import numpy as np
from PyMuTT import _force_pass_arguments
from PyMuTT.io_.jsonio import json_to_PyMuTT, remove_class

class Reaction:
    """Represents a chemical reaction

    Attributes
    ----------
        reactants : list of PyMuTT model objects
            Reactants
        reactants_stoich : list of float
            Stoichiometric quantities of reactants
        products : list of PyMuTT model objects
            Products
        products_stoich : list of float
            Stoichiometric quantities of products
        transition_state : PyMuTT model object
            Transition state specie, optional
    """

    def __init__(self, reactants, reactants_stoich, products, products_stoich,
                 transition_state=None):
        self.reactants = reactants
        self.reactants_stoich = reactants_stoich
        self.products = products
        self.products_stoich = products_stoich
        self.transition_state = transition_state


    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict


    def count_elements(self, species, stoich):
        """Total the number of elements

        Parameters
        ----------
            species : list of PyMuTT model objects
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


    def check_element_balance(self):
        """Checks the reactants, products and transition state elemental 
        composition

        Raises
        ------
            ValueError
                Raised if the reactants, products and/or transition state 
                element composition does not agree.
        """
        reactant_elements = self.count_elements(self.reactants, 
                                                self.reactants_stoich)
        product_elements = self.count_elements(self.products, 
                                               self.products_stoich)
        if reactant_elements != product_elements:
            raise ValueError('Number of elements in reactants and products do '
                             'not agree.\nReactant count: {}\n'
                             'Product count: {}'.format(reactant_elements, \
                                                        product_elements))

        if self.transition_state is not None:
            TS_elements = self.count_elements((self.transition_state,), (1.,))
            if reactant_elements != TS_elements:
                raise ValueError('Number of elements in reactants and '
                                 'transition state do not agree.\n'
                                 'Reactant count: {}\n'
                                 'Product count: {}'.format(reactant_elements, \
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
        delta_UoRT =  _get_UoRT_rxn(initial_state=self.reactants, 
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
        delta_HoRT =  _get_HoRT_rxn(initial_state=self.reactants, 
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
        delta_SoR =  _get_SoR_rxn(initial_state=self.reactants, 
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
        delta_GoRT =  _get_GoRT_rxn(initial_state=self.reactants, 
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
        """Gets change in partition function between reactants (or products) and 
        transition state

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
                              final_state=(self.transition_state,), 
                              final_state_stoich=(1.,),
                              **kwargs)
        else:
            return _get_q_rxn(initial_state=self.reactants, 
                              initial_state_stoich=self.reactants_stoich,
                              final_state=(self.transition_state,), 
                              final_state_stoich=(1.,),
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
                                 final_state=(self.transition_state,), 
                                 final_state_stoich=(1.,),
                                 **kwargs)
        else:
            return _get_CvoR_rxn(initial_state=self.reactants, 
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=(self.transition_state,), 
                                 final_state_stoich=(1.,),
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
                                 final_state=(self.transition_state,), 
                                 final_state_stoich=(1.,),
                                 **kwargs)
        else:
            return _get_CpoR_rxn(initial_state=self.reactants, 
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=(self.transition_state,), 
                                 final_state_stoich=(1.,),
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
                                 final_state=(self.transition_state,), 
                                 final_state_stoich=(1.,),
                                 **kwargs)
        else:
            return _get_UoRT_rxn(initial_state=self.reactants, 
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=(self.transition_state,), 
                                 final_state_stoich=(1.,),
                                 **kwargs)


    def get_HoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless enthalpy between reactants (or products) 
        and transition state

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
                                 final_state=(self.transition_state,), 
                                 final_state_stoich=(1.,),
                                 **kwargs)
        else:
            return _get_HoRT_rxn(initial_state=self.reactants, 
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=(self.transition_state,), 
                                 final_state_stoich=(1.,),
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
                Change in dimensionless entropy between reactants and transition
                state
        """
        if rev:
            return _get_SoR_rxn(initial_state=self.products, 
                                initial_state_stoich=self.products_stoich,
                                final_state=(self.transition_state,), 
                                final_state_stoich=(1.,),
                                **kwargs)
        else:
            return _get_SoR_rxn(initial_state=self.reactants, 
                                initial_state_stoich=self.reactants_stoich,
                                final_state=(self.transition_state,), 
                                final_state_stoich=(1.,),
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
                                 final_state=(self.transition_state,), 
                                 final_state_stoich=(1.,),
                                 **kwargs)
        else:
            return _get_AoRT_rxn(initial_state=self.reactants, 
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=(self.transition_state,), 
                                 final_state_stoich=(1.,),
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
                                 final_state=(self.transition_state,), 
                                 final_state_stoich=(1.,),
                                 **kwargs)
        else:
            return _get_GoRT_rxn(initial_state=self.reactants, 
                                 initial_state_stoich=self.reactants_stoich,
                                 final_state=(self.transition_state,), 
                                 final_state_stoich=(1.,),
                                 **kwargs)

    @classmethod
    def from_string(cls, reaction_string, species):
        pass

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
        json_obj['reactants'] = [json_to_PyMuTT(reactant) 
                                 for reactant in json_obj['reactants']]
        json_obj['products'] = [json_to_PyMuTT(product) 
                                 for product in json_obj['products']]
        json_obj['transition_state'] = json_to_PyMuTT(
                json_obj['transition_state'])
        return cls(**json_obj)

def _get_q_rxn(initial_state, initial_state_stoich, final_state, 
               final_state_stoich, **kwargs):
    """Helper function to calculate partition function

    Parameters
    ----------
        initial_state : list of PyMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of PyMuTT model objects
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
        initial_state : list of PyMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of PyMuTT model objects
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
        initial_state : list of PyMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of PyMuTT model objects
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
        initial_state : list of PyMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of PyMuTT model objects
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
        initial_state : list of PyMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of PyMuTT model objects
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
        initial_state : list of PyMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of PyMuTT model objects
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
        initial_state : list of PyMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of PyMuTT model objects
            Final state
        final_state_stoich : list of float
            Stoichiometry of final state
        kwargs : keyword arguments
            Parameters to calculate AoRT
    Returns
    -------
        AoRT : float
            Dimensionless Helmholtz energy between initial state and final state
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
        initial_state : list of PyMuTT model objects
            Initial state
        initial_state_stoich : list of float
            Stoichiometry of initial state
        final_state : list of PyMuTT model objects
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