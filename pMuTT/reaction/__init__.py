# -*- coding: utf-8 -*-
from collections import Counter
import inspect
import re
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
from pMuTT import (_force_pass_arguments, _pass_expected_arguments, 
    _is_iterable, _get_specie_kwargs, _apply_numpy_operation)
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
        bep : :class:`~pMuTT.reaction.bep.BEP` object, optional
            Bronsted Evans Polyani relationship. Default is None
        kwargs : keyword arguments
            BEP parameters. See :class:`~pMuTT.reaction.bep.BEP` documentation
            for expected parameters

    Notes
    -----
        Specie specific parameters can be passed by having a key named
        'specie' mapping onto a dictionary whose keys are the species names.

        e.g. For the reaction: H2 + 0.5O2 = H2O, the pressures can be specified
        independently using the following dictionary. In this example arbitrary
        values were used so the quantity evaluated may be meaningless.

        .. code:: python

            kwargs = {
                'T': 298.,
                'specie': {
                    'H2': {
                        'P': 2.,
                    },
                    'O2': {
                        'P': 1.,
                    },
                    'H2O': {
                        'P': 1.,
                    }
                }
            }
    """

    def __init__(self, reactants, reactants_stoich, products, products_stoich,
                 transition_state=None, transition_state_stoich=None,
                 bep=None, **kwargs):
        # If any of the entries were not iterable, assign them to a list
        if not _is_iterable(reactants) and reactants is not None:
            reactants = [reactants]
        if not _is_iterable(reactants_stoich) and reactants_stoich is not None:
            reactants_stoich = [reactants_stoich]
        if not _is_iterable(products) and products is not None:
            products = [products]
        if not _is_iterable(products_stoich) and products_stoich is not None:
            products_stoich = [products_stoich]
        if not _is_iterable(transition_state) and transition_state is not None:
            transition_state = [transition_state]
        if not _is_iterable(transition_state_stoich) \
           and transition_state_stoich is not None:
            transition_state_stoich = [transition_state_stoich]
        self.reactants = reactants
        self.reactants_stoich = reactants_stoich
        self.products = products
        self.products_stoich = products_stoich
        self.transition_state = transition_state
        self.transition_state_stoich = transition_state_stoich
        if inspect.isclass(bep):
            kwargs['reaction'] = self
            self.bep = _pass_expected_arguments(bep, **kwargs)
        else:
            self.bep = bep

    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict

    def __str__(self):
        return self.to_string()

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

        if self.transition_state is not None:
            TS_elements = _count_elements(self.transition_state,
                                          self.transition_state_stoich)
            if reactant_elements != TS_elements:
                raise ValueError('Number of elements in reactants and '
                                 'transition state do not agree.\n'
                                 'Reactant count: {}\n'
                                 'Product count: {}'.format(reactant_elements,
                                                            TS_elements))

    def get_species(self, include_TS=True, key='name'):
        """Returns the unique species included in the reaction.

        Parameters
        ----------
            include_TS : bool, optional
                Whether transition states should be included. Default is True
            key : str, optional
                Attribute to use as the key in the output dictionary. Default
                is name
        Returns
        -------
            species : dict
                Unique species in the reaction
        """
        species = {}
        # Add reactants
        for specie in self.reactants:
            species[getattr(specie, key)] = specie
        # Add products
        for specie in self.products:
            species[getattr(specie, key)] = specie
        # Add transition state if desired
        if include_TS and self.transition_state is not None:
            for specie in self.transition_state:
                species[getattr(specie, key)] = specie

        return species

    def get_q_state(self, state, **kwargs):
        """Gets partition function at a state

        Parameters
        ----------
            state : str
                state to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts'
            kwargs : keyword arguments
                Parameters with the conditions to calculate partition function.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            q : float
                Partition function of the reaction state
        """
        return self.get_state_quantity(state=state, method_name='get_q',
                                       **kwargs)

    def get_CvoR_state(self, state, **kwargs):
        """Gets dimensionless heat capacity at constant volume at a state

        Parameters
        ----------
            state : str
                state to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            kwargs : keyword arguments
                Parameters required to calculate dimensionless heat capacity at
                constant volume. See class docstring to see how to pass
                specific parameters to different species.
        Returns
        -------
            CvoR : float
                Dimensionless heat capacity at constant volume of the reaction
                state
        """
        return self.get_state_quantity(state=state, method_name='get_CvoR',
                                       **kwargs)

    def get_Cv_state(self, state, units, **kwargs):
        """Gets the heat capacity at constant volume at a state

        Parameters
        ----------
            state : str
                State to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            kwargs : keyword arguments
                Parameters required to calculate heat capacity at constant
                volume. See class docstring to see how to pass specific
                parameters to different species.
        Returns
        -------
            Cv : float
                Heat capacity at constant volume of the reaction state
        """
        return self.get_CvoR_state(state=state, **kwargs)*c.R(units)

    def get_CpoR_state(self, state, **kwargs):
        """Gets dimensionless heat capacity at constant pressure at a state

        Parameters
        ----------
            state : str
                state to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            kwargs : keyword arguments
                Parameters required to calculate dimensionless heat capacity at
                constant pressure. See class docstring to see how to pass
                specific parameters to different species.

        Returns
        -------
            CpoR : float
                Dimensionless heat capacity at constant pressure
                of the reaction state
        """
        return self.get_state_quantity(state=state, method_name='get_CpoR',
                                       **kwargs)

    def get_Cp_state(self, state, units, **kwargs):
        """Gets the heat capacity at constant pressure at a state

        Parameters
        ----------
            state : str
                State to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            kwargs : keyword arguments
                Parameters required to calculate heat capacity at constant
                pressure. See class docstring to see how to pass specific
                parameters to different species.
        Returns
        -------
            Cp : float
                Heat capacity at constant pressure of the reaction state
        """
        return self.get_CpoR_state(state=state, **kwargs)*c.R(units)

    def get_UoRT_state(self, state, **kwargs):
        """Gets dimensionless internal energy at a state

        Parameters
        ----------
            state : str
                state to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            kwargs : keyword arguments
                Parameters required to calculate dimensionless internal energy.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            UoRT : float
                Dimensionless internal energy of the reaction state.

        """
        return self.get_state_quantity(state=state, method_name='get_UoRT',
                                       **kwargs)

    def get_U_state(self, state, units, T, **kwargs):
        """Gets the internal energy at a state

        Parameters
        ----------
            state : str
                State to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            kwargs : keyword arguments
                Parameters required to calculate internal energy.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            U : float
                Internal energy of the reaction state
        """
        return self.get_UoRT_state(state=state, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

    def get_EoRT_state(self, state, include_ZPE=False, **kwargs):
        """Gets dimensionless electronic energy at a state

        Parameters
        ----------
            state : str
                state to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            include_ZPE : bool, optional
                If True, includes the zero point energy. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless electronic
                energy. See class docstring to see how to pass specific
                parameters to different species.
        Returns
        -------
            EoRT : float
                Dimensionless electronic energy of the reaction state.
        """
        return self.get_state_quantity(state=state, method_name='get_EoRT',
                                       include_ZPE=include_ZPE, **kwargs)

    def get_E_state(self, state, units, T=c.T0('K'), include_ZPE=False,
                    **kwargs):
        """Gets the electronic energy at a state

        Parameters
        ----------
            state : str
                State to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            include_ZPE : bool, optional
                If True, includes the zero point energy. Default is False
            kwargs : keyword arguments
                Parameters required to calculate electronic energy.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            E : float
                Electronic energy of the reaction state
        """
        return self.get_EoRT_state(state=state, T=T, include_ZPE=include_ZPE, 
                                   **kwargs) \
               *T*c.R('{}/K'.format(units))

    def get_HoRT_state(self, state, **kwargs):
        """Gets dimensionless enthalpy at a state

        Parameters
        ----------
            state : str
                state to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            kwargs : keyword arguments
                Parameters required to calculate dimensionless enthalpy. See
                class docstring to see how to pass specific parameters to
                different species.

        Returns
        -------
            HoRT : float
                Dimensionless heat capacity at constant pressure
                of the reaction state.
        """
        return self.get_state_quantity(state=state, method_name='get_HoRT',
                                       **kwargs)

    def get_H_state(self, state, units, T, **kwargs):
        """Gets the enthalpy at a state

        Parameters
        ----------
            state : str
                State to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            kwargs : keyword arguments
                Parameters required to calculate enthalpy.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            H : float
                Enthalpy of the reaction state
        """
        return self.get_HoRT_state(state=state, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

    def get_SoR_state(self, state, **kwargs):
        """Gets dimensionless entropy at a state

        Parameters
        ----------
            state : str
                state to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            kwargs : keyword arguments
                Parameters required to calculate dimensionless entropy. See
                class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            SoR : float
                Dimensionless entropy of the reaction state
        """
        return self.get_state_quantity(state=state, method_name='get_SoR',
                                       **kwargs)

    def get_S_state(self, state, units, **kwargs):
        """Gets the entropy at a state

        Parameters
        ----------
            state : str
                State to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            kwargs : keyword arguments
                Parameters required to calculate entropy.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            S : float
                Entropy of the reaction state
        """
        return self.get_SoR_state(state=state, **kwargs)*c.R(units)

    def get_FoRT_state(self, state, **kwargs):
        """Gets dimensionless Helmholtz energy at a state

        Parameters
        ----------
            state : str
                state to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            kwargs : keyword arguments
                Parameters required to calculate dimensionless
                Helmholtz energy.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            FoRT : float
                Dimensionless Helmoltz energy of the reaction state
        """
        return self.get_state_quantity(state=state, method_name='get_FoRT',
                                       **kwargs)

    def get_F_state(self, state, units, T, **kwargs):
        """Gets the Helholtz energy at a state

        Parameters
        ----------
            state : str
                State to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            kwargs : keyword arguments
                Parameters required to calculate Helholtz energy.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            F : float
                Helmholtz energy of the reaction state
        """
        return self.get_FoRT_state(state=state, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

    def get_GoRT_state(self, state, **kwargs):
        """Gets dimensionless Gibbs energy at a state

        Parameters
        ----------
            state : str
                state to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            kwargs : keyword arguments
                Parameters required to calculate dimensionless Gibbs energy.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            GoRT : float
                Dimensionless Gibbs energy of the reaction state
        """
        return self.get_state_quantity(state=state, method_name='get_GoRT',
                                       **kwargs)

    def get_G_state(self, state, units, T, **kwargs):
        """Gets the Gibbs energy at a state

        Parameters
        ----------
            state : str
                State to calculate quantity. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
                - 'ts' (same as transition state)
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            kwargs : keyword arguments
                Parameters required to calculate Gibbs energy.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            G : float
                Gibbs energy of the reaction state
        """
        return self.get_GoRT_state(state=state, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

    def get_delta_q(self, rev=False, act=False, **kwargs):
        """Gets change in partition function between reactants and products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate partition function. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_q : float
                Change in partition function between reactants and products
        """
        initial_state, final_state = _get_states(rev=rev, act=act)
        delta_q = self.get_delta_quantity(initial_state=initial_state,
                                          final_state=final_state,
                                          method_name='get_q', **kwargs)
        return delta_q

    def get_delta_CvoR(self, rev=False, act=False, **kwargs):
        """Gets change in dimensionless heat capacity (constant V)
        between reactants and products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate heat capacity. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_CvoR : float
                Change in heat capacity between reactants and products
        """
        initial_state, final_state = _get_states(rev=rev, act=act)
        delta_CvoR = self.get_delta_quantity(initial_state=initial_state,
                                             final_state=final_state,
                                             method_name='get_CvoR',
                                             **kwargs)
        return delta_CvoR

    def get_delta_Cv(self, units, rev=False, act=False, **kwargs):
        """Gets change in heat capacity (constant V) between reactants and
        products

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate heat capacity. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_Cv : float
                Change in heat capacity between reactants and products
        """
        return self.get_delta_CvoR(rev=rev, act=act, **kwargs) \
               *c.R(units)

    def get_delta_CpoR(self, rev=False, act=False, **kwargs):
        """Gets change in dimensionless heat capacity between reactants and
        products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate heat capacity. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_CvoR : float
                Change in heat capacity between reactants and products
        """
        initial_state, final_state = _get_states(rev=rev, act=act)
        delta_CpoR = self.get_delta_quantity(initial_state=initial_state,
                                             final_state=final_state,
                                             method_name='get_CpoR',
                                             **kwargs)
        return delta_CpoR

    def get_delta_Cp(self, units, rev=False, act=False, **kwargs):
        """Gets change in heat capacity (constant P) between reactants and
        products

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate heat capacity. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_Cp : float
                Change in heat capacity between reactants and products
        """
        return self.get_delta_CpoR(rev=rev, act=act, **kwargs) \
               *c.R(units)

    def get_delta_UoRT(self, rev=False, act=False, **kwargs):
        """Gets change in dimensionless internal energy between reactants and
        products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate internal energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_UoRT : float
                Change in internal energy between reactants and products
        """
        initial_state, final_state = _get_states(rev=rev, act=act)
        delta_UoRT = self.get_delta_quantity(initial_state=initial_state,
                                             final_state=final_state,
                                             method_name='get_UoRT',
                                             **kwargs)
        return delta_UoRT

    def get_delta_U(self, units, T, rev=False, act=False, **kwargs):
        """Gets change in internal energy between reactants and products

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate internal energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_U : float
                Change in internal energy between reactants and products
        """
        return self.get_delta_UoRT(rev=rev, T=T, act=act,
                                   **kwargs)*T*c.R('{}/K'.format(units))

    def get_delta_EoRT(self, rev=False, act=False, **kwargs):
        """Gets change in dimensionless electronic energy between reactants and
        products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate electronic energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_EoRT : float
                Change in electronic energy between reactants and products
        """
        initial_state, final_state = _get_states(rev=rev, act=act)
        delta_EoRT = self.get_delta_quantity(initial_state=initial_state,
                                             final_state=final_state,
                                             method_name='get_EoRT',
                                             **kwargs)
        return delta_EoRT

    def get_delta_E(self, units, T, rev=False, act=False, **kwargs):
        """Gets change in electronic energy between reactants and products

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate electronic energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_E : float
                Change in electronic energy between reactants and products
        """
        return self.get_delta_EoRT(rev=rev, act=act, T=T, 
                                   **kwargs)*T*c.R('{}/K'.format(units))

    def get_delta_HoRT(self, rev=False, act=False, **kwargs):
        """Gets change in dimensionless enthalpy between reactants and products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate enthalpy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_HoRT : float
                Change in enthalpy between reactants and products
        """
        initial_state, final_state = _get_states(rev=rev, act=act)
        delta_HoRT = self.get_delta_quantity(initial_state=initial_state,
                                             final_state=final_state,
                                             method_name='get_HoRT',
                                             **kwargs)
        return delta_HoRT

    def get_delta_H(self, units, T, rev=False, act=False, **kwargs):
        """Gets change in enthalpy between reactants and products

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate enthalpy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_H : float
                Change in enthalpy between reactants and products
        """
        return self.get_delta_HoRT(rev=rev, T=T, act=act, 
                                   **kwargs)*T*c.R('{}/K'.format(units))

    def get_delta_SoR(self, rev=False, act=False, **kwargs):
        """Gets change in dimensionless entropy between reactants and products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate entropy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_SoR : float
                Change in entropy between reactants and products
        """
        initial_state, final_state = _get_states(rev=rev, act=act)
        delta_SoR = self.get_delta_quantity(initial_state=initial_state,
                                            final_state=final_state,
                                            method_name='get_SoR',
                                            **kwargs)
        return delta_SoR

    def get_delta_S(self, units, rev=False, act=False, **kwargs):
        """Gets change in entropy between reactants and products

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate entropy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_S : float
                Change in entropy between reactants and products
        """
        return self.get_delta_SoR(rev=rev, act=act, **kwargs) \
               *c.R(units)

    def get_delta_FoRT(self, rev=False, act=False, **kwargs):
        """Gets change in dimensionless Helmholtz energy between reactants and
        products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate Helmholtz energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_FoRT : float
                Change in Helmholtz energy between reactants and products
        """
        initial_state, final_state = _get_states(rev=rev, act=act)
        delta_FoRT = self.get_delta_quantity(initial_state=initial_state,
                                             final_state=final_state,
                                             method_name='get_FoRT',
                                             **kwargs)
        return delta_FoRT

    def get_delta_F(self, units, T, rev=False, act=False, **kwargs):
        """Gets change in Helmholtz energy between reactants and products

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate Helmholtz energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_F : float
                Change in Helmholtz energy between reactants and products
        """
        return self.get_delta_FoRT(rev=rev, T=T, act=act, 
                                   **kwargs)*T*c.R('{}/K'.format(units))

    def get_delta_GoRT(self, rev=False, act=False, **kwargs):
        """Gets change in dimensionless Gibbs energy between reactants and
        products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate Gibbs energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_GoRT : float
                Change in Gibbs energy between reactants and products
        """
        initial_state, final_state = _get_states(rev=rev, act=act)
        delta_GoRT = self.get_delta_quantity(initial_state=initial_state,
                                             final_state=final_state,
                                             method_name='get_GoRT',
                                             **kwargs)
        return delta_GoRT

    def get_delta_G(self, units, T, rev=False, act=False, **kwargs):
        """Gets change in Gibbs energy between reactants and products

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate Gibbs energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_G : float
                Change in Gibbs energy between reactants and products
        """
        return self.get_delta_GoRT(rev=rev, T=T, act=act,
                                   **kwargs)*T*c.R('{}/K'.format(units))

    def get_q_act(self, rev=False, **kwargs):
        """Gets change in partition function between reactants/products and the 
        transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate partition function. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            q_act : float
                Change in partition function between reactants/products and the
                transition state
        """
        return self.get_delta_q(rev=rev, act=True, **kwargs)

    def get_CvoR_act(self, rev=False, **kwargs):
        """Gets change in dimensionless heat capacity (constant V)
        between reactants/products and the transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate heat capacity. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            CvoR_act : float
                Change in heat capacity between reactants/products and
                transition state
        """
        return self.get_delta_CvoR(rev=rev, act=True, **kwargs)

    def get_Cv_act(self, units, rev=False, **kwargs):
        """Gets change in heat capacity (constant V) between reactants/products
        and the transition state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate heat capacity. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            Cv_act : float
                Change in heat capacity between reactants and products
        """
        return self.get_delta_Cv(units=units, rev=rev, act=True, **kwargs)

    def get_CpoR_act(self, rev=False, **kwargs):
        """Gets change in dimensionless heat capacity between reactants/products
        and the transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate heat capacity. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            CvoR_act : float
                Change in heat capacity between reactants/products and the
                transition state
        """
        return self.get_delta_CpoR(rev=rev, act=True, **kwargs)

    def get_Cp_act(self, units, rev=False, **kwargs):
        """Gets change in heat capacity (constant P) between reactants/products
        and the transition state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate heat capacity. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_Cp : float
                Change in heat capacity between reactants/products and the
                transition state
        """
        return self.get_delta_Cp(units=units, rev=rev, act=True, **kwargs)

    def get_UoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless internal energy between 
        reactants/products and the transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate internal energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            UoRT_act : float
                Change in internal energy between reactants/products and the
                transition state
        """
        return self.get_delta_UoRT(rev=rev, act=True, **kwargs)

    def get_U_act(self, units, T, rev=False, **kwargs):
        """Gets change in internal energy between reactants/products and the
        transition state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate internal energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            U_act : float
                Change in internal energy between reactants/products and the 
                transition state
        """
        return self.get_delta_U(units=units, T=T, rev=rev, act=True, **kwargs)

    def get_HoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless enthalpy between reactants/products and
        the transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate enthalpy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            HoRT_act : float
                Change in enthalpy between reactants/products and the 
                transition state
        """
        return self.get_delta_HoRT(rev=rev, act=True, **kwargs)

    def get_H_act(self, units, T, rev=False, **kwargs):
        """Gets change in enthalpy between reactants/products and the transition
        state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate enthalpy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            H_act : float
                Change in enthalpy between reactants/products and the transition
                state
        """
        return self.get_delta_H(units=units, T=T, rev=rev, act=True, **kwargs)

    def get_SoR_act(self, rev=False, **kwargs):
        """Gets change in dimensionless entropy between reactants/products and
        the transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate entropy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            SoR_act : float
                Change in entropy between reactants/products and the transition
                state
        """
        return self.get_delta_SoR(rev=rev, act=True, **kwargs)

    def get_S_act(self, units, rev=False, **kwargs):
        """Gets change in entropy between reactants/products and the transition
        state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate entropy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            S_act : float
                Change in entropy between reactants/products and the transition
                state
        """
        return self.get_delta_S(units=units, rev=rev, act=True, **kwargs)

    def get_FoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless Helmholtz energy between
        reactants/products and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate Helmholtz energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            FoRT_act : float
                Change in Helmholtz energy between reactants/products and the 
                transition state
        """
        return self.get_delta_FoRT(rev=rev, act=True, **kwargs)

    def get_F_act(self, units, T, rev=False, **kwargs):
        """Gets change in Helmholtz energy between reactants/products and the
        transition state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate Helmholtz energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            F_act : float
                Change in Helmholtz energy between reactants/products and the
                transition state
        """
        return self.get_delta_F(units=units, T=T, rev=rev, act=True, **kwargs)

    def get_GoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless Gibbs energy between reactants/products
        and the transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate Gibbs energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            GoRT_act : float
                Change in Gibbs energy between reactants/products and the
                transition state
        """
        return self.get_delta_GoRT(rev=rev, act=True, **kwargs)

    def get_G_act(self, units, T, rev=False, **kwargs):
        """Gets change in Gibbs energy between reactants/products and the
        transition state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate Gibbs energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            G_act : float
                Change in Gibbs energy between reactants/products and the
                transition state
        """
        return self.get_delta_G(units=units, T=T, rev=rev, act=True, **kwargs)

    def get_Keq(self, rev=False, act=False, **kwargs):
        """Gets equilibrium constant between reactants and products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate equilibrium constant. See
                class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            Keq : float
                Equilibrium constant between reactants and products
        """
        return np.exp(-self.get_delta_GoRT(rev=rev, act=act,
                                           **kwargs))

    def get_EoRT_act(self, rev=False, method='any', **kwargs):
        """Gets dimensionless act energy between reactants
        (or products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            method : str, optional
                Method to use to calculate dimensionless act energy.
                Accepted options:

                - 'any' (uses whichever is available. First uses transition
                  state theory, then BEPs, then the reaction enthalpy)
                - 'bep' (uses ``self.bep``)
                - 'ts' or 'transition_state' (uses ``self.transition_state``)
                - 'enthalpy' (uses the enthalpy of the reaction. If the reaction
                  is exothermic, returns 0)

                Default is 'any'.
            kwargs : keyword arguments
                Parameters required to calculate dimensionless act
                energy
        Returns
        -------
            EoRT_act : float
                Dimensionless act energy between reactants (or products)
                and transition state
        """
        method = method.lower()
        # Assign the appropriate method is 'any' was specified
        if method == 'any':
            if self.transition_state is not None:
                method = 'transition_state'
            elif self.bep is not None:
                method = 'bep'
            else:
                method = 'enthalpy'
        
        if method == 'transition_state' or method == 'ts':
            EoRT = self.get_delta_HoRT(rev=rev, act=True, **kwargs)
        elif method == 'bep':
            EoRT = self.bep.get_EoRT_act(rev=rev, **kwargs)
        elif method == 'enthalpy':
            EoRT = np.max([0., self.get_delta_HoRT(rev=rev, act=False,
                                                   **kwargs)])
        else:
            raise ValueError(('Method "{}" not supported. See documentation '
                              'of '
                              '``pMuTT.reaction.Reaction.get_EoRT_act`` '
                              'for supported options.'.format(method)))
        return EoRT

    def get_E_act(self, units, T, rev=False, method='any', **kwargs):
        """Gets act energy between reactants (or products)
        and transition state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float
                Temperature in K
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            method : str, optional
                Method to use to calculate act energy.
                Accepted options:

                - 'any' (uses whichever is available. ``self.transition_state``
                  is used preferentially over ``self.BEP``)
                - 'bep' (uses ``self.bep``)
                - 'ts' or 'transition_state' (uses ``self.transition_state``)

                Default is 'any'.
            kwargs : keyword arguments
                Parameters required to calculate act energy. See class
                docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            E_act : float
                act energy between reactants (or products) and
                transition state
        """
        return self.get_EoRT_act(rev=rev, method=method, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

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
                Parameters required to calculate pre-exponential factor. See
                class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            A : float
                Pre-exponential factor
        """
        # Calculate molecularity (e.g. unimolecular, bimolecular)
        if rev:
            m = _get_molecularity(self.products_stoich)
        else:
            m = _get_molecularity(self.reactants_stoich)

        return c.kb('J/K')*T/c.h('J s') \
               *np.exp(self.get_delta_SoR(rev=rev, act=True, T=T, 
                                          **kwargs)+m)

    def _parse_state(self, state):
        """Helper method to get the relevant species and stoichiometry

        Parameters
        ----------
            state : str
                Thermodynamic state. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
        Returns
        -------
            species : list of ``pMuTT`` specie objects
                Species for the specified state
            species_stoich : list of float
                Stoichiometry corresponding to the species
        """
        state = state.lower()
        if state == 'reactants':
            species = self.reactants
            species_stoich = self.reactants_stoich
        elif state == 'products':
            species = self.products
            species_stoich = self.products_stoich
        elif state == 'transition state' or state == 'ts':
            species = self.transition_state
            species_stoich = self.transition_state_stoich
        else:
            raise ValueError(
                    'Thermodynamic state {} not supported'.format(state))
        return (species, species_stoich)

    def get_state_quantity(self, state, method_name, **kwargs):
        """Helper method to calculate the thermodynamic quantity of the state

        Parameters
        ----------
            state : str
                Thermodynamic state. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
            method_name : str
                Name of method to use to calculate quantity. Calculates any
                quantity as long as the relevant objects have the same method
                name
            kwargs : keyword arguments
                Arguments required to calculate the thermodynamic quantity of
                interest.
        Returns
        -------
            state_quantity : float
                Thermodynamic quantity of particular state
        """
        species, stoich = self._parse_state(state=state)

        if method_name == 'get_q':
            state_quantity = 1.
        else:
            state_quantity = 0.

        for specie, stoich in zip(species, stoich):
            # Process the inputs and methods for each specie
            specie_kwargs = _get_specie_kwargs(specie.name, **kwargs)
            try:
                specie_kwargs['specie'] = kwargs['specie']
            except KeyError:
                pass

            method = getattr(specie, method_name)
            if method_name == 'get_q':
                state_quantity *= \
                        _force_pass_arguments(method, **specie_kwargs)**stoich
            else:
                state_quantity += \
                        _force_pass_arguments(method, **specie_kwargs)*stoich
        return state_quantity

    def get_delta_quantity(self, initial_state, final_state, method_name,
                            **kwargs):
        """Helper method to calculate the change in thermodynamic quantity
        between states

        Parameters
        ----------
            initial_state : str
                Thermodynamic state. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
            final_state : str
                Thermodynamic state. Supported options:

                - 'reactants'
                - 'products'
                - 'transition state'
            method_name : str
                Name of method to use to calculate quantity. Calculates any
                quantity as long as the relevant objects have the same method
                name
            kwargs : keyword arguments
                Arguments passed to evaluate the quantity of the reactants and 
                products
        Returns
        -------
            delta_quantity : float
                Change in thermodynamic quantity between particular states
        """
        initial_quantity = self.get_state_quantity(state=initial_state,
                                                   method_name=method_name,
                                                   **kwargs)
        final_quantity = self.get_state_quantity(state=final_state,
                                                 method_name=method_name,
                                                 **kwargs)
        if method_name == 'get_q':
            return final_quantity/initial_quantity
        else:
            return final_quantity - initial_quantity

    @classmethod
    def from_string(cls, reaction_str, species, species_delimiter='+',
                    reaction_delimiter='=', **kwargs):
        """Create a reaction object using the reaction string

        Parameters
        ----------
            reaction_str : str
                Reaction string.
            species : dict
                Dictionary using the names as keys. If you have a list of
                species, use pMuTT.pMuTT_list_to_dict to make a dict.
            species_delimiter : str, optional
                Delimiter that separate species. Leading and trailing spaces
                will be trimmed. Default is '+'
            reaction_delimiter : str, optional
                Delimiter that separate states of the reaction. Leading and
                trailing spaces will be trimmed. Default is '='
            kwargs : keyword arguments
                BEP parameters. See :class:`~pMuTT.reaction.bep.BEP`
                documentation for expected parameters
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
                   transition_state=ts, transition_state_stoich=ts_stoich,
                   **kwargs)

    def to_string(self, species_delimiter='+', reaction_delimiter='=',
                  stoich_format='.2f', include_TS=True, key='name'):
        """Writes the Reaction object as a stoichiometric reaction

        Parameters
        ----------
            species_delimiter : str, optional
                Separates species. Default is '+'
            reaction_delimiter : str, optional
                Separates reaction states. Default is '='
            stoich_format : float, optional
                Format to write stoichiometric numbers. Default is '.2f' (float
                rounded to 2 decimal places)
            include_TS : bool, optional
                If True, includes transition states in output. Default is True
            key : str, optional
                Attribute to use to print out species. Default is name
        Returns
        -------
            reaction_str : str
                Reaction string
        """
        # Write reactants
        reaction_str = _write_reaction_state(species=self.reactants,
                                             stoich=self.reactants_stoich,
                                             species_delimiter=species_delimiter,
                                             stoich_format=stoich_format,
                                             key=key)
        reaction_str += reaction_delimiter

        # Write transition state if any
        if include_TS and self.transition_state is not None:
            reaction_str += _write_reaction_state(
                    species=self.transition_state,
                    stoich=self.transition_state_stoich,
                    species_delimiter=species_delimiter,
                    stoich_format=stoich_format,
                    key=key)
            reaction_str += reaction_delimiter

        # Write products
        reaction_str += _write_reaction_state(
                species=self.products, stoich=self.products_stoich,
                species_delimiter=species_delimiter,
                stoich_format=stoich_format,
                key=key)
        return reaction_str

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = {
            'class': str(self.__class__)}
        # Reactants
        if _is_iterable(self.reactants):
            obj_dict['reactants'] = \
                    [reactant.to_dict() for reactant in self.reactants]
        else:
            obj_dict['reactants'] = self.reactants
        # Reactants stoichiometry
        if _is_iterable(self.reactants_stoich):
            obj_dict['reactants_stoich'] = list(self.reactants_stoich)
        else:
            obj_dict['reactants_stoich'] = self.reactants_stoich
        # Products
        if _is_iterable(self.products):
            obj_dict['products'] = \
                    [product.to_dict() for product in self.products]
        else:
            obj_dict['products'] = self.products
        # Product stoichiometry
        if _is_iterable(self.products_stoich):
            obj_dict['products_stoich'] = list(self.products_stoich)
        else:
            obj_dict['products_stoich'] = self.products_stoich
        # Transition state
        if _is_iterable(self.transition_state):
            obj_dict['transition_state'] = \
                    [ts.to_dict() for ts in self.transition_state]
        else:
            obj_dict['transition_state'] = self.transition_state
        # Transition state stoich
        if _is_iterable(self.transition_state_stoich):
            obj_dict['transition_state_stoich'] = \
                    list(self.transition_state_stoich)
        else:
            obj_dict['transition_state_stoich'] = self.transition_state_stoich
        # BEP
        try:
            obj_dict['bep'] = self.bep.to_dict()
        except (AttributeError, TypeError):
            obj_dict['bep'] = self.bep

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
        json_obj['bep'] = json_to_pMuTT(json_obj['bep'])

        reaction = cls(**json_obj)
        if reaction.bep is not None:
            reaction.bep.set_descriptor(reaction=reaction)
        return reaction


class ChemkinReaction(Reaction):
    """Chemkin reaction. Has additional attributes to support input and output

    Attributes
    ----------
        beta : float, optional
            Power to raise the temperature in the rate expression. Default is 1
        is_adsorption : bool, optional
            If True, the reaction represents an adsorption. Default is False
        sticking_coeff : float, optional
            Sticking coefficient. Only relevant if ``is_adsorption`` is True
        gas_phase : bool
            True if the reaction has only gas-phase species. This attribute is
            determined based on the reactants and products
        kwargs : keyword arguments
            Keyword arguments used to initialize the reactants, transition state
            and products
    """

    def __init__(self, beta=1., is_adsorption=False, sticking_coeff=0.5,
                 **kwargs):
        super().__init__(**kwargs)
        self.beta = beta
        self.is_adsorption = is_adsorption
        # Sticking coefficient not relevant for non-adsorption reaction
        if not self.is_adsorption:
            sticking_coeff = None
        self.sticking_coeff = sticking_coeff
        self.gas_phase = self._is_gas_phase()

    def _is_gas_phase(self):
        """Determines if a reaction is gas phase
        
        Returns
        -------
            gas_phase : bool
               True if all the reactants are gas phase
        """
        return all([specie.phase == 'G' for specie in self.reactants])

    def _get_n_surf(self):
        """Counts the number of surface reactants

        Returns
        -------
            n_surf : int
                Number of surface species
        """
        n_surf = 0
        for specie, stoich in zip(self.reactants, self.reactants_stoich):
            # Skip species without catalyst site
            if specie.cat_site is None:
                continue
            # Skip bulk species
            if specie.cat_site.bulk_specie == specie.name:
                continue
            # Skip non-surface species
            if specie.phase != 'S':
                continue
            n_surf += stoich
        return n_surf

    def get_A(self, sden_operation='sum', include_entropy=True, T=c.T0('K'), 
              **kwargs):
        """Calculates the preexponential factor in the Chemkin format
        
        Parameters
        ----------
        sden_operation : str, optional
            Site density operation to use. Default is 'sum'
        include_entropy : bool, optional
            If True, includes the act entropy. Default is True
        T : float, optional
            Temperature in K. Default is 298.15 K
        kwargs : keyword arguments
            Parameters required to calculate pre-exponential factor
        """
        if self.transition_state is None or not include_entropy:
            A = c.kb('J/K')/c.h('J s')
        else:
            A = super().get_A(T=T, **kwargs)/T

        # If this is a surface reaction, adjust for site density
        if not self.gas_phase:
            # Uses site with highest site density
            site_dens = []
            for reactant in self.reactants:
                # Skip species without a catalyst site
                try:
                    site_den = reactant.cat_site.site_density
                except AttributeError:
                    continue
                # Skip bulk species
                if reactant.name == reactant.cat_site.bulk_specie:
                    continue
                site_dens.append(site_den)             

            # Apply the operation to the site densities
            eff_site_den = _apply_numpy_operation(quantity=site_dens,
                                                  operation=sden_operation,
                                                  verbose=False)
            n_surf = self._get_n_surf()
            A = A/eff_site_den**(n_surf-1)
        return A

    def get_HoRT_act(self, rev=False, **kwargs):
        """Calculates the dimensionless enthalpy. If there is no transition
        state species, calculates the delta dimensionless Gibbs energy

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate Gibbs energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            HoRT_act : float
                Change in Gibbs energy between reactants/products and the
                transition state
        """
        if self.transition_state is None:
            act = False
        else:
            act = True
        return np.max([0., super().get_delta_HoRT(rev=rev, act=act, **kwargs)])

    def get_GoRT_act(self, rev=False, act=False, **kwargs):
        """Calculates the dimensionless Gibbs energy. If there is no transition
        state species, calculates the delta dimensionless Gibbs energy

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate Gibbs energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            GoRT_act : float
                Change in Gibbs energy between reactants/products and the
                transition state
        """
        if self.transition_state is None:
            act = False
        else:
            act = True
        return np.max([0., super().get_delta_GoRT(rev=rev, act=act, **kwargs)])

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = super().to_dict()
        obj_dict['beta'] = self.beta
        obj_dict['is_adsorption'] = self.is_adsorption
        obj_dict['sticking_coeff'] = self.sticking_coeff
        obj_dict['gas_phase'] = self.gas_phase
        return obj_dict


class Reactions:
    """Contains multiple reactions. Serves as a parent class for other objects

    Attributes
    ----------
        reactions : list of :class:`~pMuTT.reaction.Reaction` objects
    """
    def __init__(self, reactions):
        self.reactions = list(reactions)

    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict

    def __iter__(self):
        for reaction in self.reactions:
            yield reaction

    def __getitem__(self, key):
        return self.reactions[key]

    def __len__(self):
        return len(self.reactions)

    def get_species(self, include_TS=True, key='name'):
        """Returns the unique species included in the reactions.

        Parameters
        ----------
            include_TS : bool, optional
                Whether transition states should be included. Default is True
            key : str, optional
                Attribute to use as the key in the output dictionary. Default is
                name
        Returns
        -------
            species : dict
                Unique species in the reactions
        """
        species = {}
        for reaction in self.reactions:
            # Merge the old dictionary with the new dictionary
            species.update(reaction.get_species(include_TS=include_TS,
                                                key=key))
        return species

    def plot_coordinate_diagram(self, method_name, ref_index=0,
                                ref_state='reactants', x_offset=1.,
                                include_TS=True, x_scale_TS=0.5, y_scale_TS=0.5, 
                                include_TS_labels=True, y_TS_label_offset=0.1, 
                                x_TS_label_offset=0., TS_label_format='.2f',
                                figure=None, axes=None, plt_kwargs={}, 
                                **reaction_kwargs):
        """Plots a reaction coordinate diagram.

        Parameters
        ----------
            method_name : str
                Name of method to use to calculate quantity. Calculates any
                quantity as long as the relevant objects have the same method
                name. Some examples include: ``get_HoRT``, ``get_H``,
                ``get_EoRT``, ``get_E``
            ref_index : int, optional
                Reaction index to use to reference states. Default is the first 
                reaction (i.e. ref_index = 0)
            ref_state : str, optional
                State of the reference to use. Supported options include:

                - reactants (default)
                - products
                - transition state
                - ts (same as transition state)
            x_offset : float, optional
                Spacing between reaction states. Shape of curve likely does not
                change with this parameter since the x axis would rescale
                appropriately
            include_TS : bool, optional
                Whether transition states should be included. Default is True
            x_scale_TS : float, optional
                Value between 0 and 1 that controls curvature of transition 
                state peaks. Higher values produce sharper peaks. Default is
                0.5
            y_scale_TS : float, optional
                Value between 0 and 1 that controls curvature of transition 
                state peaks. Higher values produce sharper peaks. Default is
                0.5
            include_TS_labels : bool, optional
                If True, adds a label to the peaks indicating the difference 
                between the reactants and the transition state. Default is True
            y_TS_label_offset : float, optional
                Vertical value to offset TS_label from the TS position. This
                value scales with the difference between major ticks. Negative 
                values will shift the label downwards. Default is 0.10
            x_TS_label_offset : float, optional
                Horizontal value to offset TS_label from the TS_position. This
                value scales with the x_offset value. Negative values will shift
                the label rightward. Default is 0 (i.e. labels are directly
                above peaks by default)
            TS_label_format : str, optional
                String format to print TS_labels. Uses the `str.format`_ syntax.
                Default is '.2f' (i.e. a floating point value rounded to the 
                second decimal place)
            figure : `matplotlib.figure.Figure`_
                Add plot to this figure. If not specified, one will be generated
            ax : `matplotlib.axes.Axes.axis`_, optional
                Adds plot to this axis. If not specified, one will be generated
            plt_kwargs : dict, optional
                Extra arguments that will be fed to 
                `matplotlib.pyplot.subplots`_
            reaction_kwargs : keyword arguments
                Extra arguments that will be fed to the reactions
        Returns
        -------
            figure : `matplotlib.figure.Figure`_
                Figure
            axes : tuple of `matplotlib.axes.Axes.axis`_
                Axes of the plot.

        .. _`matplotlib.pyplot.subplots`: https://matplotlib.org/api/_as_gen/matplotlib.pyplot.subplots.html
        .. _`matplotlib.figure.Figure`: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html
        .. _`matplotlib.axes.Axes.axis`: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axis.html
        .. _`str.format`: https://docs.python.org/3/library/stdtypes.html#str.format
        """
        # Amount separating species.
        x_offset = 1.
        # Values that will be plot
        x_plot = []
        y_plot = []
        # Labels and locations of the species on the x axis
        x_labels = []
        x_label_pos = []
        # Transition state positions and labels
        if include_TS and include_TS_labels:
            TS_labels = []
            x_TS_label_pos = []
            y_TS_label_pos = []

        x = 0.
        # Calculate reference value
        ref = self.reactions[ref_index].get_state_quantity(
                method_name=method_name, state=ref_state, **reaction_kwargs)

        for i, reaction in enumerate(self.reactions):
            '''First state is the reactants of reaction 0.'''
            if i == 0:
                y_react = reaction.get_state_quantity(method_name=method_name,
                                                      state='reactants',
                                                      **reaction_kwargs) \
                    - ref
                x_plot.extend([x, x+x_offset])
                y_plot.extend([y_react, y_react])
                x_labels.append(_write_reaction_state(
                        species=reaction.reactants,
                        stoich=reaction.reactants_stoich))
                x_label_pos.append(x + x_offset/2.)                
                x += x_offset

            # Calculating y value for product in case we need it to fit the 
            # transition state curve
            y_product = reaction.get_state_quantity(method_name=method_name,
                                                    state='products',
                                                    **reaction_kwargs) \
                        - ref

            '''Calculate properties for TS if necessary'''
            if include_TS and reaction.transition_state is not None:
                x += x_offset
                y_TS = reaction.get_state_quantity(method_name=method_name,
                                                   state='transition state',
                                                   **reaction_kwargs) \
                       - ref
                '''Calculate data to fit spline'''
                x_fit = np.array([x-x_offset,
                                  x-x_offset*(1.-x_scale_TS),
                                  x,
                                  x+x_offset*(1.-x_scale_TS),
                                  x+x_offset])
                y_fit = np.array([y_plot[-1],
                                  (y_TS-y_plot[-1])*y_scale_TS+y_plot[-1],
                                  y_TS,
                                  (y_TS-y_product)*y_scale_TS+y_product,
                                  y_product])
                tck = interpolate.splrep(x_fit, y_fit, k=2)
                '''Calculate new x and y points from spline fit'''
                x_spline = np.linspace(x-x_offset, x+x_offset)
                y_spline = interpolate.splev(x_spline, tck)

                '''Add new data to the appropriate lists'''
                x_plot.extend(x_spline)
                y_plot.extend(y_spline)
                x_labels.append(_write_reaction_state(
                    species=reaction.transition_state, 
                    stoich=reaction.transition_state_stoich))
                x_label_pos.append(x)

                '''Record transition state labels if necessary'''
                if include_TS_labels:
                    TS_labels.append(reaction.get_delta_quantity(
                            method_name=method_name,
                            initial_state='reactants',
                            final_state='transition state',
                            **reaction_kwargs))
                    x_TS_label_pos.append(x + x_TS_label_offset*x_offset)
                    # Add TS y values for now. Correct later we know the graph's
                    # scale
                    y_TS_label_pos.append(y_TS)

            x += x_offset

            '''Calculate properties for product'''
            x_plot.extend([x, x+x_offset])
            y_plot.extend([y_product, y_product])
            x_labels.append(_write_reaction_state(
                    species=reaction.products,
                    stoich=reaction.products_stoich))
            x_label_pos.append(x + x_offset/2.)                
            x += x_offset
        
        ''' Plot the diagram'''
        # Create a new plot if necessary
        if axes is None:
            figure, axes = plt.subplots(1, 1, **plt_kwargs)
        axes.plot(x_plot, y_plot)
        axes.set_xticks(x_label_pos)
        axes.set_xticklabels(x_labels, rotation='vertical')
        axes.set_xlabel('Reaction Coordinate')

        # Setting the y label
        y_label = method_name.replace('get_', '')
        try:
            units = reaction_kwargs['units']
        except KeyError:
            axes.set_ylabel(y_label)
        else:
            axes.set_ylabel('{} ({})'.format(y_label, units))

        # Adding transition state labels
        if include_TS and include_TS_labels:
            y_ticks = axes.get_yticks()
            y_perb = np.ptp(y_ticks)/len(y_ticks)*y_TS_label_offset
            y_TS_label_pos = [y + y_perb for y in y_TS_label_pos]
            label_field = '{:%s}' % TS_label_format
            for TS_label, x_pos, y_pos in zip(TS_labels, x_TS_label_pos, 
                                              y_TS_label_pos):
                axes.text(x=x_pos, y=y_pos, s=label_field.format(TS_label),
                          horizontalalignment='center')

        plt.tight_layout()
        return figure, axes

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = {
            'class': str(self.__class__),
            'reactions': [reaction.to_dict() for reaction in self.reactions],
        }
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
            Reactions : Reactions object
        """
        json_obj = remove_class(json_obj)
        json_obj['reactions'] = [json_to_pMuTT(reaction)
                                 for reaction in json_obj['reactions']]
        return cls(**json_obj)


def _parse_reaction_state(reaction_str, species_delimiter='+'):
    """Takes the reactants/products state of a reaction string and parse it
    into species and stoichiometric amounts

    Parameters
    ----------
        reaction_str : str
            Reactant or product state of reaction
        species_delimiter : str
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
            specie = specie.strip()
            specie_stoich = 1.
            # species.append(specie.strip())
            # stoichiometry.append(1.)
        else:
            # Stoichiometric coefficient present
            specie_stoich = stoich_search.group()
            trim_len = len(specie_stoich)
            specie = specie[trim_len:].strip()
            specie_stoich = float(specie_stoich)
        # If the specie already exists, add its coefficient to existing amount.
        # Otherwise, append the new specie and its stoichiometry
        try:
            i = species.index(specie)
        except ValueError:
            species.append(specie)
            stoichiometry.append(specie_stoich)
        else:
            stoichiometry[i] += specie_stoich
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
            Delimiter that separate states of the reaction. Leading and
            trailing spaces will be trimmed. Default is '='
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
    # Separate states of reaction
    reaction_states = reaction_str.split(reaction_delimiter)

    reactants_state = reaction_states[0]
    products_state = reaction_states[-1]

    # Separate the species in each state
    reactants, reactants_stoich = _parse_reaction_state(
            reaction_str=reactants_state, species_delimiter=species_delimiter)
    products, products_stoich = _parse_reaction_state(
            reaction_str=products_state, species_delimiter=species_delimiter)

    # Check for transition state
    if len(reaction_states) > 2:
        transition_state_state = reaction_states[1]
        transition_state, transition_state_stoich = _parse_reaction_state(
                reaction_str=transition_state_state,
                species_delimiter=species_delimiter)
    else:
        transition_state = None
        transition_state_stoich = None

    return (reactants, reactants_stoich, products, products_stoich,
            transition_state, transition_state_stoich)


def _write_reaction_state(species, stoich, species_delimiter='+',
                          stoich_format='.2f', key='name'):
    """Writes one section of the reaction string

    Parameters
    ----------
        species : list of ``pMuTT`` specie objects
            Species to write
        stoich : list of float
            Stoichiometry corresponding to species
        species_delimiter : str, optional
            Delimiter that separates species. Default is '+'
        stoich_format : float, optional
            Format to write stoichiometric numbers. Default is '.2f' (float
            rounded to 2 decimal places)
        key : str, optional
            Attribute to display species. Default is name
    Returns
    -------
        reaction_str : str
            One section of the reaction string
    """
    # Return blank string if no species present for the reaction state
    if species is None:
        return ''

    stoich_field = '{:%s}' % stoich_format
    for i, (specie, stoich_val) in enumerate(zip(species, stoich)):
        # If the coefficient is 1, just write the specie name
        if np.isclose(stoich_val, 1.):
            specie_str = getattr(specie, key)
        else:
            stoich_val = stoich_field.format(stoich_val)
            specie_str = '{}{}'.format(stoich_val, getattr(specie, key))

        # Specie delimiter only written after the first specie
        if i == 0:
            reaction_str = specie_str
        else:
            reaction_str += '{}{}'.format(species_delimiter, specie_str)
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


def _get_molecularity(stoich):
    """Calculate the molecularity (e.g. unimolecular, biomolecular)

    Parameters
    ----------
        stoich : list of float
            Stoichiometric coefficients of each specie
    Returns
    -------
        m : int
            Molecularity of reaction
    """
    return np.sum(stoich)

def _get_states(rev, act):
    """Determines the initial state and the final state based on boolean inputs

    Parameters
    ----------
        rev : bool
            True if the reaction is in the reverse direction.
        act : bool
            True if the transition state is the final state.
    Returns
    -------
        initial_state : str
            Initial state of the reaction. Either 'reactants', 'products' or
            'transition state'
        final_state : str
            Final state of the reaction. Either 'reactants', 'products' or
            'transition state'
    """
    if rev:
        initial_state = 'products'
        final_state = 'reactants'
    else:
        initial_state = 'reactants'
        final_state = 'products'

    # Overwrites the final state if necessary
    if act:
        final_state = 'transition state'

    return initial_state, final_state