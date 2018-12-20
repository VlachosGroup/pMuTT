# -*- coding: utf-8 -*-
from collections import Counter
import inspect
import re
import numpy as np
from pMuTT import _force_pass_arguments, _pass_expected_arguments, _is_iterable
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

        if None not in self.transition_state:
            TS_elements = _count_elements(self.transition_state,
                                          self.transition_state_stoich)
            if reactant_elements != TS_elements:
                raise ValueError('Number of elements in reactants and '
                                 'transition state do not agree.\n'
                                 'Reactant count: {}\n'
                                 'Product count: {}'.format(reactant_elements,
                                                            TS_elements))

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
        return self._get_state_quantity(state=state, method_name='get_q',
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
        return self._get_state_quantity(state=state, method_name='get_CvoR',
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
        return self._get_state_quantity(state=state, method_name='get_CpoR',
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
        return self._get_state_quantity(state=state, method_name='get_UoRT',
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
        return self._get_state_quantity(state=state, method_name='get_HoRT',
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
        return self._get_state_quantity(state=state, method_name='get_SoR',
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
        return self._get_state_quantity(state=state, method_name='get_FoRT',
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
        return self._get_state_quantity(state=state, method_name='get_GoRT',
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

    def get_delta_q(self, rev=False, **kwargs):
        """Gets change in partition function between reactants and products

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
            delta_q : float
                Change in partition function between reactants and products
        """
        if rev:
            delta_q = self._get_delta_quantity(initial_state='products',
                                               final_state='reactants',
                                               method_name='get_q', **kwargs)
        else:
            delta_q = self._get_delta_quantity(initial_state='reactants',
                                               final_state='products',
                                               method_name='get_q', **kwargs)
        return delta_q

    def get_delta_CvoR(self, rev=False, **kwargs):
        """Gets change in dimensionless heat capacity (constant V)
        between reactants and products

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
            delta_CvoR : float
                Change in heat capacity between reactants and products
        """
        if rev:
            delta_CvoR = self._get_delta_quantity(initial_state='products',
                                                  final_state='reactants',
                                                  method_name='get_CvoR',
                                                  **kwargs)
        else:
            delta_CvoR = self._get_delta_quantity(initial_state='reactants',
                                                  final_state='products',
                                                  method_name='get_CvoR',
                                                  **kwargs)
        return delta_CvoR

    def get_delta_Cv(self, units, rev=False, **kwargs):
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
            kwargs : keyword arguments
                Parameters required to calculate heat capacity. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_Cv : float
                Change in heat capacity between reactants and products
        """
        return self.get_delta_CvoR(rev=rev, **kwargs)*c.R(units)

    def get_delta_CpoR(self, rev=False, **kwargs):
        """Gets change in dimensionless heat capacity between reactants and
        products

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
            delta_CvoR : float
                Change in heat capacity between reactants and products
        """
        if rev:
            delta_CpoR = self._get_delta_quantity(initial_state='products',
                                                  final_state='reactants',
                                                  method_name='get_CpoR',
                                                  **kwargs)
        else:
            delta_CpoR = self._get_delta_quantity(initial_state='reactants',
                                                  final_state='products',
                                                  method_name='get_CpoR',
                                                  **kwargs)
        return delta_CpoR

    def get_delta_Cp(self, units, rev=False, **kwargs):
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
            kwargs : keyword arguments
                Parameters required to calculate heat capacity. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_Cp : float
                Change in heat capacity between reactants and products
        """
        return self.get_delta_CpoR(rev=rev, **kwargs)*c.R(units)

    def get_delta_UoRT(self, rev=False, **kwargs):
        """Gets change in dimensionless internal energy between reactants and
        products

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
            delta_UoRT : float
                Change in internal energy between reactants and products
        """
        if rev:
            delta_UoRT = self._get_delta_quantity(initial_state='products',
                                                  final_state='reactants',
                                                  method_name='get_UoRT',
                                                  **kwargs)
        else:
            delta_UoRT = self._get_delta_quantity(initial_state='reactants',
                                                  final_state='products',
                                                  method_name='get_UoRT',
                                                  **kwargs)
        return delta_UoRT

    def get_delta_U(self, units, T, rev=False, **kwargs):
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
            kwargs : keyword arguments
                Parameters required to calculate internal energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_U : float
                Change in internal energy between reactants and products
        """
        return self.get_delta_UoRT(rev=rev, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

    def get_delta_HoRT(self, rev=False, **kwargs):
        """Gets change in dimensionless enthalpy between reactants and products

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
            delta_HoRT : float
                Change in enthalpy between reactants and products
        """
        if rev:
            delta_HoRT = self._get_delta_quantity(initial_state='products',
                                                  final_state='reactants',
                                                  method_name='get_HoRT',
                                                  **kwargs)
        else:
            delta_HoRT = self._get_delta_quantity(initial_state='reactants',
                                                  final_state='products',
                                                  method_name='get_HoRT',
                                                  **kwargs)
        return delta_HoRT

    def get_delta_H(self, units, T, rev=False, **kwargs):
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
            kwargs : keyword arguments
                Parameters required to calculate enthalpy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_H : float
                Change in enthalpy between reactants and products
        """
        return self.get_delta_HoRT(rev=rev, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

    def get_delta_SoR(self, rev=False, **kwargs):
        """Gets change in dimensionless entropy between reactants and products

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
            delta_SoR : float
                Change in entropy between reactants and products
        """
        if rev:
            delta_SoR = self._get_delta_quantity(initial_state='products',
                                                 final_state='reactants',
                                                 method_name='get_SoR',
                                                 **kwargs)
        else:
            delta_SoR = self._get_delta_quantity(initial_state='reactants',
                                                 final_state='products',
                                                 method_name='get_SoR',
                                                 **kwargs)
        return delta_SoR

    def get_delta_S(self, units, rev=False, **kwargs):
        """Gets change in entropy between reactants and products

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
            delta_S : float
                Change in entropy between reactants and products
        """
        return self.get_delta_SoR(rev=rev, **kwargs)*c.R(units)

    def get_delta_FoRT(self, rev=False, **kwargs):
        """Gets change in dimensionless Helmholtz energy between reactants and
        products

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
            delta_FoRT : float
                Change in Helmholtz energy between reactants and products
        """
        if rev:
            delta_FoRT = self._get_delta_quantity(initial_state='products',
                                                  final_state='reactants',
                                                  method_name='get_FoRT',
                                                  **kwargs)
        else:
            delta_FoRT = self._get_delta_quantity(initial_state='reactants',
                                                  final_state='products',
                                                  method_name='get_FoRT',
                                                  **kwargs)
        return delta_FoRT

    def get_delta_F(self, units, T, rev=False, **kwargs):
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
            kwargs : keyword arguments
                Parameters required to calculate Helmholtz energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_F : float
                Change in Helmholtz energy between reactants and products
        """
        return self.get_delta_FoRT(rev=rev, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

    def get_delta_GoRT(self, rev=False, **kwargs):
        """Gets change in dimensionless Gibbs energy between reactants and
        products

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
            delta_GoRT : float
                Change in Gibbs energy between reactants and products
        """
        if rev:
            delta_GoRT = self._get_delta_quantity(initial_state='products',
                                                  final_state='reactants',
                                                  method_name='get_GoRT',
                                                  **kwargs)
        else:
            delta_GoRT = self._get_delta_quantity(initial_state='reactants',
                                                  final_state='products',
                                                  method_name='get_GoRT',
                                                  **kwargs)
        return delta_GoRT

    def get_delta_G(self, units, T, rev=False, **kwargs):
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
            kwargs : keyword arguments
                Parameters required to calculate Gibbs energy. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            delta_G : float
                Change in Gibbs energy between reactants and products
        """
        return self.get_delta_GoRT(rev=rev, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

    def get_Keq(self, rev=False, **kwargs):
        """Gets equilibrium constant between reactants and products

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate equilibrium constant. See
                class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            Keq : float
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
                Parameters required to calculate partition function. See class
                docstring to see how to pass specific parameters to different
                species.
        Returns
        -------
            q_act : float
                Change in partition function between reactants (or products)
                and transition state
        """
        if rev:
            q_act = self._get_delta_quantity(initial_state='products',
                                             final_state='transition state',
                                             method_name='get_q',
                                             **kwargs)
        else:
            q_act = self._get_delta_quantity(initial_state='reactants',
                                             final_state='transition state',
                                             method_name='get_q',
                                             **kwargs)
        return q_act

    def get_CvoR_act(self, rev=False, **kwargs):
        """Gets change in dimensionless heat capacity between reactants (or
        products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless heat capacity.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            CvoR_act : float
                Change in dimensionless heat capacity between reactants (or
                products) and transition state
        """
        if rev:
            CvoR_act = self._get_delta_quantity(initial_state='products',
                                                final_state='transition state',
                                                method_name='get_CvoR',
                                                **kwargs)
        else:
            CvoR_act = self._get_delta_quantity(initial_state='reactants',
                                                final_state='transition state',
                                                method_name='get_CvoR',
                                                **kwargs)
        return CvoR_act

    def get_Cv_act(self, units, rev=False, **kwargs):
        """Gets change in heat capacity between reactants (or
        products) and transition state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate heat capacity.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            Cv_act : float
                Change in heat capacity between reactants (or products) and
                transition state
        """
        return self.get_CvoR_act(rev=rev, **kwargs)*c.R(units)

    def get_CpoR_act(self, rev=False, **kwargs):
        """Gets change in dimensionless heat capacity between reactants (or
        products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless heat capacity.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            CpoR_act : float
                Change in dimensionless heat capacity between reactants (or
                products) and transition state
        """
        if rev:
            CpoR_act = self._get_delta_quantity(initial_state='products',
                                                final_state='transition state',
                                                method_name='get_CpoR',
                                                **kwargs)
        else:
            CpoR_act = self._get_delta_quantity(initial_state='reactants',
                                                final_state='transition state',
                                                method_name='get_CpoR',
                                                **kwargs)
        return CpoR_act

    def get_Cp_act(self, units, rev=False, **kwargs):
        """Gets change in heat capacity between reactants (or products) and
        transition state

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
                docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            Cp_act : float
                Change in heat capacity between reactants (or products) and
                transition state
        """
        return self.get_CpoR_act(rev=rev, **kwargs)*c.R(units)

    def get_UoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless internal energy between reactants (or
        products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless internal energy.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            UoRT_act : float
                Change in dimensionless internal energy between reactants
                (or products) and transition state
        """
        if rev:
            UoRT_act = self._get_delta_quantity(initial_state='products',
                                                final_state='transition state',
                                                method_name='get_UoRT',
                                                **kwargs)
        else:
            UoRT_act = self._get_delta_quantity(initial_state='reactants',
                                                final_state='transition state',
                                                method_name='get_UoRT',
                                                **kwargs)
        return UoRT_act

    def get_U_act(self, units, T, rev=False, **kwargs):
        """Gets change in internal energy between reactants (or products) and
        transition state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate internal energy. See class
                docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            U_act : float
                Change in internal energy between reactants (or products) and
                transition state
        """
        return self.get_UoRT_act(rev=rev, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

    def get_HoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless enthalpy between reactants
        (or products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless enthalpy. See
                class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            HoRT_act : float
                Change in dimensionless enthalpy between reactants
                (or products) and transition state
        """
        if rev:
            HoRT_act = self._get_delta_quantity(initial_state='products',
                                                final_state='transition state',
                                                method_name='get_HoRT',
                                                **kwargs)
        else:
            HoRT_act = self._get_delta_quantity(initial_state='reactants',
                                                final_state='transition state',
                                                method_name='get_HoRT',
                                                **kwargs)
        return HoRT_act

    def get_H_act(self, units, T, rev=False, **kwargs):
        """Gets change in enthalpy between reactants (or products) and
        transition state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate enthalpy. See class
                docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            H_act : float
                Change in enthalpy between reactants (or products) and
                transition state
        """
        return self.get_HoRT_act(rev=rev, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

    def get_EoRT_act(self, rev=False, method='any', **kwargs):
        """Gets dimensionless activation energy between reactants
        (or products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            method : str, optional
                Method to use to calculate dimensionless activation energy.
                Accepted options:

                - 'any' (uses whichever is available. ``self.transition_state``
                  is used preferentially over ``self.BEP``)
                - 'bep' (uses ``self.bep``)
                - 'ts' or 'transition_state' (uses ``self.transition_state``)

                Default is 'any'.
            kwargs : keyword arguments
                Parameters required to calculate dimensionless activation
                energy
        Returns
        -------
            EoRT_act : float
                Dimensionless activation energy between reactants (or products)
                and transition state
        """
        method = method.lower()
        if method == 'any':
            if None not in self.transition_state:
                return self.get_HoRT_act(rev=rev, **kwargs)
            else:
                return self.bep.get_EoRT_act(rev=rev, **kwargs)
        elif method == 'bep':
            return self.bep.get_EoRT_act(rev=rev, **kwargs)
        elif method == 'transition state' or method == 'ts':
            return self.get_HoRT_act(rev=rev, **kwargs)
        else:
            raise ValueError(('Method "{}" not supported. See documentation '
                              'of '
                              '``pMuTT.reaction.Reaction.get_EoRT_act`` '
                              'for supported options.'.format(method)))

    def get_E_act(self, units, T, rev=False, method='any', **kwargs):
        """Gets change in activation energy between reactants (or products)
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
                Method to use to calculate activation energy.
                Accepted options:

                - 'any' (uses whichever is available. ``self.transition_state``
                  is used preferentially over ``self.BEP``)
                - 'bep' (uses ``self.bep``)
                - 'ts' or 'transition_state' (uses ``self.transition_state``)

                Default is 'any'.
            kwargs : keyword arguments
                Parameters required to calculate activation energy. See class
                docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            E_act : float
                Change in activation energy between reactants (or products) and
                transition state
        """
        return self.get_EoRT_act(rev=rev, method=method, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

    def get_SoR_act(self, rev=False, **kwargs):
        """Gets change in dimensionless entropy between reactants (or products)
        and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless entropy. See
                class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            SoR_act : float
                Change in dimensionless entropy between reactants (or products)
                and transition state
        """
        if rev:
            SoR_act = self._get_delta_quantity(initial_state='products',
                                               final_state='transition state',
                                               method_name='get_SoR',
                                               **kwargs)
        else:
            SoR_act = self._get_delta_quantity(initial_state='reactants',
                                               final_state='transition state',
                                               method_name='get_SoR',
                                               **kwargs)
        return SoR_act

    def get_S_act(self, units, rev=False, **kwargs):
        """Gets change in entropy between reactants (or products) and
        transition state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate entropy. See class
                docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            S_act : float
                Change in entropy between reactants (or products) and
                transition state
        """
        return self.get_SoR_act(rev=rev, **kwargs)*c.R(units)

    def get_FoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless Helmholtz energy between reactants (or
        products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless Helmholtz
                energy. See class docstring to see how to pass specific
                parameters to different species.
        Returns
        -------
            FoRT_act : float
                Change in dimensionless Helmholtz energy between reactants
                (or products) and transition state
        """
        if rev:
            FoRT_act = self._get_delta_quantity(initial_state='products',
                                                final_state='transition state',
                                                method_name='get_FoRT',
                                                **kwargs)
        else:
            FoRT_act = self._get_delta_quantity(initial_state='reactants',
                                                final_state='transition state',
                                                method_name='get_FoRT',
                                                **kwargs)
        return FoRT_act

    def get_F_act(self, units, T, rev=False, **kwargs):
        """Gets change in Helmholtz energy between reactants (or products) and
        transition state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate Helmholtz energy. See class
                docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            F_act : float
                Change in Helmholtz energy between reactants (or products) and
                transition state
        """
        return self.get_FoRT_act(rev=rev, T=T, **kwargs)*T \
            * c.R('{}/K'.format(units))

    def get_GoRT_act(self, rev=False, **kwargs):
        """Gets change in dimensionless Gibbs energy between reactants (or
        products) and transition state

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate dimensionless Gibbs energy.
                See class docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            GoRT_act : float
                Change in dimensionless Gibbs energy between reactants (or
                products) and transition state
        """
        if rev:
            GoRT_act = self._get_delta_quantity(initial_state='products',
                                                final_state='transition state',
                                                method_name='get_GoRT',
                                                **kwargs)
        else:
            GoRT_act = self._get_delta_quantity(initial_state='reactants',
                                                final_state='transition state',
                                                method_name='get_GoRT',
                                                **kwargs)
        return GoRT_act

    def get_G_act(self, units, T, rev=False, **kwargs):
        """Gets change in Gibbs energy between reactants (or products) and
        transition state

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate Gibbs energy. See class
                docstring to see how to pass specific parameters to
                different species.
        Returns
        -------
            G_act : float
                Change in Gibbs energy between reactants (or products) and
                transition state
        """
        return self.get_GoRT_act(rev=rev, T=T, **kwargs)*T \
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

        return c.kb('J/K')*T/c.h('J s')\
            * np.exp(self.get_SoR_act(rev=rev, T=T, **kwargs)+m)

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

    def _get_state_quantity(self, state, method_name, **kwargs):
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
                quantity as long as the relevant objects' have the same method
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
            method = getattr(specie, method_name)
            if method_name == 'get_q':
                state_quantity *= \
                        _force_pass_arguments(method, **specie_kwargs)**stoich
            else:
                state_quantity += \
                        _force_pass_arguments(method, **specie_kwargs)*stoich
        return state_quantity

    def _get_delta_quantity(self, initial_state, final_state, method_name,
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
                quantity as long as the relevant objects' have the same method
                name
        Returns
        -------
            delta_quantity : float
                Change in thermodynamic quantity between particular states
        """
        initial_quantity = self._get_state_quantity(state=initial_state,
                                                    method_name=method_name,
                                                    **kwargs)
        final_quantity = self._get_state_quantity(state=final_state,
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

    def to_string(self, species_delimiter='+', reaction_delimiter='='):
        """Writes the Reaction object as a stoichiometric reaction

        Parameters
        ----------
            species_delimiter : str, optional
                Separates species. Default is '+'
            reaction_delimiter : str, optional
                Separates reaction states. Default is '='
        Returns
        -------
            reaction_str : str
                Reaction string
        """
        # Write reactants
        reaction_str = _write_reaction_state(species=self.reactants,
                                             stoich=self.reactants_stoich,
                                             species_delimiter=species_delimiter)
        reaction_str += reaction_delimiter

        # Write transition state if any
        reaction_str += _write_reaction_state(
                species=self.transition_state,
                stoich=self.transition_state_stoich,
                species_delimiter=species_delimiter)
        reaction_str += reaction_delimiter

        # Write products
        reaction_str += _write_reaction_state(species=self.products,
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
            'products_stoich': list(self.products_stoich),
            }
        try:
            obj_dict['transition_state'] = [ts.to_dict()
                                            for ts in self.transition_state]
        except (AttributeError, TypeError):
            obj_dict['transition_state'] = self.transition_state
            obj_dict['transition_state_stoich'] = self.transition_state_stoich
        else:
            obj_dict['transition_state_stoich'] = \
                    list(self.transition_state_stoich)
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
        stoich_search = re.search('^\d+\.?\d*', specie)
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


def _write_reaction_state(species, stoich, species_delimiter='+'):
    """Writes one section of the reaction string

    Parameters
    ----------
        species : list of ``pMuTT`` specie objects
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
        if specie is None:
            reaction_str = ''
            break
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


def _get_specie_kwargs(specie_name, **kwargs):
    """Gets the keyword arguments specific to a specie

    Parameters
    ----------
        specie_name : str
            Name of the specie
        kwargs : keyword arguments
            Parameters with the conditions. Specie specific parameters can be
            passed by having a key named 'specie' mapping onto a dictionary
            whose keys are the species names.

            e.g. For the reaction: H2 + 0.5O2 = H2O
            kwargs = {
                'T': 298.,
                'specie': {
                    'H2': {
                        'P': 1.,
                    },
                    'O2': {
                        'P': 0.5,
                    },
                }
            }
    Returns
    -------
        specie_kwargs : dict
            Dictionary containing the specie-specific kwargs
    """
    specie_kwargs = kwargs.copy()
    specie_specific = specie_kwargs.pop('specie', None)
    # See if there was an entry for the specific species
    try:
        specie_kwargs.update(specie_specific[name])
    except (KeyError, TypeError, NameError):
        pass
    return specie_kwargs
