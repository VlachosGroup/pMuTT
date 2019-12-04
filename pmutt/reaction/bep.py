# -*- coding: utf-8 -*-
from pmutt import _ModelBase
from pmutt import constants as c
from pmutt.io.json import remove_class


class BEP(_ModelBase):
    """Represents a Bronsted Evans Polyani relationship. Intended to represent
    a species. Inherits from :class:`~pmutt._ModelBase`.

    :math:`E_a = \\alpha H + \\beta`

    Attributes
    ----------
        slope : float
            Slope of BEP relationship.
        intercept : float
            Intercept of BEP relationship in kcal/mol.
        name : str, optional
            Name of the BEP. Default is None
        descriptor : str, optional
            Descriptor to calculate the activation energy. Supported options:

            ===========  ===========================================================================
            Descriptor   Description
            ===========  ===========================================================================
            delta_H      H_products - H_reactants (change in enthalpy)
            rev_delta_H  H_reactants - H_products (change in enthalpy in reverse direction)
            reactants_H  H_reactants (enthalpy of reactants)
            products_H   H_products (enthalpy of products)
            delta_E      E_products - E_reactants (change in electronic energy)
            rev_delta_E  E_reactants - E_products (change in electronic energy in reverse direction)
            reactants_E  E_reactants (electronic energy of reactants)
            products_E   E_products (electronic energy of products)
            ===========  ===========================================================================

            Default is 'delta_H'.
        elements : dict
            Composition of the species.
            Keys of dictionary are elements, values are stoichiometric
            values in a formula unit.
            e.g. CH3OH can be represented as:
            {'C': 1, 'H': 4, 'O': 1,}.
    
        notes : str or dict
            Notes relevant to BEP relationship such as its source. If using a
            dictionary, the keys and values must be simple types supported by
            JSON
    """

    def __init__(self, slope, intercept, name=None, reaction=None,
                 descriptor='delta_H', elements=None, notes=None):
        self.name = name
        self.slope = slope
        self.intercept = intercept
        self.descriptor = descriptor
        self.elements = elements
        self.notes = notes

    def _get_descriptor_val(self, reaction, **kwargs):
        """Sets the appropriate method handle to the BEP object

        Parameters
        ----------
            reaction : :class:`~pmutt.reaction.Reaction` object
                Reaction related to BEP.
            kwargs : keyword arguments
                Parameters required to calculate descriptor value such as
                temperature and pressure.
        Returns
        -------
            val : float
                Descriptor value in kcal/mol. Note that kcal/mol is the standard
                unit since the intercept is also in kcal/mol.
        """
        # Assign the descriptor to the reaction
        if self.descriptor == 'delta_H':
            val = reaction.get_delta_H(rev=False, units='kcal/mol', **kwargs)
        elif self.descriptor == 'rev_delta_H':
            val = reaction.get_delta_H(rev=True, units='kcal/mol', **kwargs)
        elif self.descriptor == 'reactants_H':
            val = reaction.get_H_state(units='kcal/mol', state='reactants',
                                       **kwargs)
        elif self.descriptor == 'products_H':
            val = reaction.get_H_state(units='kcal/mol', state='products',
                                       **kwargs)
        elif self.descriptor == 'delta_E':
            val = reaction.get_delta_E(rev=False, units='kcal/mol', **kwargs)
        elif self.descriptor == 'rev_delta_E':
            val = reaction.get_delta_E(rev=True, units='kcal/mol', **kwargs)
        elif self.descriptor == 'reactants_E':
            val = reaction.get_E_state(units='kcal/mol', state='reactants',
                                       **kwargs)
        elif self.descriptor == 'products_E':
            val = reaction.get_E_state(units='kcal/mol', state='products',
                                       **kwargs)
        else:
            err_msg = ('Descriptor "{}" not supported. See documentation of '
                       'pmutt.reaction.bep.BEP for supported options.'
                       ''.format(self.descriptor))
            raise ValueError(err_msg)
        return val

    def _get_adjusted_slope(self, rev):
        """Calculates the adjusted slope based on the descriptor and the
        direction of the reaction

        Parameters
        ----------
            rev : bool
                If True, the reaction is reversible.
        Returns
        -------
            adj_slope : float
                Adjusted slope
        """
        # If the descriptor is for the reverse reaction, the slope has to
        # be modified
        if 'rev_delta' in self.descriptor:
            if rev:
                adj_slope = self.slope
            else:
                adj_slope = self.slope - 1.
        else:
            if rev:
                adj_slope = self.slope - 1.
            else:
                adj_slope = self.slope
        return adj_slope


    def get_E_act(self, units, reaction, rev=False, **kwargs):
        """Calculate Arrhenius activation energy using BEP relationship

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            reaction : :class:`~pmutt.reaction.Reaction` object
                Reaction related to BEP.
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate the descriptor
        Returns
        -------
            E_act : float
                Dimensionless activation energy
        """
        adj_slope = self._get_adjusted_slope(rev=rev)
        descriptor_val = self._get_descriptor_val(reaction=reaction, **kwargs)
        E_act = adj_slope*descriptor_val + self.intercept
        return E_act*c.convert_unit(initial='kcal/mol', final=units)

    def get_EoRT_act(self, reaction, rev=False, T=c.T0('K'), **kwargs):
        """Calculates dimensionless Arrhenius activation energy using BEP
        relationship

        Parameters
        ----------
            reaction : :class:`~pmutt.reaction.Reaction` object
                Reaction related to BEP.
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            T : float, optional
                Temperature in K. Default is 298.15
            kwargs : keyword arguments
                Parameters required to calculate the descriptor
        Returns
        -------
            EoRT_act : float
                Dimensionless activation energy
        """
        return self.get_E_act(units='kcal/mol', reaction=reaction,
                              rev=rev, T=T, **kwargs)/c.R('kcal/mol/K')/T

    def get_UoRT(self, reaction, T=c.T0('K'), **kwargs):
        """Calculates the dimensionless internal energy using BEP relationship
        and initial state internal energy
        
        Parameters
        ----------
            reaction : :class:`~pmutt.reaction.Reaction` object
                Reaction related to BEP.
            T : float, optional
                Temperature in K. Default is 298.15
            kwargs : keyword arguments
                Parameters required to calculate the descriptor
        Returns
        -------
            UoRT : float
                Dimensionless internal energy
        """
        UoRT_reactants = reaction.get_UoRT_state(state='reactants', T=T,
                                                 **kwargs)
        return self.get_EoRT_act(rev=True, reaction=reaction, T=T,
                                 **kwargs) + UoRT_reactants

    def get_HoRT(self, reaction, T=c.T0('K'), **kwargs):
        """Calculates the dimensionless enthalpy using BEP relationship
        and reactants or products enthalpy
        
        Parameters
        ----------
            reaction : :class:`~pmutt.reaction.Reaction` object
                Reaction related to BEP.
            T : float, optional
                Temperature in K. Default is 298.15
            kwargs : keyword arguments
                Parameters required to calculate the descriptor
        Returns
        -------
            HoRT : float
                Dimensionless enthalpy
        """
        HoRT_reactants = reaction.get_HoRT_state(state='reactants', T=T,
                                                 **kwargs)
        return self.get_EoRT_act(rev=False, reaction=reaction, T=T, **kwargs) \
               + HoRT_reactants

    def get_SoR(self, reaction=None, T=c.T0('K'), entropy_state='reactants',
                **kwargs):
        """Calculates the dimensionless entropy using reactants or products
        entropy. The BEP relationship has no entropic contribution
        
        Parameters
        ----------
            reaction : :class:`~pmutt.reaction.Reaction` object, optional
                Reaction related to BEP. If `entropy_state` is None, `reaction`
                is not required.
            T : float, optional
                Temperature in K. Default is 298.15
            entropy_state : str or None, optional
                State to use to estimate entropy. Supported arguments:

                - 'reactants' (default)
                - 'products'
                - None (Entropy contribution is 0. Useful if misc_models
                  have been specified for entropy)
            kwargs : keyword arguments
                Parameters required to calculate the descriptor
        Returns
        -------
            SoR : float
                Dimensionless entropy
        """
        if entropy_state is None:
            SoR = 0.
        else:
            try:
                SoR = reaction.get_SoR_state(state=entropy_state, T=T, **kwargs)
            except AttributeError:
                err_msg = ('Unable to calculate SoR of BEP object since '
                           'a Reaction object was not supplied. Either set '
                           'entropy_state to None or supply a Reaction to '
                           'get_SoR.')
                raise ValueError(err_msg)
        return SoR

    def get_GoRT(self, reaction, T=c.T0('K'), entropy_state='reactants', 
                 **kwargs):
        """Calculates the dimensionless Gibbs energy using BEP relationship
        and reactants Gibbs energy. The BEP relationship has no entropic
        contribution
        
        Parameters
        ----------
            reaction : :class:`~pmutt.reaction.Reaction` object
                Reaction related to BEP.
            T : float, optional
                Temperature in K. Default is 298.15
            entropy_state : str or None, optional
                State to use to estimate entropy. Supported arguments:

                - 'reactants' (default)
                - 'products'
                - None (Entropy contribution is 0. Useful if misc_models
                  have been specified for entropy)
            kwargs : keyword arguments
                Parameters required to calculate the descriptor
        Returns
        -------
            GoRT : float
                Dimensionless Gibbs energy
        """
        HoRT = self.get_HoRT(T=T, reaction=reaction, **kwargs)
        SoR = self.get_SoR(T=T, reaction=reaction, entropy_state=entropy_state,
                           **kwargs)
        return HoRT - SoR

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = {
            'class': str(self.__class__),
            'name': self.name,
            'slope': self.slope,
            'intercept': self.intercept,
            'descriptor': self.descriptor,
            'notes': self.notes}
        return obj_dict