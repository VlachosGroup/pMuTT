# -*- coding: utf-8 -*-
from pMuTT import constants as c
from pMuTT.io.json import remove_class


class BEP:
    """Represents a Bronsted Evans Polyani relationship

    :math:`E_a = \\alpha H + \\beta`

    Attributes
    ----------
        slope : float
            Slope of BEP relationship.
        intercept : float
            Intercept of BEP relationship in kcal/mol.
        reaction : :class:`~pMuTT.reaction.Reaction` object, optional
            Reaction related to BEP. The Reaction does not need to be supplied
            immediately
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
        _descriptor : method
            Method taken from reaction to calculate enthalpy. This attribute is
            not supplied to constructor.
    """

    def __init__(self, slope, intercept, reaction=None, descriptor='delta_H'):
        self.slope = slope
        self.intercept = intercept
        self.set_descriptor(reaction=reaction, descriptor=descriptor)

    def set_descriptor(self, reaction=None, descriptor=None):
        """Sets the appropriate method handle to the BEP object

        Parameters
        ----------
            reaction : :class:`~pMuTT.reaction.Reaction` object, optional
                Reaction related to BEP. If specified, overwrites the value
                held by BEP object
            descriptor : str, optional
                Descriptor to calculate the activation energy. Supported
                options:

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
                
                If specified, overwites the value held by BEP object.
        """
        try:
            self.reaction
        except AttributeError:
            # Assigns reaction to default value
            self.reaction = reaction
        else:
            # Overwrites reaction if not previously set
            if reaction is not None:
                self.reaction = reaction

        try:
            self.descriptor
        except AttributeError:
            # Assigns descriptor to default value
            self.descriptor = descriptor
        else:
            # Overwrites reaction if not previously set
            if descriptor is not None:
                self.descriptor = descriptor

        if self.descriptor is None or self.reaction is None:
            self._descriptor = None
        elif self.descriptor == 'delta_H':
            self._descriptor = \
                    lambda **kwargs: self.reaction.get_delta_H(rev=False,
                                                               units='kcal/mol',
                                                               **kwargs)
        elif self.descriptor == 'rev_delta_H':
            self._descriptor = \
                    lambda **kwargs: self.reaction.get_delta_H(rev=True,
                                                               units='kcal/mol',
                                                               **kwargs)
        elif self.descriptor == 'reactants_H':
            self._descriptor = \
                    lambda **kwargs: self.reaction.get_H_state(
                            units='kcal/mol',
                            state='reactants',
                            **kwargs)
        elif descriptor == 'products_H':
            self._descriptor = \
                    lambda **kwargs: self.reaction.get_H_state(
                            units='kcal/mol',
                            state='products',
                            **kwargs)
        elif self.descriptor == 'delta_E':
            self._descriptor = \
                    lambda **kwargs: self.reaction.get_delta_E(rev=False,
                                                               units='kcal/mol',
                                                               **kwargs)
        elif self.descriptor == 'rev_delta_E':
            self._descriptor = \
                    lambda **kwargs: self.reaction.get_delta_E(rev=True,
                                                               units='kcal/mol',
                                                               **kwargs)
        elif self.descriptor == 'reactants_E':
            self._descriptor = \
                    lambda **kwargs: self.reaction.get_E_state(
                            units='kcal/mol',
                            state='reactants',
                            **kwargs)
        elif self.descriptor == 'products_E':
            self._descriptor = \
                    lambda **kwargs: self.reaction.get_E_state(
                            units='kcal/mol',
                            state='products',
                            **kwargs)
        else:
            raise ValueError(('Descriptor "{}" not supported. See '
                              'documentation of pMuTT.reaction.bep.BEP for '
                              'supported options.'.format(self.descriptor)))

    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict

    def get_E_act(self, units, rev=False, **kwargs):
        """Calculate Arrhenius activation energy using BEP relationship

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
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
        if 'rev_delta' in self.descriptor:
            # If the descriptor is for the reverse reaction, the slope has to
            # be modified
            if rev:
                E_act = self.slope*self._descriptor(**kwargs) + self.intercept
            else:
                E_act = (self.slope-1.)*self._descriptor(**kwargs) \
                        + self.intercept
        else:
            if rev:
                E_act = (self.slope-1.)*self._descriptor(**kwargs) \
                        + self.intercept
            else:
                E_act = self.slope*self._descriptor(**kwargs) + self.intercept
        return E_act*c.convert_unit(initial='kcal/mol', final=units)

    def get_EoRT_act(self, rev=False, T=c.T0('K'), **kwargs):
        """Calculate dimensionless Arrhenius activation energy using BEP
        relationship

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate the descriptor
        Returns
        -------
            EoRT_act : float
                Dimensionless activation energy
        """
        return self.get_E_act(units='kcal/mol', rev=rev, T=T, **kwargs) \
               /c.R('kcal/mol/K')/T


    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = {
            'class': str(self.__class__),
            'slope': self.slope,
            'intercept': self.intercept,
            'descriptor': self.descriptor}
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
        return cls(**json_obj)
