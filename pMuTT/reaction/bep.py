# -*- coding: utf-8 -*-
from pMuTT.io_.jsonio import remove_class

class BEP:
    """Represents a Bronsted Evans Polyani relationship

    Attributes
    ----------
        slope : float
            Slope of BEP relationship. 
        intercept : float
            Intercept of BEP relationship in dimensionless units.
        reaction : :class:`~pMuTT.reaction.Reaction` object, optional
            Reaction related to BEP. The Reaction does not need to be supplied 
            immediately
        descriptor : str, optional
            Descriptor to calculate the activation energy. Supported options:
            - 'delta' (H_products - H_reactants)
            - 'rev_delta' (H_reactants - H_products)
            - 'reactants' (H_reactants)
            - 'products' (H_products)
            Default is delta.
        _descriptor : method
            Method taken from reaction to calculate enthalpy. This attribute is 
            not supplied to constructor.
    """

    def __init__(self, slope, intercept, reaction=None, descriptor='delta'):
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
                - 'delta' (H_products - H_reactants)
                - 'rev_delta' (H_reactants - H_products)
                - 'reactants' (H_reactants)
                - 'products' (H_products)
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
        elif self.descriptor == 'delta':
            self._descriptor = \
                    lambda **kwargs: self.reaction.get_delta_HoRT(rev=False, 
                                                                  **kwargs)
        elif self.descriptor == 'rev_delta':
            self._descriptor = \
                    lambda **kwargs: self.reaction.get_delta_HoRT(rev=True, 
                                                                  **kwargs)
        elif self.descriptor == 'reactants' or descriptor == 'products':
            self._descriptor = \
                    lambda **kwargs: self.reaction.get_HoRT_state(
                            state=descriptor, 
                            **kwargs)
        else:
            raise ValueError(('Descriptor "{}" not supported. See documentation'
                              ' of pMuTT.reaction.bep.BEP for supported '
                              'options.'.format(self.descriptor)))

    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict

    def get_EoRT_act(self, rev=False, **kwargs):
        """Calculate dimensionless activation energy using BEP relationship

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
        if self.descriptor == 'rev_delta':
            # If the descriptor is for the reverse reaction, the slope has to 
            # be modified
            if rev:
                return self.slope*self._descriptor(**kwargs) + self.intercept
            else:
                return (self.slope-1.)*self._descriptor(**kwargs) \
                       + self.intercept
        else:
            if rev:
                return (self.slope-1.)*self._descriptor(**kwargs) \
                       + self.intercept
            else:
                return self.slope*self._descriptor(**kwargs) + self.intercept

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