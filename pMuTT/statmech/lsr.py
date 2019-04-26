# -*- coding: utf-8 -*-
from pMuTT import constants as c
from pMuTT import _force_pass_arguments, _pMuTTBase
from pMuTT.statmech import StatMech, ConstantMode, presets
from pMuTT.reaction import Reaction
from pMuTT.io.json import remove_class, json_to_pMuTT


class LSR(_pMuTTBase):
    """Represents a linear scaling relationship

    :math:`\\Delta E^{AH_x} = \\alpha E^{A} + \\beta`
        
    :math:`E^{AH_{x}^*}=\\Delta E^{AH_x} + E^* + E^{AH_{x(g)}}`

    Attributes
    ----------
        slope : float
            Slope of LSR relationship. Represents alpha in above equation
        intercept : float
            Intercept of LSR relationship in kcal/mol. Represents beta in above
            equation
        reaction : float or :class:`~pMuTT.reaction.Reaction` object
            Reaction to use to calculate reference binding energy. Binding
            energy calculated using `get_delta_E` and if that fails,
            `get_delta_H`. If the binding energy is specified as a float
            (in kcal/mol), it will be converted to a
            :class:`~pMuTT.reaction.Reaction` made of `~pMuTT.statmech.StatMech`
            objects
        surf_specie : float or :class:`~pMuTT.statmech.StatMech` object, optional
            Surface specie. If the surface's energy is specified as a float (in
            kcal/mol), it will be converted to a
            :class:`~pMuTT.statmech.StatMech` object with a single
            :class:`~pMuTT.statmech.ConstantMode` mode. Default is 0.
        gas_specie : float or :class:`~pMuTT.statmech.StatMech` object, optional
            Gas-phase specie. If the gas' energy is specified as a float (in
            kcal/mol), it will be converted to a
            :class:`~pMuTT.statmech.StatMech` object with a single
            :class:`~pMuTT.statmech.ConstantMode` mode. Default is 0.
        notes : str or dict, optional
            Extra notes such as the source of the LSR. Default is None
    """
    def __init__(self, slope, intercept, reaction, surf_specie=0.,
                 gas_specie=0., notes=None):
        self.slope = slope
        self.intercept = intercept
        if not isinstance(reaction, Reaction):
            reaction = self._float_to_reaction(reaction)
        self.reaction = reaction
        if not isinstance(surf_specie, _pMuTTBase):
            surf_specie = self._float_to_specie(surf_specie)
        self.surf_specie = surf_specie
        if not isinstance(gas_specie, _pMuTTBase):
            gas_specie = self._float_to_specie(gas_specie)
        self.gas_specie = gas_specie
        self.notes = notes

    def get_UoRT(self, T=c.T0('K'), **kwargs):
        """Calculates the dimensionless internal energy using the LSR
        relationship
        
        Parameters
        ----------
            T : float, optional
                Temperature in K. Default is 298.15 K
            kwargs : keyword arguments
                Parameters to calculate reference binding energy, surface
                energy and gas-phase energy
        Returns
        -------
            UoRT : float
                Dimensionless internal energy
        """
        kwargs['T'] = T
        kwargs['units'] = 'kcal/mol'
        # Calculate reference binding energy
        try:
            deltaE_ref = self.reaction.get_delta_E(**kwargs)
        except AttributeError:
            deltaE_ref = self.reaction.get_delta_H(**kwargs)
        # Calculate surface energy
        try:
            E_surf = self.surf_specie.get_E(**kwargs)
        except AttributeError:
            E_surf = self.surf_specie.get_H(**kwargs)
        # Calculate gas-phase energy
        try:
            E_gas = self.gas_specie.get_E(**kwargs)
        except AttributeError:
            E_gas = self.gas_specie.get_H(**kwargs)
        return (self.slope*deltaE_ref + self.intercept + E_surf + E_gas) \
               /c.R('kcal/mol/K')/T

    def get_HoRT(self, T=c.T0('K'), **kwargs):
        """Calculates the dimensionless enthalpy using the LSR
        relationship
        
        Parameters
        ----------
            T : float, optional
                Temperature in K. Default is 298.15 K
            kwargs : keyword arguments
                Parameters to calculate reference binding energy, surface
                energy and gas-phase energy
        Returns
        -------
            HoRT : float
                Dimensionless enthalpy
        """
        return self.get_UoRT(T=T, **kwargs)

    def get_SoR(self):
        """Calculates the dimensionless entropy using the LSR
        relationship
        
        Returns
        -------
            SoR : float
                Dimensionless entropy. Since LSR handles binding energies,
                always returns 0
        """
        return 0.

    def get_FoRT(self, T=c.T0('K'), **kwargs):
        """Calculates the dimensionless Helmholtz energy using the LSR
        relationship
        
        Parameters
        ----------
            T : float, optional
                Temperature in K. Default is 298.15 K
            kwargs : keyword arguments
                Parameters to calculate reference binding energy, surface
                energy and gas-phase energy
        Returns
        -------
            FoRT : float
                Dimensionless Helmholtz energy
        """
        return self.get_UoRT(T=T, **kwargs) - self.get_SoR()

    def get_GoRT(self, T=c.T0('K'), **kwargs):
        """Calculates the dimensionless Gibbs energy using the LSR
        relationship
        
        Parameters
        ----------
            T : float, optional
                Temperature in K. Default is 298.15 K
            kwargs : keyword arguments
                Parameters to calculate reference binding energy, surface
                energy and gas-phase energy
        Returns
        -------
            GoRT : float
                Dimensionless Gibbs energy
        """
        return self.get_HoRT(T=T, **kwargs) - self.get_SoR()

    def _float_to_specie(self, val):
        """Converts a float to a :class:`~pMuTT.statmech.StatMech` object

        Parameters
        ----------
            val : float
                Value (in kcal/mol)
        Returns
        -------
            obj : :class:`~pMuTT.statmech.StatMech` object
                :class:`~pMuTT.statmech.StatMech` object that gives the val
                when `get_E` is called
        """
        val = val*c.convert_unit(initial='kcal/mol', final='eV/molecule')
        return StatMech(U=val, H=val, F=val, G=val, **presets['constant'])

    def _float_to_reaction(self, val):
        """Converts a float to a :class:`~pMuTT.reaction.Reaction` object

        Parameters
        ----------
            val : float
                Value (in kcal/mol)
        Returns
        -------
            obj : :class:`~pMuTT.reaction.Reaction` object
                :class:`~pMuTT.reaction.Reaction` object that gives the val
                when `get_delta_E` is called
        """
        reactant = StatMech()
        product = self._float_to_specie(val=val)
        return Reaction(reactants=[reactant], reactants_stoich=[1.],
                        products=[product], products_stoich=[1.])

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = {'class': str(self.__class__),
                    'slope': self.slope,
                    'intercept': self.intercept,
                    'reaction': self.reaction.to_dict(),
                    'surf_specie': self.surf_specie.to_dict(),
                    'gas_specie': self.gas_specie.to_dict(),
                    'notes': self.notes}
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
            LSR : LSR object
        """
        json_obj = remove_class(json_obj)
        json_obj['reaction'] = json_to_pMuTT(json_obj['reaction'])
        json_obj['surf_specie'] = json_to_pMuTT(json_obj['surf_specie'])
        json_obj['gas_specie'] = json_to_pMuTT(json_obj['gas_specie'])
        return cls(**json_obj)