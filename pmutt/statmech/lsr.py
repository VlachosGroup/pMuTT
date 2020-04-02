# -*- coding: utf-8 -*-
from pmutt import constants as c
from pmutt import _force_pass_arguments, _ModelBase
from pmutt.statmech import StatMech, ConstantMode, presets
from pmutt.reaction import Reaction
from pmutt.io.json import remove_class, json_to_pmutt


class LSR(_ModelBase):
    """Represents a linear scaling relationship

    :math:`\\Delta E^{AH_x} = \\alpha \\Delta E^{A} + \\beta`
        
    :math:`E^{AH_{x}^*}=\\Delta E^{AH_x} + E^* + E^{AH_{x(g)}}`

    where:

    +--------------------------+----------------------------------------------------------+----------+
    | Symbol                   | Description                                              | Units    |
    +==========================+==========================================================+==========+
    | :math:`\\Delta E^{AH_x}`  | Binding energy of adsorbate, AHx, on interested surface. | kcal/mol |
    +--------------------------+----------------------------------------------------------+----------+
    | :math:`\\Delta E^{A}`     | Binding energy of adsorbate, A, on interested surface.   | kcal/mol |
    +--------------------------+----------------------------------------------------------+----------+
    | :math:`\\alpha`           | Slope of LSR                                             | \\-       |
    +--------------------------+----------------------------------------------------------+----------+
    | :math:`\\beta`            | Intercept of LSR                                         | kcal/mol |
    +--------------------------+----------------------------------------------------------+----------+
    | :math:`E^*`              | Electronic energy of interested surface.                 | kcal/mol |
    +--------------------------+----------------------------------------------------------+----------+
    | :math:`E^{AH_{x(g)}}`    | Electronic energy of AHx in the gas phase.               | kcal/mol |
    +--------------------------+----------------------------------------------------------+----------+

    Attributes
    ----------
        slope : float
            Slope of LSR relationship. Represents alpha in above equation
        intercept : float
            Intercept of LSR relationship in kcal/mol. Represents beta in above
            equation
        reaction : float or :class:`~pmutt.reaction.Reaction` object
            Reaction to use to calculate reference binding energy. Binding
            energy calculated using ``get_delta_E`` and if that fails,
            ``get_delta_H``. If the binding energy is specified as a float
            (in kcal/mol), it will be converted to a
            :class:`~pmutt.reaction.Reaction` made of
            :class:`~pmutt.statmech.StatMech` objects
        surf_species : float or :class:`~pmutt.statmech.StatMech` object, optional
            Surface species. If the surface's energy is specified as a float (in
            kcal/mol), it will be converted to a
            :class:`~pmutt.statmech.StatMech` object with a single
            :class:`~pmutt.statmech.ConstantMode` mode. Default is 0.
        gas_species : float or :class:`~pmutt.statmech.StatMech` object, optional
            Gas-phase species. If the gas' energy is specified as a float (in
            kcal/mol), it will be converted to a
            :class:`~pmutt.statmech.StatMech` object with a single
            :class:`~pmutt.statmech.ConstantMode` mode. Default is 0.
        notes : str or dict, optional
            Extra notes such as the source of the LSR. Default is None
    Notes
    -----
        If the LSR relationship represents the binding of an adsorbate (AHx) on
        surface (M2) relative to a reference surface (M1) (see equation below),
        you may still use this class.

        :math:`\\Delta E^{AH_x}_{M_2} = \\alpha (\\Delta E^{A}_{M_2}
        - \\Delta E^{A}_{M_1}) + \\Delta E^{AH_x}_{M_1}`

        The ``reaction`` inputted will represent the following stoichiometric
        reaction.

        :math:`*_{(M_1)} + A^*_{(M_2)} \\leftrightarrow A^*_{(M_1)} + *_{(M_2)}`

    """
    def __init__(self,
                 slope,
                 intercept,
                 reaction,
                 surf_species=0.,
                 gas_species=0.,
                 notes=None):
        self.slope = slope
        self.intercept = intercept
        self.reaction = reaction
        self.surf_species = surf_species
        self.gas_species = gas_species
        self.notes = notes

    @property
    def reaction(self):
        return self._reaction
    
    @reaction.setter
    def reaction(self, val):
        if isinstance(val, Reaction):
            reaction = val
        else:
            reaction = _float_to_reaction(val)
        self._reaction = reaction

    @property
    def gas_species(self):
        return self._gas_species
    
    @gas_species.setter
    def gas_species(self, val):
        if isinstance(val, _ModelBase):
            gas_species = val
        else:
            gas_species = _float_to_species(val)
        self._gas_species = gas_species

    @property
    def surf_species(self):
        return self._surf_species
    
    @surf_species.setter
    def surf_species(self, val):
        if isinstance(val, _ModelBase):
            surf_species = val
        else:
            surf_species = _float_to_species(val)
        self._surf_species = surf_species

    def get_CvoR(self):
        """Calculates the dimensionless heat capacity at constant volume.

        :math:`\\frac{C_V^{LSR}}{R}=0`
        """
        return 0.

    def get_CpoR(self):
        """Calculates the dimensionless heat capacity at constant pressure.

        :math:`\\frac{C_V^{LSR}}{R}=0`
        """
        return self.get_CvoR()

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
            E_surf = self.surf_species.get_E(**kwargs)
        except AttributeError:
            E_surf = self.surf_species.get_H(**kwargs)
        # Calculate gas-phase energy
        try:
            E_gas = self.gas_species.get_E(**kwargs)
        except AttributeError:
            E_gas = self.gas_species.get_H(**kwargs)
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
            'reaction': self.reaction.to_dict(),
            'surf_species': self.surf_species.to_dict(),
            'gas_species': self.gas_species.to_dict(),
            'notes': self.notes
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
            LSR : LSR object
        """
        json_obj = remove_class(json_obj)
        json_obj['reaction'] = json_to_pmutt(json_obj['reaction'])
        json_obj['surf_species'] = json_to_pmutt(json_obj['surf_species'])
        json_obj['gas_species'] = json_to_pmutt(json_obj['gas_species'])
        return cls(**json_obj)

def _float_to_species(val):
    """Helper method to convert a float to a :class:`~pmutt.statmech.StatMech`
    object

    Parameters
    ----------
        val : float
            Value (in kcal/mol)
    Returns
    -------
        obj : :class:`~pmutt.statmech.StatMech` object
            :class:`~pmutt.statmech.StatMech` object that gives the val
            when `get_E` is called
    """
    val = val * c.convert_unit(initial='kcal/mol', final='eV/molecule')
    return StatMech(U=val, H=val, F=val, G=val, **presets['constant'])

def _float_to_reaction(val):
    """Helper method to convert a float to a :class:`~pmutt.reaction.Reaction`
    object

    Parameters
    ----------
        val : float
            Value (in kcal/mol)
    Returns
    -------
        obj : :class:`~pmutt.reaction.Reaction` object
            :class:`~pmutt.reaction.Reaction` object that gives the val
            when `get_delta_E` is called
    """
    reactant = StatMech()
    product = self._float_to_specie(val=val)
    return Reaction(reactants=[reactant],
                    reactants_stoich=[1.],
                    products=[product],
                    products_stoich=[1.])
