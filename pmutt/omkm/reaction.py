import numpy as np
import more_itertools as mit

from pmutt import _apply_numpy_operation
from pmutt import constants as c
from pmutt.cantera import _get_range_CTI
from pmutt.omkm.phase import InteractingInterface, StoichSolid
from pmutt.omkm.units import Units
from pmutt.reaction import Reaction
from pmutt.reaction.bep import BEP as BEP_parent


class SurfaceReaction(Reaction):
    """Expresses OpenMKM surface reaction in Cantera CTI format reaction.
    Inherits from :class:`~pmutt.reaction.Reaction`.
    
    Attributes
    ----------
        id : str, optional
            ID of the reaction. Default is None
        is_adsorption : bool, optional
            If True, the reaction represents an adsorption. Default is False
        beta : float, optional
            Power to raise the temperature in the rate expression. Default is 1
            if ``is_adsorption`` is False, 0 if ``is_adsorption`` is True.
        sticking_coeff : float, optional
            Sticking coefficient. Only relevant if ``is_adsorption`` is True.
            Default is 0.5 if ``is_adsorption`` is True, None if
            ``is_adsorption`` is False.
        direction : str, optional
            Direction of the reaction. Used for BEP relationships.
            Accepted options are 'cleavage' and 'synthesis'.
        kwargs : keyword arguments
            Keyword arguments used to initialize the reactants, transition
            state and products
    """
    def __init__(self,
                 id=None,
                 is_adsorption=False,
                 beta=None,
                 sticking_coeff=None,
                 direction=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.id = id
        self.is_adsorption = is_adsorption
        self.beta = beta
        self.direction = direction
        self.sticking_coeff = sticking_coeff

        # Assigns BEP to reaction
        if self.transition_state is not None:
            for species in self.transition_state:
                # Skip species that are not BEP relationships
                if not isinstance(species, BEP):
                    continue

                # Assign BEP relationship for easier reference
                self.bep = species

                # Generate ID from BEP if not previously assigned
                if self.id is None:
                    new_id = species._get_new_id(direction=self.direction)
                    self.id = '{}_{}_{:04d}'.format(species.name,
                                                    self.direction[:3], new_id)

                # Add reaction to BEP
                if self.direction == 'synthesis':
                    self.bep.synthesis_reactions.append(self)
                elif self.direction == 'cleavage':
                    self.bep.cleavage_reactions.append(self)

                # Assuming only one BEP would be assigned per reaction
                break
            else:
                self.bep = None

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, val):
        if isinstance(val, int) or isinstance(val, float):
            val = '{:04d}'.format(int(val))
        self._id = val

    @property
    def beta(self):
        return self._beta

    @beta.setter
    def beta(self, val):
        if val is None:
            if self.is_adsorption:
                val = 0.
            else:
                val = 1.
        self._beta = val

    @property
    def sticking_coeff(self):
        return self._sticking_coeff

    @sticking_coeff.setter
    def sticking_coeff(self, val):
        if val is None and self.is_adsorption:
            val = 0.5
        self._sticking_coeff = val

    def _get_n_surf(self):
        """Counts the number of surface reactants

        Returns
        -------
            n_surf : int
                Number of surface species
        """
        n_surf = 0
        for species, stoich in zip(self.reactants, self.reactants_stoich):
            if isinstance(species.phase, InteractingInterface):
                n_surf += stoich
        return n_surf

    def get_A(self,
              sden_operation='min',
              include_entropy=True,
              T=c.T0('K'),
              units='molec/cm2',
              **kwargs):
        """Calculates the preexponential factor in the Cantera format

        Parameters
        ----------
        sden_operation : str, optional
            Site density operation to use. Default is 'min'
        include_entropy : bool, optional
            If True, includes the entropy of activation. Default is True
        T : float, optional
            Temperature in K. Default is 298.15 K
        units : str or :class:`~pmutt.omkm.units.Units`, optional
            Units for A. If `Units` class specified, determines the units for A.
            Default is 'molec/cm2'
        kwargs : keyword arguments
            Parameters required to calculate pre-exponential factor
        """

        if self.transition_state is None or not include_entropy:
            A = c.kb('J/K') / c.h('J s')
        else:
            A = super().get_A(T=T, **kwargs) / T

        # Uses site with highest site density
        site_dens = []
        for reactant, stoich in zip(self.reactants, self.reactants_stoich):
            # Skip species without a catalyst site
            try:
                site_den = reactant.phase.site_density
            except AttributeError:
                continue
            site_dens.extend([site_den] * int(stoich))

        # Apply the operation to the site densities
        if len(site_dens) == 0:
            err_msg = ('At least one species requires a catalytic site with '
                       'site density to calculate A.')
            raise ValueError(err_msg)
        eff_site_den = _apply_numpy_operation(quantity=site_dens,
                                              operation=sden_operation,
                                              verbose=False)
        # Convert site density to appropriate unit
        if isinstance(units, Units):
            quantity_unit = units.quantity
            area_unit = '{}2'.format(units.length)
        else:
            quantity_unit, area_unit = units.split('/')
        eff_site_den = eff_site_den\
                       *c.convert_unit(initial='mol', final=quantity_unit)\
                       /c.convert_unit(initial='cm2', final=area_unit)

        n_surf = self._get_n_surf()
        A = A / eff_site_den**(n_surf - 1)
        return A

    def get_HoRT_act(self, rev=False, **kwargs):
        """Calculates the dimensionless enthalpy of activation. If there is no
        transition state species, calculates the delta dimensionless enthalpy

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            kwargs : keyword arguments
                Parameters required to calculate enthalpy of activation.
        Returns
        -------
            HoRT_act : float
                Dimensionless activation enthalpy. Returns the max of the
                following to ensure stable MKM performance:
                - Difference between reactants/products and the transition state
                - Difference between the reactants and the products
                - 0
        """
        act = self.transition_state is not None
        return np.max([
            0.,
            super().get_delta_HoRT(rev=rev, act=act, **kwargs),
            super().get_delta_HoRT(rev=rev, act=False, **kwargs)
        ])

    def get_H_act(self, units, T, rev=False, **kwargs):
        """Calculates the enthalpy of activation. If there is no transition
        state species, calculates the delta enthalpy

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            T : float
                Temperature in K
            kwargs : keyword arguments
                Parameters required to calculate enthalpy of activation.
        Returns
        -------
            HoRT_act : float
                Dimensionless activation enthalpy. Returns the max of the
                following to ensure stable MKM performance:
                - Difference between reactants/products and the transition state
                - Difference between the reactants and the products
                - 0
        """
        R_units = '{}/K'.format(units)
        return self.get_HoRT_act(rev=rev, T=T, **kwargs)*T*c.R(R_units)

    def get_GoRT_act(self, rev=False, act=False, **kwargs):
        """Calculates the dimensionless Gibbs energy of activation. If there is
        no transition state species, calculates the delta dimensionless
        Gibbs energy

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            act : bool, optional
                If True, uses the transition state as the final state. Default
                is False
            kwargs : keyword arguments
                Parameters required to calculate Gibbs energy.
        Returns
        -------
            GoRT_act : float
                Dimensionless Gibbs energy of activation. Returns the max of the
                following to ensure stable MKM performance:
                - Difference between reactants/products and the transition state
                - Difference between the reactants and the products
                - 0
        """
        act = self.transition_state is not None
        return np.max([
            0.,
            super().get_delta_GoRT(rev=rev, act=act, **kwargs),
            super().get_delta_GoRT(rev=rev, act=False, **kwargs)
        ])

    def get_G_act(self, units, T, P=1., rev=False, **kwargs):
        """Calculates the Gibbs energy of activation. If there is no transition
        state species, calculates the delta Gibbs energy

        Parameters
        ----------
            rev : bool, optional
                Reverse direction. If True, uses products as initial state
                instead of reactants. Default is False
            T : float
                Temperature in K
            P : float, optional
                Pressure in bar. Default is 1 bar
            kwargs : keyword arguments
                Parameters required to calculate enthalpy of activation.
        Returns
        -------
            HoRT_act : float
                Dimensionless activation enthalpy. Returns the max of the
                following to ensure stable MKM performance:
                - Difference between reactants/products and the transition state
                - Difference between the reactants and the products
                - 0
        """
        R_units = '{}/K'.format(units)
        return self.get_GoRT_act(rev=rev, T=T, P=P, **kwargs)*T*c.R(R_units)

    @classmethod
    def from_string(cls,
                    reaction_str,
                    species,
                    species_delimiter='+',
                    reaction_delimiter='=',
                    notes=None,
                    beta=1,
                    is_adsorption=False,
                    sticking_coeff=0.5,
                    direction=None,
                    id=None):
        """Create a reaction object using the reaction string

        Parameters
        ----------
            reaction_str : str
                Reaction string.
            species : dict
                Dictionary using the names as keys. If you have a list of
                species, use pmutt.pmutt_list_to_dict to make a dict.
            species_delimiter : str, optional
                Delimiter that separate species. Leading and trailing spaces
                will be trimmed. Default is '+'
            reaction_delimiter : str, optional
                Delimiter that separate states of the reaction. Leading and
                trailing spaces will be trimmed. Default is '='
            notes : str or dict, optional
                Other notes such as the source of the reaction. Default is None
            beta : float, optional
                Power to raise the temperature in the rate expression.
                Default is 1
            is_adsorption : bool, optional
                If True, the reaction represents an adsorption. Default is False
            sticking_coeff : float, optional
                Sticking coefficient. Only relevant if ``is_adsorption`` is
                True. Default is 0.5
            gas_phase : bool
                True if the reaction has only gas-phase species. This attribute
                is determined based on the reactants and products
        Returns
        -------
            SurfaceReaction : :class:`~pmutt.omkm.SurfaceReaction` object
        """
        rxn = Reaction.from_string(reaction_str=reaction_str,
                                   species=species,
                                   species_delimiter=species_delimiter,
                                   reaction_delimiter=reaction_delimiter)
        return cls(reactants=rxn.reactants,
                   reactants_stoich=rxn.reactants_stoich,
                   products=rxn.products,
                   products_stoich=rxn.products_stoich,
                   transition_state=rxn.transition_state,
                   transition_state_stoich=rxn.transition_state_stoich,
                   notes=notes,
                   beta=beta,
                   is_adsorption=is_adsorption,
                   sticking_coeff=sticking_coeff,
                   direction=direction,
                   id=id)

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
        return obj_dict

    def to_CTI(self,
               T=c.T0('K'),
               P=c.P0('bar'),
               quantity_unit='molec',
               length_unit='cm',
               act_energy_unit='cal/mol',
               ads_act_method='get_H_act',
               units=None):
        """Writes the object in Cantera's CTI format.

        Parameters
        ----------
            T : float, optional
                Temperature in K. Default is 298.15 K
            P : float, optional
                Pressure in bar. Default is 1 bar
            quantity_unit : str, optional
                Quantity unit to calculate A. Default is 'molec'
            length_unit : str, optional
                Length unit to calculate A. Default is 'cm'
            act_energy_unit : str, optional
                Unit to use for activation energy. Default is 'cal/mol'
            ads_act_method : str, optional
                Activation method to use for adsorption reactions. Accepted 
                options include 'get_H_act' and 'get_G_act'. Default is
                'get_H_act'.
            units : :class:`~pmutt.omkm.units.Units` object
                If specified, `quantity_unit`, `length_unit`, and
                `act_energy_unit` are overwritten. Default is None.
        Returns
        -------
            cti_str : str
                Surface reaction string in CTI format
        """
        if units is not None:
            quantity_unit = units.quantity
            length_unit = units.length
            act_energy_unit = units.act_energy

        reaction_str = self.to_string(stoich_space=True,
                                      species_delimiter=' + ',
                                      reaction_delimiter=' <=> ',\
                                      include_TS=False)
        # Determine the reaction IDs
        try:
            id = self.id
        except AttributeError:
            id_str = ''
        else:
            if id is None:
                id_str = ''
            else:
                id_str = ',\n                 id="{}"'.format(self.id)

        if self.is_adsorption:
            act_method = getattr(self, ads_act_method)
            act_val = act_method(units=act_energy_unit, T=T, P=P)
            cti_str = ('surface_reaction("{}",\n'
                       '                 stick({: .5e}, {}, {: .5e}){})'
                       ''.format(reaction_str, self.sticking_coeff, self.beta,
                                 act_val, id_str))
        else:
            A_units = '{}/{}2'.format(quantity_unit, length_unit)
            cti_str = ('surface_reaction("{}",\n'
                       '                 [{: .5e}, {}, {: .5e}]{})'
                       ''.format(reaction_str,
                                 self.get_A(T=T, P=P, include_entropy=False,
                                            units=A_units),
                                 self.beta,
                                 self.get_G_act(units=act_energy_unit, T=T,
                                                P=P),
                                 id_str))
        return cti_str


class BEP(BEP_parent):
    """Represents BEP relationships used by OpenMKM. Contains other attributes
    to aid in writing CTI file. Inherits from :class:`~pmutt.reaction.bep.BEP`
    
    Attributes
    ----------
        direction : str, optional
            Direction of the BEP. Accepted options are 'cleavage' and
            'synthesis'
        reactions : list of str or :class:`~pmutt.omkm.reaction.SurfaceReaction`
            Reactions associated with BEP relationship. Used for writing OpenMKM
            CTI files.
    """
    def __init__(self,
                 direction=None,
                 synthesis_reactions=[],
                 cleavage_reactions=[],
                 **kwargs):
        super().__init__(**kwargs)
        self.direction = direction
        self.synthesis_reactions = list(synthesis_reactions)
        self.cleavage_reactions = list(cleavage_reactions)

    def _get_bep_template(self, direction):
        """Get the BEP reation ID template

        Parameters
        ----------
            direction : str
                Direction of the new reaction. Either 'synthesis' or 'cleavage'.
        Returns
        -------
            bep_template : str
                BEP template
        """
        return '{}_{}'.format(self.name, direction[:3])

    def _get_reactions(self, direction):
        """Get the reactions in a particular direction

        Parameters
        ----------
            direction : str
                Direction of the new reaction. Either 'synthesis' or 'cleavage'.
        Returns
        -------
            reactions : list of Reaction
                Reactions in a particular direction.
        """
        # Get reactions in relevant direction
        direction = direction.lower()
        if direction == 'synthesis':
            reactions = self.synthesis_reactions
        elif direction == 'cleavage':
            reactions = self.cleavage_reactions
        else:
            err_msg = (
                'Invalid direction, {}, provided. See documentation for '
                'pmutt.omkm.reaction.BEP._get_reactions for '
                'supported options.'.format(direction))
            raise ValueError(err_msg)
        return reactions

    def _get_bep_reaction_ids(self, direction, format='int'):
        """Get the indices of the reactions that match the BEP template.

        Parameters
        ----------
            direction : str
                Direction of the new reaction. Either 'synthesis' or 'cleavage'.
            format : int
                Format to return ids. Accepted options are int and str. Default
                is int.
        Returns
        -------
            ids : list of int
                Integers of reactions that match BEP template.
        """
        reactions = self._get_reactions(direction=direction)
        bep_template = self._get_bep_template(direction=direction)
        reaction_ids = []
        for reaction in reactions:
            # Skip reactions that do not have the BEP template
            if not isinstance(reaction.id, str):
                continue
            if bep_template not in reaction.id:
                continue
            # Separate and record index from reactions with BEP ID template
            reaction_ids.append(reaction.id)
        if format == 'int':
            reaction_ids = [int(id.split('_')[-1]) for id in reaction_ids]
        return reaction_ids

    def _get_new_id(self, direction):
        """Gets the ID for the new reaction being assigned.
        
        Parameters
        ----------
            direction : str
                Direction of the new reaction. Either 'synthesis' or 'cleavage'.
        Returns
        -------
            new_id : int
                New ID to assign to reaction.
        """
        # Separate reactions with BEP-related ID
        reaction_ids = self._get_bep_reaction_ids(direction=direction)

        # New ID is incremented
        if len(reaction_ids) == 0:
            new_id = 1
        else:
            new_id = np.max(reaction_ids) + 1
        return new_id

    def _get_reactions_CTI(self, direction):
        """Returns the reactions names using the BEP using the following format:
        <bep_id>_<direction>_<reaction_id>
        
        Parameters
        ----------
            direction : str
                Direction of reaction. Accepted options are 'synthesis' or
                'cleavage'
            delimiter : str, optional
                Delimiter to separate reaction id. Default is '_'
        """
        reactions = self._get_reactions(direction=direction)
        if len(reactions) == 0:
            CTI_out = '[]'
        else:
            CTI_out = '['
            # Add reactions ranges with BEP template
            bep_template = self._get_bep_template(direction=direction)
            reaction_int_ids = self._get_bep_reaction_ids(direction=direction,
                                                          format='int')
            reaction_int_ids.sort()
            reaction_id_ranges = mit.consecutive_groups(reaction_int_ids)
            for id_range in reaction_id_ranges:
                id_range = list(id_range)
                if len(id_range) == 1:
                    CTI_range = '"{}_{:04d}", '.format(bep_template,
                                                       id_range[0])
                else:
                    CTI_range = ('"{0}_{1:04d} to {0}_{2:04d}", '
                                 ''.format(bep_template, id_range[0],
                                           id_range[-1]))
                CTI_out += CTI_range

            # Get reactions that do not have the BEP template
            for reaction in reactions:
                if isinstance(reaction.id,
                              str) and bep_template in reaction.id:
                    continue
                CTI_out += '"{}", '.format(reaction.id)

            CTI_out = '{}]'.format(CTI_out[:-2])
        return CTI_out

    def to_CTI(self, act_energy_unit=None, units=None, delimiter='_'):
        """Writes the object in Cantera's CTI format.

        Parameters
        ----------
            act_energy_unit : str, optional
                Unit to use for energy. Default is 'cal/mol'
            units : :class:`~pmutt.omkm.units.Units` object
                If specified, `energy_unit` is overwritten. Default is None.
        Returns
        -------
            cti_str : str
                Surface reaction string in CTI format
        """
        if units is not None:
            act_energy_unit = units.act_energy
        # synthesis_reactions = self._get_reactions_CTI(direction='synthesis')
        # cleavage_reactions = self._get_reactions_CTI(direction='cleavage')
        synthesis_reactions = _get_range_CTI(objs=self.synthesis_reactions,
                                             parent_obj=self,
                                             delimiter=delimiter)
        cleavage_reactions = _get_range_CTI(objs=self.cleavage_reactions,
                                            parent_obj=self,
                                            delimiter=delimiter)
        intercept = c.convert_unit(self.intercept, 'kcal/mol', act_energy_unit)
        cti_str = ('bep(id="{}",\n'
                   '    slope={},\n'
                   '    intercept={},\n'
                   '    direction="{}",\n'
                   '    cleavage_reactions={},\n'
                   '    synthesis_reactions={})\n'
                   ''.format(self.name, self.slope, intercept, self.direction,
                             cleavage_reactions, synthesis_reactions))
        return cti_str
