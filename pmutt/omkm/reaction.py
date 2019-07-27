import numpy as np
from pmutt import _apply_numpy_operation
from pmutt import constants as c
from pmutt.reaction import Reaction
from pmutt.omkm.phase import InteractingInterface, StoichSolid

class SurfaceReaction(Reaction):
    """Cantera reaction. Has additional attributes to support input and output

    Attributes
    ----------
        beta : float, optional
            Power to raise the temperature in the rate expression. Default is 1
        is_adsorption : bool, optional
            If True, the reaction represents an adsorption. Default is False
        sticking_coeff : float, optional
            Sticking coefficient. Only relevant if ``is_adsorption`` is True.
            Default is 0.5
        kwargs : keyword arguments
            Keyword arguments used to initialize the reactants, transition
            state and products
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

    def get_A(self, sden_operation='min', include_entropy=True, T=c.T0('K'),
              units='molecule/cm2', **kwargs):
        """Calculates the preexponential factor in the Cantera format

        Parameters
        ----------
        sden_operation : str, optional
            Site density operation to use. Default is 'min'
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

        # Uses site with highest site density
        site_dens = []
        for reactant, stoich in zip(self.reactants, self.reactants_stoich):
            # Skip species without a catalyst site
            try:
                site_den = reactant.phase.site_density
            except AttributeError:
                continue
            site_dens.extend([site_den]*int(stoich))
        
        # Apply the operation to the site densities
        if len(site_dens) == 0:
            raise ValueError('At least one species requires a catalytic site '
                             'with site density to calculate A.')
        eff_site_den = _apply_numpy_operation(quantity=site_dens,
                                              operation=sden_operation,
                                              verbose=False)
        # Convert site density to appropriate unit
        quantity_unit, area_unit = units.split('/')
        eff_site_den = eff_site_den\
                       *c.convert_unit(initial='mol', final=quantity_unit)\
                       /c.convert_unit(initial='cm2', final=area_unit)

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
        act = self.transition_state is not None
        return np.max([0.,
                       super().get_delta_HoRT(rev=rev, act=act, **kwargs),
                       super().get_delta_HoRT(rev=rev, act=False, **kwargs)])

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
        act = self.transition_state is not None
        return np.max([0., 
                       super().get_delta_GoRT(rev=rev, act=act, **kwargs),
                       super().get_delta_GoRT(rev=rev, act=False, **kwargs)])

    @classmethod
    def from_string(cls, reaction_str, species, species_delimiter='+',
                    reaction_delimiter='=', bep_descriptor=None, notes=None,
                    beta=1, is_adsorption=False, sticking_coeff=0.5):
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
            bep_descriptor : str, optional
                If the transition state is a :class:`~pmutt.reaction.bep.BEP`
                object, the descriptor can be set using this parameter. See
                :class:`~pmutt.reaction.bep.BEP` documentation for accepted
                descriptors
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
        rxn = super().from_string(reaction_str=reaction_str,
                                  species=species,
                                  species_delimiter=species_delimiter,
                                  reaction_delimiter=reaction_delimiter,
                                  bep_descriptor=bep_descriptor)
        return cls(reactants=rxn.reactants,
                   reactants_stoich=rxn.reactants_stoich,
                   products=rxn.products, products_stoich=rxn.products_stoich,
                   transition_state=rxn.transition_state,
                   transition_state_stoich=rxn.transition_state_stoich,
                   bep_descriptor=bep_descriptor,
                   notes=notes, beta=beta, is_adsorption=is_adsorption,
                   sticking_coeff=sticking_coeff)

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

    def to_CTI(self, T=c.T0('K'), P=c.P0('bar'), quantity_unit='molecule',
               length_unit='cm', act_energy_unit='J/mol'):
        """Writes the string for a single surface reaction reaction

        Parameters
        ----------
            T : float, optional
                Temperature in K. Default is 298.15 K
            P : float, optional
                Pressure in bar. Default is 1 bar
            act_energy_unit : str, optional
                Unit to use for activation energy. Default is 'J/mol'
        Returns
        -------
            cti_str : str
                Surface reaction string in CTI format
        """
        reaction_str = self.to_string(stoich_space=True,
                                      species_delimiter=' + ',
                                      reaction_delimiter=' <=> ',\
                                      include_TS=False)

        if self.is_adsorption:
            cti_str = 'surface_reaction("{}", stick({: .5e}, {}, {: .5e}))'.format(
                            reaction_str, self.sticking_coeff, self.beta,
                            self.get_G_act(units=act_energy_unit, T=T, P=P))
        else:
            A_units = '{}/{}2'.format(quantity_unit, length_unit)
            cti_str = 'surface_reaction("{}", [{: .5e}, {}, {: .5e}])'.format(
                            reaction_str,
                            self.get_A(T=T, P=P, include_entropy=False,
                                       units=A_units),
                            self.beta,
                            self.get_G_act(units=act_energy_unit, T=T, P=P))
        return cti_str