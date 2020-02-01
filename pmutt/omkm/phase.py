from pmutt import constants as c
from pmutt.cantera import _get_range_CTI
import pmutt.cantera.phase as phase_cantera


class IdealGas(phase_cantera.IdealGas):
    """OpenMKM implementation of the ideal gas phase. Currently there are no
    differences between this class and :class:`~pmutt.cantera.phase.IdealGas`
    but one could add changes here in the future if necessary.
    """
    pass


class StoichSolid(phase_cantera.StoichSolid):
    """OpenMKM implementation of the stoichiometric solid phase. Currently there
    are no differences between this class and 
    :class:`~pmutt.cantera.phase.StoichSolid` but one could add changes here in
    the future if necessary.
    """
    pass


class InteractingInterface(phase_cantera.Phase):
    """OpenMKM implementation of interacting interface. Inherits from
    :class:`~pmutt.cantera.phase.Phase`.

    Attributes
    ----------
        name : str
            Name of the phase
        species : list of :class:`~pmutt._ModelBase` objects
            Species present in Phase
        site_density : float, optional
            Site density in g/cm2. Default is None
        phases : list of :class:`~pmutt.cantera.phase.Phase` objects or str
            Phases associated with this interface
        interactions : str, optional
            Source of lateral interactions. If any lateral interactions in CTI
            file occur in this phase, specify 'all'. Default is None.
        note : str, optional
            Comment field for users. Default is None.
        initial_state : None, optional
            Currently not supported. Gives ability to set initial temperature,
            pressure or other operating variables of the phase. Default is None.
        transport : None, optional
            Currently not supported. Gives ability to specify transport model
            to use for phase. Default is None.
        options : None, optional
            Currently not supported. Specify special options to the phase.
            Default is None.
    """
    def __init__(self,
                 name,
                 species=[],
                 site_density=None,
                 phases=None,
                 initial_state=None,
                 kinetics=None,
                 reactions=None,
                 transport=None,
                 interactions=None,
                 options=None,
                 note=None):
        super().__init__(name=name,
                         species=species,
                         kinetics=kinetics,
                         transport=transport,
                         options=options,
                         note=note,
                         reactions=reactions,
                         initial_state=initial_state)
        self.site_density = site_density
        self.phases = phases
        self.interactions = interactions

    @property
    def beps(self):
        if self.reactions is None:
            beps = None
        else:
            beps = []
            for reaction in self.reactions:
                # Skip if the reaction does not have a BEP
                try:
                    bep = reaction.bep
                except AttributeError:
                    continue
                # Skip if the reaction does not have a BEP
                if bep is None:
                    continue
                # Skip if this is not a unique BEP
                if bep.name in beps:
                    continue
                # Append if this is a valid BEP
                beps.append(bep.name)
        return beps

    def to_CTI(self,
               max_line_len=80,
               quantity_unit='molec',
               length_unit='cm',
               units=None,
               delimiter='_'):
        """Writes the object in Cantera's CTI format.

        Parameters
        ----------
            max_line_len : int, optional
                Maximum number of characters in the line. Default is 80.
            quantity_unit : str, optional
                Quantity unit to use to calculate A. Default is 'molec'
            length_unit : str, optional
                Length unit to use to calculate A. Default is 'cm'
            units : :class:`~pmutt.omkm.units.Units` object
                If specified, `quantity_unit` and `length_unit` are overwritten.
                Default is None.
            delimiter : str, optional
                Delimiter used to separate header from footer of reaction and
                lateral interaction IDs. Default is '_'.
        Returns
        -------
            CTI_str : str
                Object represented as a CTI string.
        """
        if units is not None:
            quantity_unit = units.quantity
            length_unit = units.length

        species_names = [species.name for species in self.species]
        area_unit = '{}2'.format(length_unit)
        site_den = self.site_density\
                   *c.convert_unit(initial='mol', final=quantity_unit)\
                   /c.convert_unit(initial='cm2', final=area_unit)

        phases_names = []
        for phase in self.phases:
            try:
                phases_names.append(phase.name)
            except AttributeError:
                phases_names.append(phase)
        # Add required fields
        cti_str = ('interacting_interface(name={},\n'
                   '                      elements={},\n'
                   '                      species={},\n'
                   '                      phases={},\n'
                   '                      site_density={},\n'.format(
                       phase_cantera.obj_to_CTI(self.name,
                                                line_len=max_line_len - 27,
                                                max_line_len=max_line_len - 28),
                       phase_cantera.obj_to_CTI(self.elements,
                                                line_len=max_line_len - 31,
                                                max_line_len=max_line_len),
                       phase_cantera.obj_to_CTI(species_names,
                                                line_len=max_line_len - 30,
                                                max_line_len=max_line_len),
                       phase_cantera.obj_to_CTI(phases_names,
                                                line_len=max_line_len - 29,
                                                max_line_len=max_line_len),
                       site_den))
        # Fields with ranges
        range_fields = ('interactions', 'reactions')
        for range_field in range_fields:
            val = getattr(self, range_field)
            # Skip empty fields
            if val is None:
                continue
            cti_str += ('                      {}={},\n'
                        ''.format(
                            range_field,
                            _get_range_CTI(objs=val,
                                           parent_obj=self,
                                           delimiter=delimiter)))

        # Add optional fields
        optional_fields = ('beps', 'transport', 'options', 'note')
        for field in optional_fields:
            val = getattr(self, field)
            # Skip empty fields
            if val is None:
                continue
            # Skip blank lists
            if len(val) == 0:
                continue

            cti_str += '                      {}={},\n'.format(
                field,
                phase_cantera.obj_to_CTI(val,
                                         max_line_len=max_line_len,
                                         line_len=max_line_len - len(field) -
                                         23))

        # Terminate the string
        cti_str = '{})\n'.format(cti_str[:-2])
        return cti_str
