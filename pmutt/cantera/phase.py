import more_itertools as mit

from pmutt import constants as c
from pmutt.io.cantera import obj_to_CTI

class Phase:
    """Parent class for Cantera phases

    Attributes
    ----------
        name : str
            Name of the phase
        species : list of :class:`~pmutt._ModelBase` objects
            Species present in Phase
    """
    def __init__(self, name, species=[], initial_state=None, kinetics=None,
                 transport=None, reactions=None, options=None, note=None):
        self.name = name
        self.species = species
        self.initial_state = initial_state
        self.kinetics = kinetics
        self.transport = transport
        self.reactions = reactions
        self.options = options
        self.note = note
    
    @property
    def elements(self):
        elements = set()
        for ind_species in self.species:
            elements |= set(ind_species.elements.keys())
        return elements

    @property
    def species_names(self):
        return [ind_species.name for ind_species in self.species]

    @property
    def species(self):
        return self._species

    @species.setter
    def species(self, val):
        for i in range(len(val)):
            val[i].phase = self
        self._species = val

    def append_species(self, val):
        self._species.append(val)
        self._species[-1].phase = self

    def extend_species(self, val):
        for i in range(len(val)):
            val[i].phase = self
        self._species.extend(val)

    def remove_species(self, name):
        i = self.index_species(name)
        self.pop_species(i)

    def index_species(self, name):
        return self.species_names.index(name)

    def pop_species(self, i):
        self._species.pop(i)

    def clear_species(self):
        self._species.clear()

    def copy_species(self):
        return self._species.copy()

class IdealGas(Phase):
    """Expresses ideal gas as Cantera CTI file. Inherits from
    :class:`~pmutt.cantera.phase.Phase`.

    Attributes
    ----------
        name : str
            Name of the phase
        species : list of :class:`~pmutt._ModelBase` objects
            Species present in Phase
        reactions : str, optional
            Source of reactions. If any reactions in CTI file occur in this
            phase, specify 'all'. Default is None.
        note : str, optional
            Comment field for users. Default is None.
        initial_state : None, optional
            Currently not supported. Gives ability to set initial temperature,
            pressure or other operating variables of the phase. Default is None.
        kinetics : None, optional
            Currently not supported. Gives ability to specify reaction kinetics
            to use for phase. Default is None.
        transport : None, optional
            Currently not supported. Gives ability to specify transport model
            to use for phase. Default is None.
        options : None, optional
            Currently not supported. Specify special options to the phase.
            Default is None.
    """
    def __init__(self, name, species=[], initial_state=None, kinetics=None,
                 reactions=None, transport=None, options=None, note=None):
        super().__init__(name=name, species=species, kinetics=kinetics,
                         transport=transport, options=options, note=note,
                         reactions=reactions, initial_state=initial_state)

    def to_CTI(self, max_line_len=80):
        """Writes the object in Cantera's CTI format.

        Parameters
        ----------
            max_line_len : int, optional
                Maximum number of characters in the line. Default is 80.
        Returns
        -------
            CTI_str : str
                Object represented as a CTI string.
        """
        species_names = [species.name for species in self.species]
        # Add required fields
        cti_str = ('ideal_gas(name={},\n'
                   '          elements={},\n'
                   '          species={},\n'.format(
                       obj_to_CTI(self.name, line_len=max_line_len-15,
                                  max_line_len=max_line_len-16),
                       obj_to_CTI(self.elements, line_len=max_line_len-19,
                                  max_line_len=max_line_len),
                       obj_to_CTI(species_names, line_len=max_line_len-18,
                                  max_line_len=max_line_len)))
        # Add optional fields
        optional_fields = ('kinetics', 'transport', 'options', 'note',
                           'reactions')
        for field in optional_fields:
            val = getattr(self, field)
            # Skip empty fields
            if val is None:
                continue

            cti_str += '          {}={},\n'.format(
                    field, obj_to_CTI(val, line_len=max_line_len-len(field)-11,
                    max_line_len=max_line_len))

        # Terminate the string
        cti_str = '{})\n'.format(cti_str[:-2])
        return cti_str

class StoichSolid(Phase):
    """Expresses stoichiometric solid as Cantera CTI file. Inherits from
    :class:`~pmutt.cantera.phase.Phase`.

    Attributes
    ----------
        name : str
            Name of the phase
        species : list of :class:`~pmutt._ModelBase` objects
            Species present in Phase
        density : float, optional
            Bulk density in g/cm3. Default is None
        initial_state : None, optional
            Currently not supported. Gives ability to set initial temperature,
            pressure or other operating variables of the phase. Default is None.
        transport : None, optional
            Currently not supported. Gives ability to specify transport model
            to use for phase. Default is None.
        options : None, optional
            Currently not supported. Specify special options to the phase.
            Default is None.
        note : str, optional
            Comment field for users. Default is None.
    """
    def __init__(self, name, species=[], initial_state=None, transport=None,
                 options=None, density=None, note=None):
        super().__init__(name=name, species=species, transport=transport,
                         options=options, note=note,
                         initial_state=initial_state)
        self.density = density

    def to_CTI(self, max_line_len=80, mass_unit='g', length_unit='cm',
               units=None):
        """Writes the object in Cantera's CTI format.

        Parameters
        ----------
            max_line_len : int, optional
                Maximum number of characters in the line. Default is 80.
            mass_unit : str, optional
                Mass unit for `density`. Default is 'g'
            length_unit : str, optional
                Length unit for `density`. Default is 'cm'
            units : :class:`~pmutt.cantera.units.Units` object, optional
                If specified, `mass_unit` and `length_unit` are overwritten.
                Default is None.
        Returns
        -------
            CTI_str : str
                Object represented as a CTI string.
        """
        if units is not None:
            length_unit = units.length
            mass_unit = units.mass

        species_names = [species.name for species in self.species]
        volume_unit = '{}3'.format(length_unit)
        density = self.density*c.convert_unit(initial='g', final=mass_unit)\
                  /c.convert_unit(initial='cm3', final=volume_unit)
        # Add required fields
        cti_str = ('stoichiometric_solid(name={},\n'
                   '                     elements={},\n'
                   '                     species={},\n'
                   '                     density={},\n'.format(
                       obj_to_CTI(self.name, line_len=max_line_len-26,
                                  max_line_len=max_line_len-27),
                       obj_to_CTI(self.elements,line_len=max_line_len-30,
                                  max_line_len=max_line_len),
                       obj_to_CTI(species_names, line_len=max_line_len-29,
                                  max_line_len=max_line_len),
                       density))
        # Add optional fields
        optional_fields = ('transport', 'options', 'note', 'initial_state')
        for field in optional_fields:
            val = getattr(self, field)
            # Skip empty fields
            if val is None:
                continue

            cti_str += '                     {}={},\n'.format(
                    field, obj_to_CTI(val, line_len=max_line_len-len(field)-22,
                    max_line_len=max_line_len))

        # Terminate the string
        cti_str = '{})\n'.format(cti_str[:-2])
        return cti_str