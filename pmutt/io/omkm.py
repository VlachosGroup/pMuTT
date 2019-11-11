import yaml

from pmutt import _force_pass_arguments
from pmutt.io import _get_file_timestamp
from pmutt.io.cantera import obj_to_CTI
from pmutt.io.ctml_writer import convert
from pmutt.omkm.units import Units


def write_cti(phases=None, species=None, reactions=None,
              lateral_interactions=None, units=None, filename=None,
              T=300., P=1., newline='\n', use_motz_wise=False):
    """Writes the units, phases, species, lateral interactions, reactions and 
    additional options in the CTI format for OpenMKM
    
    Parameters
    ----------
        phases : list of :class:`~pmutt.omkm.phase.Phase` objects
            Phases to write in CTI file. The species should already be assigned.
        species : list of :class:`~pmutt.empirical.nasa.Nasa`, :class:`~pmutt.empirical.nasa.Nasa9` or :class:`~pmutt.empirical.shomate.Shomate`
            Species to write in CTI file.
        reactions : list of :class:`~pmutt.omkm.reaction.SurfaceReaction`
            Reactions to write in CTI file.
        lateral_interactions : list of :class:`~pmutt.mixture.cov.PiecewiseCovEffect` objects, optional
            Lateral interactions to include in CTI file. Default is None.
        units : dict or :class:`~pmutt.omkm.units.Unit` object, optional
            Units to write file. If a dict is inputted, the key is the quantity and
            the value is the unit. If not specified, uses the default units of
            :class:`~pmutt.omkm.units.Unit`
        filename: str, optional
            Filename for the input.cti file. If not specified, returns file
            as str
        T : float, optional
            Temperature in K. Default is 300 K.
        P : float, optional
            Pressure in atm. Default is 1 atm.
    Returns
    -------
        lines_out : str
            If ``filename`` is None, CTI file is returned
    """
    lines = [_get_file_timestamp(comment_char='# '),
             '# See documentation for OpenMKM CTI file here:',
             '# https://vlachosgroup.github.io/openmkm/input']

    '''Write units'''
    lines.extend(['', '#' + '-'*79, '# UNITS', '#' + '-'*79])
    if units is None:
        units = Units()
    elif isinstance(units, dict):
        units = Units(**units)
    lines.append(units.to_CTI())

    '''Pre-assign IDs for lateral interactions so phases can be written'''
    if lateral_interactions is not None:
        lat_inter_lines = []
        i = 0
        if lateral_interactions is not None:
            for lat_interaction in lateral_interactions:
                if lat_interaction.name is None:
                    lat_interaction.name = '{:04d}'.format(i)
                    i += 1

                lat_inter_CTI = _force_pass_arguments(lat_interaction.to_CTI,
                                                      units=units)
                lat_inter_lines.append(lat_inter_CTI)

    '''Pre-assign IDs for reactions so phases can be written'''
    if reactions is not None:
        beps = []
        reaction_lines = []
        i = 0
        for reaction in reactions:
            # Assign reaction ID if not present
            if reaction.id is None:
                reaction.id = '{:04d}'.format(i)
                i += 1
            # Write reaction
            reaction_CTI = _force_pass_arguments(reaction.to_CTI, units=units)
            reaction_lines.append(reaction_CTI)

            # Add unique BEP relationship if any
            try:
                bep = reaction.bep
            except AttributeError:
                pass
            else:
                if bep is not None and bep not in beps:
                    beps.append(bep)

    '''Write phases'''
    if phases is not None:
        lines.extend(['', '#' + '-'*80, '# PHASES', '#' + '-'*80])
        for phase in phases:
            phase_CTI = _force_pass_arguments(phase.to_CTI, units=units)
            lines.append(phase_CTI)

    '''Write species'''
    if species is not None:
        lines.extend(['', '#' + '-'*80, '# SPECIES', '#' + '-'*80])
        for ind_species in species:
            ind_species_CTI = _force_pass_arguments(ind_species.to_CTI,
                                                    units=units)
            lines.append(ind_species_CTI)

    '''Write lateral interactions'''
    if lateral_interactions is not None:
        lines.extend(['', '#' + '-'*80, '# LATERAL INTERACTIONS', '#' + '-'*80])
        lines.extend(lat_inter_lines)

    if reactions is not None:
        '''Write reaction options'''
        lines.extend(['', '#' + '-'*80, '# REACTION OPTIONS', '#' + '-'*80])
        if use_motz_wise:
            lines.extend(['enable_motz_wise()\n'])
        else:
            lines.extend(['disable_motz_wise()\n'])

        '''Write reactions'''
        lines.extend(['', '#' + '-'*80, '# REACTIONS', '#' + '-'*80])
        lines.extend(reaction_lines)

        '''Write BEP Relationships'''
        if len(beps) > 0:
            lines.extend(['', '#' + '-'*80, '# BEP Relationships', '#' + '-'*80])
            # Only write each BEP once
            i = 0
            for bep in beps:
                bep_CTI = _force_pass_arguments(bep.to_CTI, units=units)
                # Increment counter if necessary
                if bep.name is None:
                    i += 1
                lines.append(bep_CTI)

    '''Write to file'''
    lines_out = '\n'.join(lines)
    if filename is not None:
        with open(filename, 'w', newline=newline) as f_ptr:
            f_ptr.write(lines_out)

        '''Write XML file'''
        xml_filename = '{}.xml'.format(filename.replace('.cti', ''))
        convert(filename=filename, outName=xml_filename)
    else:
        # Or return as string
        return lines_out

def write_yaml(reactor_type=None, mode=None, V=None, T=None, P=None, A=None,
               L=None, cat_abyv=None, flow_rate=None, end_time=None,
               transient=None, stepping=None, init_step=None, step_size=None,
               atol=None, rtol=None, phases=None, reactor={}, inlet_gas={},
               solver={}, simulation={}, misc={}, units=None, filename=None,
               yaml_options={'default_flow_style': False}):
    """Writes the reactor options in a YAML file for OpenMKM. 
    
    Parameters
    ----------
        reactor_type : str
            Type of reactor. Supported options include:
            
            - PFR
            - CSTR
            - Batch

            Value written to reactor.type.
        mode : str
            Operation of reactor. Supported options include:

            - Isothermal
            - Adiabatic

            Value written to reactor.mode.
        V : float or tuple
            Volume of reactor. Value written to reactor.volume. Units of
            length^3. See Notes section regarding unit specification.
        T : float or tuple
            Temperature of reactor. Value written to reactor.temperature. Units
            of temperature. See Notes section regarding unit specification.
        P : float or tuple
            Pressure of reactor. Value written to reactor.pressure. Units of
            pressure. See Notes section regarding unit specification.
        A : float or tuple
            Surface area of reactor. Value written to reactor.area. Units of
            length^2. See Notes section regarding unit specification.
        L : float or tuple
            Length of reactor. Value written to reactor.length. Units of length.
            See Notes section regarding unit specification.
        cat_abyv : float or tuple
            Catalyst surface area to volume ratio. Value written to
            reactor.cat_abyv. Units of 1/length. See Notes section regarding
            unit specification.
        flow_rate : float or tuple
            Volumetric flow rate of inlet. Value written to inlet_gas.flow_rate.
            Units of length^3/time. See Notes section regarding unit
            specification.
        end_time : float or tuple
            Reactor simulation time. For continuous reactors, the system is
            assumed to reach steady state by this time. Value written to
            simulation.end_time. Units of time. See Notes section regarding
            unit specification.
        transient : bool
            If True, transient operation results are saved. Otherwise,
            transient files are left blank. Value written to
            simulation.transient.
        stepping : str
            Time steps taken to simulate reactor. Supported options include:

            - logarithmic
            - regular

            Value written to simulation.stepping.
        init_step : float or tuple
            Initial step to take. Value written to simulation.init_step.
        step_size : float or tuple
            If ``stepping`` is logarithmic, represents the ratio between the
            next step and the current step. If ``stepping`` is regular,
            represents the time between the next step and the current step in
            units of time. Value written to simulation.step_size. See Notes
            section regarding unit specification.
        atol : float
            Absolute tolerance for solver. Value written to
            simulation.solver.atol.
        rtol : float
            Relative tolerance for solver. Value written to
            simulation.solver.rtol.
        phases : list of ``Phase`` objects
            Phases present in reactor. Each phase should have the ``name``
            and ``initial_state`` attribute.
        reactor : dict
            Generic dictionary for reactor to specify values not supported by
            ``write_yaml``.
        inlet_gas : dict
            Generic dictionary for inlet_gas to specify values not supported by
            ``write_yaml``.
        solver : dict
            Generic dictionary for solver to specify values not supported by
            ``write_yaml``.            
        simulation : dict
            Generic dictionary for simultaion to specify values not supported by
            ``write_yaml``.            
        units : dict or :class:`~pmutt.omkm.units.Unit` object, optional
            Units used for file. If a dict is inputted, the key is the quantity
            and the value is the unit. If not specified, uses the default units
            of :class:`~pmutt.omkm.units.Unit`.
        filename: str, optional
            Filename for the input.cti file. If not specified, returns file
            as str.
    Returns
    -------
        lines_out : str
            If ``filename`` is None, CTI file is returned
    Notes
    -----
        **Units**
        If ``units`` is not specified, all values in file are assumed to be SI
        units. If ``units`` is specified, all values inputted are assumed to
        be in the units specified by ``units``. If a particular unit set is
        desired for a value, a str can be inputted with the form
        quantity="<<value>> <<desired units>>" where ``value`` is a float-like
        and ``desired units`` is a string. See
        https://vlachosgroup.github.io/openmkm/input for the most up-to-date
        supported values.

        **Generic Dictionaries**
        Also, values in generic dictionaries (i.e. ``reactor``, ``inlet_gas``,
        ``simulation``) will ben written preferentially over arguments passed.
        i.e. The value in ``reactor['temperature']`` will be written instead
        of ``T``.


    """
    lines = [_get_file_timestamp(comment_char='# '),
            '# See documentation for OpenMKM YAML file here:',
            '# https://vlachosgroup.github.io/openmkm/input']

    '''Organize reactor parameters'''
    _assign_yaml_val('type', reactor, reactor_type)
    _assign_yaml_val('mode', reactor, mode)
    _assign_yaml_val('volume', reactor, V)
    _assign_yaml_val('temperature', reactor, T)
    _assign_yaml_val('pressure', reactor, P)
    _assign_yaml_val('area', reactor, A)
    _assign_yaml_val('length', reactor, L)
    _assign_yaml_val('cat_abyv', reactor, cat_abyv)

    '''Organize inlet gas parameters'''
    _assign_yaml_val('flow_rate', inlet_gas, flow_rate)

    '''Organize solver parameters'''
    _assign_yaml_val('atol', solver, atol)
    _assign_yaml_val('rtol', solver, rtol)

    '''Organize simulation parameters'''
    _assign_yaml_val('end_time', simulation, end_time)
    _assign_yaml_val('transient', simulation, transient)
    _assign_yaml_val('stepping', simulation, stepping)
    _assign_yaml_val('step_size', simulation, step_size)
    if len(solver) > 0:
        _assign_yaml_val('solver', simulation, solver)
    
    '''Organize phase parameters'''
    # TODO Iterate through phases and check for type to separate gas, bulk, surfaces
    # For each phase, assign name and initial conditions attribute
    # Sample dictionary:
    # phases = {'gas': [
    #                {'name': 'gas',
    #                 'initial_state': 'NH3: 1'}],
    #           'bulk': [
    #                {'name': 'bulk'}],
    #           'surfaces': [
    #               {'name': 'terrace',
    #                'initial_state': 'RU(S1):1'},
    #               {'name': 'step',
    #                'initial_state': 'Ru(S2):1'}]
    # }
    
    '''Assign misc values'''
    yaml_dict = misc.copy()
    
    '''Assign values to overall YAML dict'''
    headers_str = ('reactor', 'inlet_gas', 'simulation', 'phases')
    headers_dict = (reactor, inlet_gas, simulation, phases)
    for header_str, header_dict in zip(headers_str, headers_dict):
        if len(header_dict) > 0:
            yaml_dict[header_str] = header_dict

    '''Convert dictionary to YAML str'''
    yaml_str = yaml.dump(yaml_dict, **yaml_options)
    lines.append(yaml_str)

    '''Write to file'''
    lines_out = '\n'.join(lines)
    if filename is not None:
        with open(filename, 'w', newline=newline) as f_ptr:
            f_ptr.write(lines_out)
    else:
        # Or return as string
        return lines_out


def _assign_yaml_val(label, header, val=None):
    """Helper method to assign label to header
    
    Parameters
    ----------
    label : str
        Label name for attribute being assigned.
    header : dict
        Upper level dictionary that ``label`` will be nested under.
    val : obj, optional
        Value that is assigned to header[label] if the key does not exist or
        val is None.
    """
    if label not in header and val is not None:
        header[label] = val