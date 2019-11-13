import yaml

from pmutt import _force_pass_arguments, _is_iterable
from pmutt.io import _get_file_timestamp
from pmutt.io.cantera import obj_to_CTI
from pmutt.io.ctml_writer import convert
from pmutt.cantera.phase import IdealGas, StoichSolid
from pmutt.omkm.phase import InteractingInterface
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
               atol=None, rtol=None, phases=None, reactor=None, inlet_gas=None,
               multi_T=None, multi_P=None, multi_flow_rate=None,
               output_format=None, solver=None, simulation=None,
               multi_input=None, misc=None, units=None, filename=None, 
               yaml_options={'default_flow_style': False, 'indent': 4},
               newline='\n'):
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
        V : float or str
            Volume of reactor. Value written to reactor.volume. Units of
            length^3. See Notes section regarding unit specification.
        T : float
            Temperature (in K) of reactor. Value written to reactor.temperature.
        P : float or str
            Pressure of reactor. Value written to reactor.pressure. Units of
            pressure. See Notes section regarding unit specification.
        A : float or str
            Surface area of reactor. Value written to reactor.area. Units of
            length^2. See Notes section regarding unit specification.
        L : float or str
            Length of reactor. Value written to reactor.length. Units of length.
            See Notes section regarding unit specification.
        cat_abyv : float or str
            Catalyst surface area to volume ratio. Value written to
            reactor.cat_abyv. Units of 1/length. See Notes section regarding
            unit specification.
        flow_rate : float or str
            Volumetric flow rate of inlet. Value written to inlet_gas.flow_rate.
            Units of length^3/time. See Notes section regarding unit
            specification.
        end_time : float or str
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
        init_step : float or str
            Initial step to take. Value written to simulation.init_step.
        step_size : float or str
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
        multi_T : list of float
            Multiple temperatures (in K) of reactor. Value written to
            simulation.multi_input.temperature.
        multi_P : list of float/str
            Multiple pressures of reactor. Value written to
            simulation.multi_input.pressure. Units of pressure. See Notes
            section regarding unit specification.
        multi_flow_rate : list of float/str
            Multiple flow rates to run the model. Value written to
            simulation.multi_input.flow_rate. Units of length3/time. See Notes
            section regarding unit specification.
        output_format : str
            Format for output files. Supported options include:
            
            - CSV
            - DAT

            Value written to simulation.output_format.

        misc : dict
            Generic dictionary for any parameter specified at the top level.
        solver : dict
            Generic dictionary for solver to specify values not supported by
            ``write_yaml``.            
        simulation : dict
            Generic dictionary for simultaion to specify values not supported by
            ``write_yaml``.            
        multi_input : dict
            Generic dictionary for multi_input to specify values not supported
            by ``write_yaml``.
        units : dict or :class:`~pmutt.omkm.units.Unit` object, optional
            Units used for file. If a dict is inputted, the key is the quantity
            and the value is the unit. If not specified, uses the default units
            of :class:`~pmutt.omkm.units.Unit`.
        filename: str, optional
            Filename for the input.cti file. If not specified, returns file
            as str.
        yaml_options : dict
            Options to pass when converting the parameters to YAML format. See
            `PyYAML documentation`_ for ``dump`` for available options.
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
        quantity="<<value>> <<desired units>>" where ``value`` is float-like
        and ``desired units`` is a string. For example, flow_rate="1 cm3/s".
        See https://vlachosgroup.github.io/openmkm/input for the most up-to-date
        supported values.

        **Generic Dictionaries**
        Also, values in generic dictionaries (i.e. ``reactor``, ``inlet_gas``,
        ``simulation``) will ben written preferentially over arguments passed.
        i.e. The value in ``reactor['temperature']`` will be written instead
        of ``T``.

    .. _`PyYAML Documentation`: https://pyyaml.org/wiki/PyYAMLDocumentation
    """
    lines = [_get_file_timestamp(comment_char='# '),
            '# See documentation for OpenMKM YAML file here:',
            '# https://vlachosgroup.github.io/openmkm/input']

    '''Initialize units'''
    if isinstance(units, dict):
        units = Units(**units)

    '''Organize reactor parameters'''
    if reactor is None:
        reactor = {}
    reactor_params = [('type', reactor_type, None), ('mode', mode, None),
                      ('volume', V, '_length3'),
                      ('area', A, '_length2'),
                      ('length', L, '_length'),
                      ('cat_abyv', cat_abyv, '/_length')]
    # Process temperature
    if T is not None:
        reactor_params.append(('temperature', T, None))
    elif multi_T is not None:
        reactor_params.append(('temperature', multi_T[0], None))
    # Process pressure
    if P is not None:
        reactor_params.append(('pressure', P, '_pressure'))
    elif multi_P is not None:
        reactor_params.append(('pressure', multi_P[0], '_pressure'))
    
    for parameter in reactor_params:
        _assign_yaml_val(parameter[0], parameter[1], parameter[2], reactor,
                         units)

    '''Organize inlet gas parameters'''
    if inlet_gas is None:
        inlet_gas = {}
    inlet_gas_params = []
    # Process flow rate
    if flow_rate is not None:
        inlet_gas_params.append(('flow_rate', flow_rate, '_length3/_time'))
    elif multi_flow_rate is not None:
        inlet_gas_params.append(('flow_rate', multi_flow_rate[0],
                                 '_length3/_time'))

    for parameter in inlet_gas_params:
        _assign_yaml_val(parameter[0], parameter[1], parameter[2], inlet_gas,
                         units)

    '''Organize solver parameters'''
    if solver is None:
        solver = {}
    solver_params = [('atol', atol, None), ('rtol', rtol, None)]
    for parameter in solver_params:
        _assign_yaml_val(parameter[0], parameter[1], parameter[2], solver,
                         units)

    '''Organize multi parameters'''
    if multi_input is None:
        multi_input = {}
    # Ensure multi quantities are in correct type
    if _is_iterable(multi_T):
        multi_T = list(multi_T)
    if _is_iterable(multi_P):
        multi_P = list(multi_P)
    if _is_iterable(multi_flow_rate):
        multi_flow_rate = list(multi_flow_rate)
    multi_input_params = [('temperature', multi_T, None),
                          ('pressure', multi_P, '_pressure'),
                          ('flow_rate', multi_flow_rate, '_length3/_time')]
    for parameter in multi_input_params:
        _assign_yaml_val(parameter[0], parameter[1], parameter[2], multi_input,
                         units)

    '''Organize simulation parameters'''
    if simulation is None:
        simulation = {}
    simulation_params = (('end_time', end_time, '_time'),
                         ('transient', transient, None),
                         ('stepping', stepping, None),
                         ('step_size', step_size, '_time'),
                         ('init_step', init_step, None),
                         ('output_format', output_format, None))
    for parameter in simulation_params:
        _assign_yaml_val(parameter[0], parameter[1], parameter[2], simulation,
                         units)
    if len(solver) > 0:
        _assign_yaml_val('solver', solver, None, simulation, units)
    if len(multi_input) > 0:
        _assign_yaml_val('multi_input', multi_input, None, simulation, units)
    
    '''Organize phase parameters'''
    if phases is not None:
        if isinstance(phases, dict):
            phases_dict = phases.copy()
        elif isinstance(phases, list):
            phases_dict = {}
            for phase in phases:
                phase_info = {'name': phase.name}

                # Assign intial state if available
                if phase.initial_state is not None:
                    initial_state_str = '"'
                    for species, mole_frac in phase.initial_state.items():
                        initial_state_str += '{}:{}, '.format(species,mole_frac)
                    initial_state_str = '{}"'.format(initial_state_str[:-2])
                    phase_info['initial_state'] = initial_state_str

                # Assign phase type
                if isinstance(phase, IdealGas):
                    phase_type = 'gas'
                elif isinstance(phase, StoichSolid):
                    phase_type = 'bulk'
                elif isinstance(phase, InteractingInterface):
                    phase_type = 'surfaces'

                try:
                    phases_dict[phase_type].append(phase_info)
                except KeyError:
                    phases_dict[phase_type] = [phase_info]
        # If only one entry for phase type, reassign it.
        for phase_type, phases in phases_dict.items():
            if len(phases) == 1:
                phases_dict[phase_type] = phases[0]
    '''Assign misc values'''
    if misc is None:
        misc = {}
    yaml_dict = misc.copy()
    
    '''Assign values to overall YAML dict'''
    headers = (('reactor', reactor), ('inlet_gas', inlet_gas),
               ('simulation', simulation), ('phases', phases_dict))
    for header in headers:
        if len(header[1]) > 0:
            yaml_dict[header[0]] = header[1]

    '''Convert dictionary to YAML str'''
    yaml_str = yaml.dump(yaml_dict, **yaml_options)
    # Remove redundant quotes
    yaml_str = yaml_str.replace('\'', '')
    lines.append(yaml_str)

    '''Write to file'''
    lines_out = '\n'.join(lines)
    if filename is not None:
        with open(filename, 'w', newline=newline) as f_ptr:
            f_ptr.write(lines_out)
    else:
        # Or return as string
        return lines_out


def _assign_yaml_val(label, val, val_units, header, units=None):
    """Helper method to assign label to header
    
    Parameters
    ----------
    label : str
        Label name for attribute being assigned.
    val : obj
        Value that is assigned to header[label] if the key does not exist or
        val is None.
    val_units : str
        Units for ``val`` where quantities are proceeded by '_'.
        e.g. '_length3/_time' for the volumetric flow rate.
    header : dict
        Upper level dictionary that ``label`` will be nested under.
    units : :class:`~pmutt.omkm.units.Unit` object, optional
        Units to write file.
    """
    # Do nothing if the label was previously assigned
    if label in header:
        return
    # Do nothing if value is not specified
    if val is None:
        return
    if isinstance(val, str):
        val = '\"{}\"'.format(val)
    # Assign the value as is if the unit type is None
    if val_units is None:
        header[label] = val
        return
    # Assume SI units if units is not specified
    if units is None:
        header[label] = val
        return
    # If the value is numerical and units were specified, add the units
    if isinstance(val, (int, float)):
        val_str = '\"{} {}\"'.format(val, val_units)
        for unit_type, unit in units.__dict__.items():
            val_str = val_str.replace('_{}'.format(unit_type), unit)
        header[label] = val_str
    # If the value is a list and units were specified, add units to each entry
    elif isinstance(val, list):
        vals_list = ['\"{} {}\"'.format(i, val_units) for i in val]
        for unit_type, unit in units.__dict__.items():
            old_str = '_{}'.format(unit_type)
            for i, val in enumerate(vals_list):
                vals_list[i] = val.replace(old_str, unit)
        header[label] = vals_list