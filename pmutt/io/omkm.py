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
    lines = [_get_file_timestamp(comment_char='# ')]

    '''Write units'''
    lines.extend(['', '#' + '-'*80, '# UNITS', '#' + '-'*80])
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