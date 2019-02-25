# -*- coding: utf-8 -*-
"""
pMuTT.io_.chemkin

Reads reactions lists from Chemkin gas.inp and surf.inp files
"""
import re
import numpy as np
import datetime
from pMuTT import constants as c
from pMuTT import pMuTT_list_to_dict, _force_pass_arguments

def read_reactions(filename, species=None):
    """Directly read reactions from Chemkin gas.inp or surf.inp files

    Parameters
    ----------
        filename : str
            Input filename for Chemkin surf or gas .inp file
        species : list of :class:`~pMuTT.empirical.nasa.Nasa` objects
            List of NASA objects containing thermodynamic properties for
            all Reactants and Products in Reactions
            default = None. Will not return React_obj and Prod_obj
    Returns
    -------
        Reactions   : list of reactions
        Reactants   : list of reactants found in reactions
        React_obj   : list of NASA polynomials for each Reactant
                      If species object list is supplied
        React_stoic : list of reaction stoichiometries for Reactants
        Products    : list of products found in reactions
        Prod_obj    : list of NASA polynomials for each Product
                      If species object list is supplied
        Prod_stoic  : list of reaction stoichiometries for Products
    Raises
    ------
        FileNotFoundError
            If the surf.inp or gas.inp file isn't found.
        NameError
            If the species file does not exist
        AttributeError
            If the species list is incorrect format
    """
    if species is not None:
        species_dict = pMuTT_list_to_dict(species)

    rxns = []
    with open(filename, 'r') as lines:
        for line in lines:
            if re.findall(r'(^[^\!].+)( *<*(?<![0-9][eE])[=\-]>* *)', line):
                rxns.append(line.strip())
    RHS = []
    LHS = []
    for rxn in rxns:
        LHS.append(re.split(r' *<*(?<![0-9][eE])[=\-]>* *', rxn)[0])
        RHS.append(re.split(r' *<*(?<![0-9][eE])[=\-]>* *', rxn)[1])
    Reactants = []
    Products = []
    React_obj = []
    Prod_obj = []
    React_stoic = []
    Prod_stoic = []
    for Reacs, Prods in zip(LHS, RHS):
        Reactants.append(re.split(r' *\+ *| +', Reacs))
        Products.append(re.split(r' *\+ *| +', Prods)[0:-3])
        R = []
        RS = []
        Rx = []
        for RR in Reactants[-1]:
            stoic = re.findall(r'^[0-9]*', RR)[0]
            if stoic == '':
                stoic = 1
            else:
                RR = RR.replace(stoic, "")
                stoic = int(stoic)
            Rx.append(RR)
            RS.append(stoic)
            if species is not None:
                R.append(species_dict[RR])
        Reactants[-1] = Rx
        React_stoic.append(RS)
        P = []
        PS = []
        Px = []
        for PP in Products[-1]:
            stoic = re.findall(r'^[0-9]*', PP)[0]
            if stoic == '':
                stoic = 1
            else:
                PP = PP.replace(stoic, "")
                stoic = int(stoic)
            Px.append(PP)
            PS.append(stoic)
            if species is not None:
                P.append(species_dict[PP])
        Products[-1] = Px
        Prod_stoic.append(PS)
        React_obj.append(R)
        Prod_obj.append(P)
    Reactions = []
    for rxn, Prods in zip(rxns, Products):
        Reactions.append(rxn[0:rxn.rindex(Prods[-1]) + len(Prods[-1])])
    if species is not None:
        return(Reactions, Reactants, React_obj, React_stoic,
               Products, Prod_obj, Prod_stoic)
    else:
        return(Reactions, Reactants, React_stoic, Products, Prod_stoic)

def write_EA(reactions, conditions, write_gas_phase=False, 
             filename='EAs.inp', act_method_name='get_EoRT_act', 
             float_format=' .2E', species_delimiter='+', 
             reaction_delimiter='<=>', stoich_format='.0f', newline='\n',
             column_delimiter='  '):
    """Writes the EAs.inp or EAg.inp file for Chemkin

    Parameters
    ---------- 
        reactions : list of :class:`~pMuTT.reaction.ChemkinReaction` objects
            Reactions to write
        conditions : list of dicts
            Conditions to evaluate each reaction. The key of the dictionaries
            should be a relevant quantity to evaluate the reaction (e.g. T, P)
        write_gas_phase : bool, optional
            If True, only gas phase reactions are written (including
            adsorption). If False, only surface phase reactions are written.
            Default is False
        filename : str, optional
            Filename for the EAs.inp file. Default is 'EAs.inp'
        method_name : str, optional
            Method to use to calculate values. Typical values are:

            - 'get_EoRT_act' (default)
            - 'get_HoRT_act'
            - 'get_GoRT_act'

        float_format : float, optional
            Format to write numbers. Default is ' .2E' (scientific notation
            rounded to 2 decimal places with a preceding space if the value is
            positive)
        stoich_format : float, optional
            Format to write stoichiometric numbers. Default is '.0f' (integer)
        newline : str, optional
            Newline character to use. Default is the Linux newline character
        column_delimiter : str, optional
            Delimiter for columns. Default is '  '
    """
    valid_reactions = []
    for reaction in reactions:
        # Skip gas-phase reactions if we want surface phase
        if not write_gas_phase and reaction.gas_phase:
            continue
        # Skip surface-phase reactions if we want gas phase
        if write_gas_phase and not reaction.gas_phase:
            continue
        valid_reactions.append(reaction)
        

    n_reactions = len(valid_reactions)
    # Add initial comment help line and number of reactions
    lines = [
        _get_filestamp(),
        ('!The first line is the number of reactions. Subsequent lines follow ' 
         'the format'),
        ('!of rxn (from surf.out) followed by the EA/RT value at each run '
         'condition.'),
        ('!There may be one slight deviation from surf.out: any repeated '
         'species should'),
        ('!be included in the reaction string with a stoichiometric '
         'coefficient equal to'),
        ('!the number of times the species appears in the reaction. '
         'If not using'),
        '!MultiInput, then only the first value is used.',
        '  {}  !Number of reactions'.format(n_reactions)]

    # Find length to pad reactions
    max_rxn_len = _get_max_reaction_len(reactions=valid_reactions,
                                        species_delimiter=species_delimiter,
                                        reaction_delimiter=reaction_delimiter,
                                        stoich_format=stoich_format,
                                        include_TS=False)
    rxn_padding = max_rxn_len + len(column_delimiter)

    # Define string formats to use
    float_field = '{:%s}' % float_format
    str_field = '{:%d}' % rxn_padding

    column_line = _write_column_line(padding=rxn_padding, 
                                     column_delimiter=column_delimiter,
                                     float_format=float_format,
                                     n_conditions=len(conditions))
    lines.append(column_line)

    # Add line for each reaction step
    for reaction in valid_reactions:
        line = [str_field.format(
                    reaction.to_string(species_delimiter=species_delimiter,
                                       reaction_delimiter=reaction_delimiter,
                                       stoich_format=stoich_format,
                                       include_TS=False))]
        for condition in conditions:
            method = getattr(reaction, act_method_name)
            quantity = _force_pass_arguments(method, **condition)
            line.append(float_field.format(quantity))
        lines.append(column_delimiter.join(line))
    lines.append('EOF')
    with open(filename, 'w', newline=newline) as f_ptr:
        f_ptr.write('\n'.join(lines))

def write_gas(nasa_species, filename='gas.inp', T=c.T0('K'), reactions=[],
              species_delimiter='+', reaction_delimiter='=',
              act_method_name='get_E_act', act_unit='kcal/mol', 
              float_format=' .3E', stoich_format='.0f', newline='\n',
              column_delimiter='  ', **kwargs):
    """Writes the gas.inp Chemkin file.

    Parameters
    ----------
        nasa_species : list of :class:`~pMuTT.empirical.nasa.Nasa objects
            Surface and gas species used in Chemkin mechanism. Used to write
            elements section
        filename : str, optional
            File name for gas.inp file. Default is 'gas.inp'
        reactions : :class:`~pMuTT.reaction.Reactions` object, optional
            Reactions in mechanism. Reactions with only gas-phase species will
            be written to this file
    """
    # Get unique elements and gas-phase species
    unique_elements = set()
    gas_species = []
    for specie in nasa_species:
        for element in specie.elements.keys():
            unique_elements.add(element)
        if specie.phase.upper() == 'G':
            gas_species.append(specie.name)
    unique_elements = list(unique_elements)

    # Get gas-phase reactions
    gas_reactions = [reaction for reaction in reactions if reaction.gas_phase]
    reaction_lines = _write_reaction_lines(
            reactions=gas_reactions, species_delimiter=species_delimiter,
            reaction_delimiter=reaction_delimiter, include_TS=False,
            stoich_format=stoich_format, act_method_name=act_method_name,
            act_unit=act_unit, float_format=float_format, 
            column_delimiter=column_delimiter, T=T, sden_operation=None,
            **kwargs)
    # Collect all the lines into a list
    lines = [
        _get_filestamp(),
        '!Elements present in gas and surface species',
        'ELEMENTS']
    lines.extend(unique_elements)
    lines.extend(['END',
                  '',
                  '!Gas-phase species',
                  'SPECIES'])
    lines.extend(gas_species)
    lines.extend(['END',
                  '',
                  '!Gas-phase reactions. The rate constant expression is:',
                  '!k = kb/h * (T)^beta * exp(-Ea/RT)',
                  '!Each line has 4 columns:',
                  '!- Reaction reactants and products separated by <=>',
                  '!- Preexponential factor, kb/h',
                  '!- Beta (power to raise T in rate constant expression)',
                  ('!- Ea (Activation Energy or Gibbs energy of activation in '
                   'kcal/mol'),
                  'REACTIONS'])
    lines.extend(reaction_lines)
    lines.append('END')

    with open(filename, 'w', newline=newline) as f_ptr:
        f_ptr.write('\n'.join(lines))

def write_surf(nasa_species, sden_operation='sum', 
               filename='surf.inp', T=c.T0('K'), reactions=[],
               species_delimiter='+', reaction_delimiter='=',
               act_method_name='get_E_act', act_unit='kcal/mol', 
               float_format=' .3E', stoich_format='.0f', newline='\n', 
               column_delimiter='  ', use_mw_correction=True, **kwargs):
    """Writes the surf.inp Chemkin file
    
    Parameters
    ----------
        nasa_species : list of :class:`~pMuTT.empirical.nasa.Nasa` objects
            Surface used in Chemkin mechanism. If gas-phase species are present,
            they will be ignored
        filename : str, optional
            Filename for surf.inp file. Default is 'surf.inp'
        T : float, optional
            Temperature to calculate activation energy. Default is 298.15 K
        reactions : list of :class:`~pMuTT.reaction.ChemkinReaction` objects
            Chemkin reactions to write in surf.inp file. Purely gas-phase
            reactions will be ignored
        species_delimiter : str, optional
            Delimiter to separate species when writing reactions. Default is '+'
        reaction_delimiter : str, optional
            Delimiter to separate reaction sides. Default is '='
        act_method_name : str, optional
            Name of method to use to calculate activation function
        act_unit : str, optional
            Units to calculate activation energy. Default is 'kcal/mol'
        float_format : str, optional
            String format to print floating numbers. Default is ' .3E' (i.e.
            scientific notation rounded to 3 decimal places with a leading space
            for positive numbers)
        stoich_format : str, optional
            String format to print stoichiometric coefficients. Default is '.0f'
            (i.e. integers)
        newline : str, optional
            Newline character. Default is the Linux newline character
        column_delimiter : str, optional
            Delimiter to separate columns. Default is '  '
        use_mw_correction : bool, optional
            If True, uses the Motz-Wise corrections. Default is True
        kwargs : keyword arguments
            Parameters needed to calculate activation energy and preexponential
            factor
    """
    # Organize species by their catalyst sites
    cat_adsorbates = {}
    unique_cat_sites = []
    for specie in nasa_species:
        # Skip gas phase species
        if specie.phase.upper() == 'G':
            continue
        # Skip the bulk species
        if specie.cat_site.bulk_specie == specie.name:
            continue

        cat_name = specie.cat_site.name
        try:
            cat_adsorbates[cat_name].append(specie)
        except KeyError:
            cat_adsorbates[cat_name] = [specie]
            unique_cat_sites.append(specie.cat_site)
    
    cat_site_lines = []
    for cat_site in unique_cat_sites:
        # Add catalyst site header
        cat_site_name = '{}/'.format(cat_site.name)
        cat_site_lines.append('SITE/{:<14}SDEN/{:.5E}/'.format(
                cat_site_name, cat_site.site_density))
        cat_site_lines.append('')
        # Add species adsorbed on that site
        for specie in cat_adsorbates[cat_site.name]:
            cat_site_lines.append('{}{}/{}/'.format(column_delimiter, 
                                                    specie.name, 
                                                    specie.n_sites))
        cat_site_lines.append('')
    # Write bulk species
    for cat_site in unique_cat_sites:
        # Add bulk line
        cat_site_lines.append('BULK {}/{:.1f}/'.format(cat_site.bulk_specie, 
                                                       cat_site.density))

    # Get surface reaction lines
    surf_reactions = \
            [reaction for reaction in reactions if not reaction.gas_phase]
    reaction_lines = _write_reaction_lines(
                reactions=surf_reactions, species_delimiter=species_delimiter,
                reaction_delimiter=reaction_delimiter, include_TS=False,
                stoich_format=stoich_format, act_method_name=act_method_name,
                act_unit=act_unit, float_format=float_format, 
                column_delimiter=column_delimiter, T=T, 
                sden_operation=sden_operation, **kwargs)

    # Preparing reaction line header
    mw_field = '{:<5}'
    if use_mw_correction:
        mw_str = mw_field.format('MWON')
    else:
        mw_str = mw_field.format('MWOFF')

    if 'oRT' in act_method_name:
        act_unit_str = ''
    else:
        act_unit_str = act_unit.upper()

    lines = [
        _get_filestamp(),
        '!Surface species',
        '!Each catalyst site has the following format:',
        '!SITE/[Site name]/      SDEN/[Site density in mol/cm2]/',
        '![Adsorbate Name]/[# of Sites occupied]/ (for every adsorbate)',
        '!BULK [Bulk name]/[Bulk density in g/cm3]',
    ]
    lines.extend(cat_site_lines)
    lines.extend(['END',
                  '',
                  '!Gas-phase reactions.',
                  '!The reaction line has the following format:',
                  '!REACTIONS  MW[ON/OFF]   [Ea units]',
                  '!where MW stands for Motz-Wise corrections and if the Ea',
                  ('!units are left blank, then the activation energy should '
                   'be dimensionless (i.e. E/RT)'),
                  '!The rate constant expression is:',
                  '!k = kb/h/site_den^(n-1) * (T)^beta * exp(-Ea/RT)',
                  ('!where site_den is the site density and is the number '
                   'of surface species (including empty sites)'),
                  '!Each line has 4 columns:',
                  '!- Reaction reactants and products separated by =',
                  '!- Preexponential factor, kb/h/site_den^(n-1), or ',
                  '!  sticking coefficient if adsorption reaction',
                  '!- Beta (power to raise T in rate constant expression)',
                  ('!- Ea (Activation Energy or Gibbs energy of activation in '
                   'specified units'),
                  ('!Adsorption reactions can be represented using the STICK '
                   'keyword'),
                  'REACTIONS{2}{0}{2}{1}'.format(mw_str, act_unit_str, 
                                                 column_delimiter)])
    lines.extend(reaction_lines)
    lines.append('END')

    # Write the file
    with open(filename, 'w', newline=newline) as f_ptr:
        f_ptr.write('\n'.join(lines))

def write_T_flow(T=None, P=None, Q=None, abyv=None, conditions=None, 
                 filename='T_flow.inp', float_format='.3E', newline='\n',
                 column_delimiter='  '):
    """Writes the T_flow.inp Chemkin file

    Parameters
    ----------
    T : list of float
        Temperatures in K. If not specified
    P : list of float
        Pressures in atm. If not specified, data will be taken from
        conditions
    Q : list of float
        Volumetric flow rate in cm3/s
    filename : str, optional
        Name of the file. Default is T_flow.inp
    float_format : str, optional
        Format to write floating point numbers. Default is '.3E' (i.e.
        scientific notation rounded to 3 decimal places)
    newline : str, optional
        Newline character. Default is the Linux newline character
    column_delimiter : str, optional
        Delimiter for columns. Default is '  '
    """
    line_field = '{:%s}{}{:%s}{}{:%s}{}{:%s}  !{:<3}' % (float_format,
                                                         float_format,
                                                         float_format,
                                                         float_format)
    lines = [
        _get_filestamp(),
        '!Conditions for each reaction run',
        '!Only used when MultiInput in tube.inp is set to "T"',
        '!T[K]    {0}P[atm]   {0}Q[cm3/s] {0}abyv[cm-1]{0}Run #'.format(
            column_delimiter)]
    for i, (T_i, P_i, Q_i, abyv_i) in enumerate(zip(T, P, Q, abyv)):
        lines.append(line_field.format(T_i, column_delimiter, 
                                       P_i, column_delimiter,
                                       Q_i, column_delimiter,
                                       abyv_i, i+1))
    with open(filename, 'w', newline=newline) as f_ptr:
        f_ptr.write('\n'.join(lines))

def write_tube_mole(mole_frac_conditions, nasa_species, 
                    filename='tube_mole.inp', float_format=' .3f', newline='\n',
                    column_delimiter='  '):
    """Write tube_mole.inp Chemkin file
    
    Parameters
    ----------
    mole_frac_conditions : list of dict
        Each dictionary should have the keys as the name of the species and the
        value as the initial mole fraction
    nasa_species : list of :class:`~pMuTT.empirical.nasa.Nasa` objects
        Nasa species to find phase information
    filename : str, optional
        Name of the file. Default is 'tube_mole.inp'
    float_format : str, optional
        Format to write floating point numbers. Default is '.3E' (i.e.
        scientific notation rounded to 3 decimal places)
    newline : str, optional
        Newline character. Default is the Linux newline character
    column_delimiter : str, optional
        Delimiter for columns. Default is '  '
    """
    # Prepare the float field for printing mole fractions
    float_field = '{:%s}' % float_format

    # Get relevant species
    unique_specie_names = set()
    for mole_fracs in mole_frac_conditions:
        for specie_name in mole_fracs.keys():
            unique_specie_names.add(specie_name)
    # Convert nasa_species to a dict for quicker lookup
    unique_species = [specie for specie in nasa_species 
                      if specie.name in unique_specie_names]

    # Find length to pad species
    max_species_len = _get_max_species_len(species=unique_species,
                                           ignore_gas_phase=False,
                                           include_phase=True)
    species_padding = max_species_len + len(column_delimiter)

    column_line = _write_column_line(padding=species_padding, 
                                     column_delimiter=column_delimiter,
                                     float_format=float_format,
                                     n_conditions=len(mole_frac_conditions))
    species_lines = []
    for specie in unique_species:
        specie_line = _get_specie_str(specie=specie, 
                                      include_phase=True).ljust(species_padding)
        print(specie_line)
        for condition in mole_frac_conditions:
            # If the mole fraction was not specified, assumed to be 0
            try:
                float_str = float_field.format(condition[specie.name])
            except KeyError:
                float_str = float_field.format(0.)
            specie_line = '{}{}{}'.format(specie_line, column_delimiter, 
                                          float_str)
        species_lines.append(specie_line)


    lines = [
        _get_filestamp(),
        "!Specify the 'species/phase/' pair /(in quotes!)/ and the associated",
        ("!composition values. If the composition does not sum to 1 for each "
         "phase or"),
        ("!site type, it will be renormalized to 1. At the end of a " 
         "calculation, a"),
        ("!file containing the complete composition and mass flux "
         "(the last entry) will"),
        ("!be generated. This file's format is completely compatible with the "
         "current"),
        "!input file and can be used to restart that calculation.",
        ("0       itube_restart -- will be >0 if a restart file is used or 0 "
         "for the first run"),
        '{:<3}    Number of nonzero species'.format(len(unique_species)),
        column_line,
    ]
    lines.extend(species_lines)
    lines.append('EOF')

    with open(filename, 'w', newline=newline) as f_ptr:
        f_ptr.write('\n'.join(lines))


def _get_filestamp():
    """Writes a comment indicating when the file was written by pMuTT"""
    return '!File generated by pMuTT on {}'.format(datetime.datetime.now())

def _get_max_reaction_len(reactions, species_delimiter, reaction_delimiter,
                          stoich_format, include_TS):
    """Returns the maximum length of the reactions

    Parameters
    ----------
        reactions : :class:`~pMuTT.reaction.Reactions` object
            Reactions to check
        species_delimiter : str
            Separates species
        reaction_delimiter : str
            Separates reaction states
        stoich_format : float
            Format to write stoichiometric numbers
        include_TS : bool
            If True, includes transition states in output
    Returns
    -------
        max_len : int
            Maximum length of reactions
    """
    rxns_len = []
    for reaction in reactions:
        rxn_len = len(reaction.to_string(species_delimiter=species_delimiter,
                                         reaction_delimiter=reaction_delimiter,
                                         stoich_format=stoich_format,
                                         include_TS=include_TS))
        rxns_len.append(rxn_len)
    if len(rxns_len)==0:
        max_len = 0
    else:
        max_len = max(rxns_len)
    return max_len

def _get_specie_str(specie, include_phase):
    """Gets the specie string

    Parameters
    ----------
        specie : :class:`~pMuTT.empirical.nasa.Nasa` object
            Specie to use
        include_phase : bool, optional
            If True, includes the reaction phase. Default is True
    Returns
    -------
        specie_str : str
            Species string
    """
    if include_phase:
        if specie.phase.upper() == 'G':
            phase = '/GAS/'
        elif specie.cat_site is not None:
            phase = '/{}/'.format(specie.cat_site.name)
        else:
            raise ValueError('Must specify phase or catalyst site to include '
                             'phase')
    else:
        phase = ''
    return "'{}{}'".format(specie.name, phase)

def _get_max_species_len(species, ignore_gas_phase=False, include_phase=True):
    """Calculate the maximum species length

    Parameters
    ----------
        species : list of :class:`~pMuTT.empirical.nasa.Nasa` objects
            Species to consider
        ignore_gas_phase : bool, optional
            If True, ignores the gas-phase species. Default is False
        include_phase : bool, optional
            If True, includes the reaction phase. Default is True
    Returns
    -------
        max_len : int
            Maximum length of species
    """
    species_len = []
    for specie in species:
        if ignore_gas_phase and specie.phase == 'G':
            continue
        specie_str = _get_specie_str(specie=specie, include_phase=include_phase)
        species_len.append(len(specie_str))
    return np.max(species_len)

def _write_reaction_lines(reactions, species_delimiter, reaction_delimiter, 
                          include_TS, stoich_format, act_method_name, act_unit, 
                          float_format, column_delimiter, sden_operation,
                          **kwargs):
    """Write the reaction lines in the Chemkin format

    Parameters
    ----------
        reactions : list of :class:`~pMuTT.reaction.ChemkinReaction` objects
            Chemkin reactions to write in surf.inp file
        species_delimiter : str
            Delimiter to separate species when writing reactions
        reaction_delimiter : str
            Delimiter to separate reaction sides
        act_method_name : str
            Name of method to use to calculate activation function
        act_unit : str
            Units to calculate activation energy
        float_format : str
            String format to print floating numbers
        stoich_format : str
            String format to print stoichiometric coefficients
        column_delimiter : str
            Delimiter to separate columns
        kwargs : keyword arguments
            Parameters needed to calculate activation energy and preexponential
            factor
    Returns
    -------
        reaction_lines : str
            Reactions represented in Chemkin format
    """
    max_reaction_len = _get_max_reaction_len(
            reactions=reactions, species_delimiter=species_delimiter,
            reaction_delimiter=reaction_delimiter, stoich_format=stoich_format,
            include_TS=include_TS)
    float_field = '{:%s}' % float_format

    reaction_lines = []
    for reaction in reactions:
        # Get reaction string
        reaction_str = reaction.to_string(
                species_delimiter=species_delimiter,
                reaction_delimiter=reaction_delimiter,
                stoich_format=stoich_format,
                include_TS=include_TS).ljust(max_reaction_len)
        # Calculate preexponential factor
        if reaction.is_adsorption:
            A = reaction.sticking_coeff
        else:
            # If using delta_G, take out entropic contribution in A
            if 'get_delta_G' in act_method_name:
                include_entropy = False
            else:
                include_entropy = True
            A = reaction.get_A(include_entropy=include_entropy,
                               sden_operation=sden_operation, **kwargs)
        A_str = float_field.format(A)

        # Format beta value
        beta_str = float_field.format(reaction.beta)

        # Calculate activation energy
        if act_method_name!='get_EoRT_act' and act_method_name!='get_E_act':
            kwargs['activation'] = True
        kwargs['units'] = act_unit
        Ea_method = getattr(reaction, act_method_name)
        try:
            Ea = _force_pass_arguments(Ea_method, **kwargs)
        except AttributeError:
            Ea = 0.
        Ea_str = float_field.format(Ea)

        reaction_line = '{0}{4}{1}{4}{2}{4}{3}'.format(reaction_str, A_str,
                                                       beta_str, Ea_str,
                                                       column_delimiter)
        if reaction.is_adsorption:
            reaction_line = '{}\nSTICK'.format(reaction_line)
        reaction_lines.append(reaction_line)
    return reaction_lines

def _write_column_line(padding, column_delimiter, float_format,
                       n_conditions):
    """Writes the run # in column format
    
    Parameters
    ----------
        padding : str
            Padding to offset first entry
        column_delimiter : str
            Column delimiter
        float_format : str
            Floating point format to use
        n_conditions : int
            Number of conditions
    Returns
    -------
        column_line : str
            Column line to be printed
    """
    float_field = '{:%s}' % float_format
    column_line = '!'.ljust(padding)
    float_field_len = len(float_field.format(np.pi)+column_delimiter)
    column_field = '{:> %d}' % float_field_len
    for i in range(n_conditions):
        column_line += column_field.format(i+1)
    return column_line
