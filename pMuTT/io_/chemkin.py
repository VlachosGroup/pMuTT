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

def write_EAs(reactions, conditions, filename='EAs.inp', 
              method_name='get_EoRT_act', float_format=' .2E',
              species_delimiter='+', reaction_delimiter='<=>',
              stoich_format='.0f', newline='\n', column_delimiter='  '):
    """Writes the EAs.inp file for Chemkin

    Parameters
    ---------- 
        reactions : :class:`~pMuTT.reaction.Reactions` object
            Reactions to write
        conditions : list of dicts
            Conditions to evaluate each reaction. The key of the dictionaries
            should be a relevant quantity to evaluate the reaction (e.g. T, P)
        filename : str, optional
            Filename for the EAs.inp file. Default is 'EAs.inp'
        method_name : str, optional
            Method to use to calculate values. Typical values are:

            - 'get_EoRT_act' (default)
            - 'get_delta_HoRT'
            - 'get_delta_GoRT'

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
    n_reactions = len(reactions)
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
        ('!the number of times the species appears in the reaction. ',
         'If not using'),
        '!MultiInput, then only the first value is used.',
        '  {}  !Number of reactions'.format(n_reactions)]

    # Find length to pad reactions
    max_rxn_len = _get_max_reaction_len(reactions=reactions,
                                        species_delimiter=species_delimiter,
                                        reaction_delimiter=reaction_delimiter,
                                        stoich_format=stoich_format,
                                        include_TS=False)
    rxn_padding = max_rxn_len + len(column_delimiter)

    # Define string formats to use
    float_field = '{:%s}' % float_format
    str_field = '{:%d}' % rxn_padding

    # Add comment line showing column numbers
    column_line = '!{}'.format(' '*(rxn_padding-1))
    float_field_len = len(float_field.format(np.pi))+len(column_delimiter)
    column_field = '{:> %d}' % float_field_len
    for i in range(len(conditions)):
        column_line += column_field.format(i+1)
    lines.append(column_line)

    # Add line for each reaction step
    for reaction in reactions:
        line = [str_field.format(
                    reaction.to_string(species_delimiter=species_delimiter,
                                       reaction_delimiter=reaction_delimiter,
                                       stoich_format=stoich_format,
                                       include_TS=False))]
        for condition in conditions:
            method = getattr(reaction, method_name)

            # If the transition state does not exist
            if reaction.transition_state is None:
                # If function is get_delta_HoRT or get_delta_GoRT, and the 
                # transition state does not exist, specify this is not a 
                # transition state calculation
                if 'delta' in method_name:
                    condition['activation'] = False
                quantity = _force_pass_arguments(method, **condition)
                # If using enthalpies, activation energies have to be positive
                if method_name == 'get_delta_HoRT':
                    quantity = np.max([0., quantity])
            else:
                # If function is get_delta_HoRT or get_delta_GoRT, specify this
                # is a transition state calculation
                if 'delta' in method_name:
                    condition['activation'] = True
                quantity = _force_pass_arguments(method, **condition)

            line.append(float_field.format(quantity))
        lines.append(column_delimiter.join(line))
    lines.append('EOF')
    with open(filename, 'w', newline=newline) as f_ptr:
        f_ptr.write('\n'.join(lines))

def write_gas(species, filename='gas.inp', T=c.T0('K'), reactions=[],
              species_delimiter=' + ', reaction_delimiter=' <=> ',
              act_method_name='get_E_act', act_unit='kcal/mol', 
              float_format=' .3E', stoich_format='.0f', newline='\n',
              column_delimiter='  ', **kwargs):
    """Writes the gas.inp Chemkin file.

    Parameters
    ----------
        species : list of :class:`~pMuTT.empirical.nasa.Nasa objects
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
    for specie in species:
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
            column_delimiter=column_delimiter, T=T, **kwargs)
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

def write_surf(species, filename='surf.inp', T=c.T0('K'), reactions=[],
               species_delimiter=' + ', reaction_delimiter=' <=> ',
               act_method_name='get_E_act', act_unit='kcal/mol', 
               float_format=' .3E', stoich_format='.0f', newline='\n',
               column_delimiter='  ', use_mw_correction=True, **kwargs):
    # Organize species by their catalyst sites
    cat_adsorbates = {}
    unique_cat_sites = []
    for specie in species:
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
                column_delimiter=column_delimiter, T=T, **kwargs)

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
                   'kcal/mol'),
                  ('!Adsorption reactions can be represented using the STICK '
                   'keyword'),
                  'REACTIONS{2}{0}{2}{1}'.format(mw_str, act_unit_str, 
                                                 column_delimiter)])
    lines.extend(reaction_lines)
    lines.append('END')

    # Write the file
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

def _write_reaction_lines(reactions, species_delimiter, reaction_delimiter, 
                          include_TS, stoich_format, act_method_name, act_unit, 
                          float_format, column_delimiter,
                          **kwargs):
    """Write the reaction lines in the Chemkin format

    Parameters
    ----------

    Returns
    -------
        reaction_lines : str
            Reactions represented in Chemkin format
    """
    max_reaction_len = _get_max_reaction_len(
            reactions=reactions, species_delimiter=species_delimiter,
            reaction_delimiter=reaction_delimiter, stoich_format=stoich_format,
            include_TS=include_TS)
    reaction_field = '{:<%s}' % max_reaction_len
    float_field = '{:%s}' % float_format

    reaction_lines = []
    for reaction in reactions:
        # Get reaction string
        reaction_str = reaction_field.format(
                reaction.to_string(species_delimiter=species_delimiter,
                                   reaction_delimiter=reaction_delimiter,
                                   stoich_format=stoich_format,
                                   include_TS=include_TS))
        # Calculate preexponential factor
        if reaction.is_adsorption:
            A = reaction.sticking_coeff
        else:
            # If using delta_G, take out entropic contribution in A
            if 'get_delta_G' in act_method_name:
                include_entropy = False
            else:
                include_entropy = True
            A = reaction.get_A(include_entropy=include_entropy, **kwargs)
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