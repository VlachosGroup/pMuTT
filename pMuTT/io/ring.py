from pMuTT.reaction import Reaction, Reactions

def read_reactions(filename, species, species_delimiter='.',
                   reaction_delimiter='>>', raise_error=True,
                   raise_warning=True):
    """Reads reactions from RING output file.

    Parameters
    ----------
        filename : str
            Input filename
        species : dict
            Dictionary using the names as keys. If you have a list of
            species, use pMuTT.pMuTT_list_to_dict to make a dict.
        species_delimiter : str, optional
            Delimiter that separate species. Leading and trailing spaces
            will be trimmed. Default is '.'
        reaction_delimiter : str, optional
            Delimiter that separate states of the reaction. Leading and
            trailing spaces will be trimmed. Default is '>>'
        raise_error : bool, optional
            If True, raises an error if the transition state is not located in
            species. Default is True
        raise_warning : bool, optional
            Only relevant if raise_error is False. Raises a warning if the
            transition state is not located in species. Default is True
    Returns
    -------
        reactions : list of :class:`~pMuTT.reaction.Reaction` objects
            Reactions
    """
    rxns = []
    with open(filename, 'r') as f_ptr:
        for line in f_ptr:
            # Skip lines that do not have a reaction
            if reaction_delimiter not in line:
                continue
            reaction_str = line.replace('\n', '')
            rxn = Reaction.from_string(reaction_str=reaction_str,
                                       species=species,
                                       species_delimiter=species_delimiter,
                                       reaction_delimiter=reaction_delimiter,
                                       raise_error=raise_error,
                                       raise_warning=raise_warning)
            rxns.append(rxn)
    return Reactions(reactions=rxns)