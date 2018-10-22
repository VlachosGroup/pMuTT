# -*- coding: utf-8 -*-
"""
pMuTT.io_.chemkin

Reads reactions lists from Chemkin gas.inp and surf.inp files
"""
import re


def read_reactions(filename):
    """Directly read reactions from Chemkin gas.inp or surf.inp files

    Parameters
    ----------
        filename : str
            Input filename
    Returns
    -------
        Reactions : list of reactions
        Reactants : list of reactants found in reactions
        Products  : list of products found in reactions
    Raises
    ------
        FileNotFoundError
            If the file isn't found.

    """
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
    for Reacs, Prods in zip(LHS, RHS):
        Reactants.append(re.split(r' *\+ *| +', Reacs))
        Products.append(re.split(r' *\+ *| +', Prods)[0:-3])
    Reactions = []
    for rxn, Prods in zip(rxns, Products):
        Reactions.append(rxn[0:rxn.index(Prods[-1]) + len(Prods[-1])])
    return(Reactions, Reactants, Products)
