# -*- coding: utf-8 -*-
"""
pMuTT.io_.chemkin

Reads reactions lists from Chemkin gas.inp and surf.inp files
"""
import re
from pMuTT.models import pMuTT_list_to_dict


def read_reactions(filename, species):
    """Directly read reactions from Chemkin gas.inp or surf.inp files

    Parameters
    ----------
        filename : str
            Input filename for Chemkin surf or gas .inp file
        species : obj
            List of NASA object containing thermodynamic properties for
            all Reactants and Products in Reactions
    Returns
    -------
        Reactions   : list of reactions
        Reactants   : list of reactants found in reactions
        React_obj   : list of NASA polynomials for each Reactant
        React_stoic : list of reaction stiociometries for Reactants
        Products    : list of products found in reactions
        Prod_obj    : list of NASA polynomials for each Product
        Prod_stoic  : list of reaction stiociometries for Products
    Raises
    ------
        FileNotFoundError
            If the surf.inp or gas.inp file isn't found.
        NameError
            If the species file does not exist
        AttributeError
            If the species list is incorrect format
    """
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
        for RR in Reactants[-1]:
            stoic = re.findall(r'^[0-9]*', RR)[0]
            if stoic == '':
                stoic = 1
            else:
                RR = RR.replace(stoic, "")
                stoic = int(stoic)
            RS.append(stoic)
            R.append(species_dict[RR])
        React_stoic.append(RS)
        P = []
        PS = []
        for PP in Products[-1]:
            stoic = re.findall(r'^[0-9]*', PP)[0]
            if stoic == '':
                stoic = 1
            else:
                PP = PP.replace(stoic, "")
                stoic = int(stoic)
            PS.append(stoic)
            P.append(species_dict[PP])
        Prod_stoic.append(PS)
        React_obj.append(R)
        Prod_obj.append(P)
    Reactions = []
    for rxn, Prods in zip(rxns, Products):
        Reactions.append(rxn[0:rxn.index(Prods[-1]) + len(Prods[-1])])
    return(Reactions, Reactants, React_obj, React_stoic,
           Products, Prod_obj, Prod_stoic)
