# -*- coding: utf-8 -*-
"""
Thermochemistry
Vlachos group code to read, write, and generate thermochemical data.
Created on Fri Jul 7 12:40:00 2018
"""

import re

def parse_formula(formula):
    """
    Parses chemical formula into its elments and returns it as a dictionary.

    Parameters    
        formula - string
            Chemical formula
            e.g. Al2O3
    Returns    
        elements - dict
            ELement composition of formula
            e.g. {'Al': 2, 'O': 3}
    """
    elements_tuples = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    elements = {}
    for (element, coefficient) in elements_tuples:
        if coefficient == '':
            elements[element] = 1
        else:
            elements[element] = int(coefficient)
    return elements