# -*- coding: utf-8 -*-
"""
Thermochemistry
Vlachos group code to read, write, and generate thermochemical data.
Created on Fri Jul 7 12:40:00 2018
"""

import re
import inspect
from Thermochemistry import constants as c

class BaseThermo:
    """
    The Thermodynamic Parent class. Holds properties of a specie, the 
    statistical-mechanical thermodynamic model.

    Attributes
        name - str
            Name of the specie
        phase - str
            Phase of the specie
                G - gas
                S - surface
        elements - dict
            Composition of the species. Keys of dictionary are elements, 
            values are stoichiometric values in a formula unit
            e.g. CH3OH can be represented as:
            {
                'C': 1,
                'H': 4,
                'O': 1,
            }
        thermo_model - Thermochemistry.thermo_model class or custom class
            Class should have the following methods:
                get_CpoR
                get_HoRT
                get_SoR
                get_GoRT
        T_ref - float
            Reference temperature. Only used for reference species.
        HoRT_ref - float
            Reference dimensionless enthalpy corresponding to T_ref. Only used 
            for reference species.
    """

    def __init__(self, name, phase, elements, thermo_model, T_ref = None, HoRT_ref = None, **kwargs):
        self.name = name
        self.phase = phase
        self.elements = elements
        self.T_ref = T_ref
        self.HoRTref = HoRT_ref
        if inspect.isclass(thermo_model):
            self.thermo_model = _pass_expected_arguments(thermo_model, **kwargs)
        else:
            #If it's an object that has already been initialized
            self.thermo_model = thermo_model

def _get_expected_arguments(fn):
    """
    Returns the arguments expected by a function. Useful for determining
    where to assign **kwargs parameters.

    Parameters
        fn - Function or Class
            Function or class you would like to find the expected arguments.
    Returns
        tuple of str
            Expected arguments. If a class is specified, returns the
            expected arguments of __init__
    """
    #If class passed, use __init__ to find expected arguments

    if inspect.isclass(fn):
        fn = fn.__init__

    fn_code = fn.__code__
    arg_count = fn_code.co_argcount
    args = fn_code.co_varnames[:arg_count]
    return args

def _pass_expected_arguments(fn, **kwargs):
    """
    Finds expected values from a function or class and passes the
    appropriate arguments.

    Arguments
        fn - Function or Class
            Function or class you would like to find the expected arguments.
        **kwargs - Keyword arguments
            Keyword arguments that contain parameters to pass to fn
    """
    expected_args = _get_expected_arguments(fn)
    expected_arg_val = {}
    for arg in expected_args:
        try:
            expected_arg_val[arg] = kwargs[arg]
        except KeyError:
            continue
    return fn(**expected_arg_val)

def parse_formula(formula):
    """
    Parses chemical formula into its elements and returns it as a dictionary.

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

def get_molecular_weight(elements):
    """
    Molecular mass (in g/mol) given the elemental composition.
    Data taken from: https://en.wikipedia.org/wiki/Standard_atomic_weight

    Parameters
        elements - dict or str
            Elemental composition of species.
            If a dictionary is passed, the keys are the element symbol, atomic number, 
            or element name and the value is the stoichiometric coefficient.
            If a string is passed, the formula will be guessed using Thermochemistry.parse_formula

    Returns
        molecular_weight - float
            Molecular weight as float
    """
    if isinstance(elements, str):
        elements = parse_formula(elements)

    molecular_weight = 0.
    for element, coefficient in elements.items():
        molecular_weight += c.atomic_weight[element] * coefficient
    return molecular_weight