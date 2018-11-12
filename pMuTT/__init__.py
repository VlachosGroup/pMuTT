# -*- coding: utf-8 -*-
"""
pMuTT
"""

####
#
# setuptools likes to see a name for the package,
# and it's best-practices to have the __version__
# present, too:
#
name = 'pMuTT'
__version__ = '1.1.1'
#
####

import re
import inspect
from warnings import warn
from pMuTT import constants as c


def _get_expected_arguments(fn):
    """Returns the arguments expected by a function. Useful for determining
    where to assign **kwargs parameters.

    Parameters
    ----------
        fn : function or class
            Function or class you would like to find the expected arguments.
    Returns
    -------
        expected_arguments : tuple of str
            Expected arguments. If a class is specified, returns the
            expected arguments of __init__
    """

    # If class passed, use __init__ to find expected arguments
    if inspect.isclass(fn):
        fn = fn.__init__

    fn_code = fn.__code__
    arg_count = fn_code.co_argcount
    args = fn_code.co_varnames[:arg_count]
    return args


def _pass_expected_arguments(fn, **kwargs):
    """Finds expected values from a function or class and passes the
    appropriate arguments.

    Parameters
    ----------
        fn : Function or class
            Function or class you would like to find the expected arguments.
        verbose : bool, Optional
            If True, warns when an argument could not be found. Default is True
        **kwargs :
            Keyword arguments that contain parameters to pass to fn
    Returns
    -------
        fn_or_class_output :
        Output of fn that has been fed the expected arguments.
    """
    expected_args = _get_expected_arguments(fn)
    expected_arg_val = {}
    for arg in expected_args:
        if arg == 'self':
            continue

        try:
            expected_arg_val[arg] = kwargs[arg]
        except KeyError:
            continue
    return fn(**expected_arg_val)


def _kwargs_allowed(fn):
    """Checks to see if kwargs is allowed

    Parameters
    ----------
        fn : Function or class
            Function or class you would like to check if kwargs is allowed.
    Returns
    -------
        kwargs_allowed : bool
            True if kwargs are allowed. False otherwise.
    """
    sig = inspect.signature(fn)
    for param in sig.parameters.values():
        if param.kind == param.VAR_KEYWORD:
            return True
    else:
        return False

def _force_pass_arguments(fn, **kwargs):
    """Checks to see if fn accepts kwargs. If it does, pass arguments using
    kwargs. If not, pass arguments using docstring

    Parameters
    ----------
        fn : Function or class
            Function or class you would like to pass the arguments.
        verbose : bool, Optional
            If True, warns when an argument could not be found. Default is True
        **kwargs :
            Keyword arguments that contain parameters to pass to fn
    Returns
    -------
        fn_or_class_output :
        Output of fn that has been fed the expected arguments.
    """
    if _kwargs_allowed(fn):
        return fn(**kwargs)
    else:
        return _pass_expected_arguments(fn, **kwargs)


def parse_formula(formula):
    """Parses chemical formula into its elements and returns it as a
    dictionary.

    Parameters
    ----------
        formula : str
            Chemical formula e.g. Al2O3 or CH3CH2CH3
    Returns
    -------
        elements : dict
            Element composition of formula e.g. {'Al': 2, 'O': 3}
            Element composition of formula e.g. {'C': 3, 'H': 8}
    """
    elements_tuples = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    elements = {}
    for (element, coefficient) in elements_tuples:
        elements[element] = elements.get(element, 0) + int(coefficient or '1')
    return elements


def get_molecular_weight(elements):
    """Molecular mass (in g/mol) given the elemental composition.
    Data taken from: https://en.wikipedia.org/wiki/Standard_atomic_weight

    Parameters
    ----------
        elements : dict or str
            Elemental composition of species.

            If a dictionary is passed, the keys are the element symbol, atomic
            number, or element name and the value is the stoichiometric
            coefficient.
            If a string is passed, the formula will be guessed using
            pMuTT.parse_formula

    Returns
    -------
        molecular_weight : float
            Molecular weight as float in kg/mol
    """
    if isinstance(elements, str):
        elements = parse_formula(elements)

    molecular_weight = 0.
    for element, coefficient in elements.items():
        molecular_weight += c.atomic_weight[element] * coefficient

    return molecular_weight
