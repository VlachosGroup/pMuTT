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
__version__ = '1.2.3'

import re
import inspect
from warnings import warn
import numpy as np
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


def _is_iterable(val):
    """
    Checks if the input if an iterable. This function will return False if a
    string is passed due to its use in pMuTT.

    Parameters
    ----------
        val : iterable or non-iterable
            Value to check if iterable
    Returns
    -------
        is_iterable : bool
            True if iterable. False if not iterable or string.
    """
    if isinstance(val, str):
        return False
    else:
        # If it's not a string, check if it's iterable
        try:
            iter(val)
        except TypeError:
            return False
        else:
            return True


def _get_mode_quantity(mode, method_name, raise_error=True, raise_warning=True,
                       default_value=0., **kwargs):
    """Calculate the quantity from that mode.

    Parameters
    ----------
        mode : ``pMuTT.statmech`` object
            Trans, Vib, Rot, Elec, or Nucl StatMech model
        method_name : str
            Name of method to use to calculate quantity.
        raise_error : bool, optional
            If True, raises an error if any of the modes do not have the
            quantity of interest. Default is True
        raise_warning : bool, optional
            Only relevant if raise_error is False. Raises a warning if any
            of the modes do not have the quantity of interest. Default is
            True
        default_value : float, optional
            Default value if the object does not contain the method. Default is
            0
        kwargs : key-word arguments
            Parameters passed to each mode
    Returns
    -------
        quantity : float
            Quantity of the mode.
    Raises
    ------
        AttributeError
            If raise_error is True and the mode does not have the method_name
    """
    try:
        method = getattr(mode, method_name)
    except AttributeError as e:
        if raise_error:
            raise e
        elif raise_warning:
            warn(e, RuntimeWarning)
        quantity = default_value
    else:
        quantity = _pass_expected_arguments(method, **kwargs)
    return quantity


def _get_specie_kwargs(specie_name, **kwargs):
    """Gets the keyword arguments specific to a specie

    Parameters
    ----------
        specie_name : str
            Name of the specie
        kwargs : keyword arguments
            Parameters with the conditions. Specie specific parameters can be
            passed by having a key named 'specie' mapping onto a dictionary
            whose keys are the species names.

            e.g. For the reaction: H2 + 0.5O2 = H2O
            kwargs = {
                'T': 298.,
                'specie': {
                    'H2': {
                        'P': 1.,
                    },
                    'O2': {
                        'P': 0.5,
                    },
                }
            }
    Returns
    -------
        specie_kwargs : dict
            Dictionary containing the specie-specific kwargs
    """
    specie_kwargs = kwargs.copy()
    specie_specific = specie_kwargs.pop('specie', None)
    # See if there was an entry for the specific species
    try:
        specie_kwargs.update(specie_specific[specie_name])
    except (KeyError, TypeError, NameError):
        pass
    return specie_kwargs


def _apply_numpy_operation(quantity, operation, verbose=False):
    """Apply operation to quantity

    Parameters
    ----------
        quantity : (N,) np.ndarray
            Array with the quantity of interest
        operation : str
            Numpy operation to perform
        verbose : bool, optional
            If True, returns quantity with no further operations. Default is
            False
    Returns
    -------
        quantity : float or (N,) np.ndarray
            Quantity of interest in the desired format
    """
    if verbose:
        out_quantity = quantity
    else:
        np_method = getattr(np, operation)
        out_quantity = np_method(quantity)
    return out_quantity


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


def pMuTT_list_to_dict(pMuTT_list, key='name'):
    """Converts a pMuTT list to a dictionary using a specified attribute. This
    allows for quicker searching.

    Parameters
    ----------
        pMuTT_list : list of objects
            List of pMuTT objects to convert
        key : str
            Name of attribute used as the keys for the dictionary
    Returns
    -------
        pMuTT_dict : dict
            Dictionary of pMuTT objects
    """
    return {getattr(obj, key): obj for obj in pMuTT_list}

def format_conditions(**kwargs):
    """Converts an arbitrary number of lists to a list of dictionaries. Useful
    for specifying the conditions in pMuTT.io.chemkin.write_EA

    Parameters
    ----------
        kwargs - keyword arguments
            Lists of the conditions where each index corresponds to a run
    Returns
    -------
        conditions - list of dict
            A list where each element is a dictionary containing the conditions
            for a specific run
    """
    conditions = []
    for cond_name, cond_values in kwargs.items():
        for i, cond_value in enumerate(cond_values):
            try:
                conditions[i][cond_name] = cond_value
            except (IndexError, KeyError):
                conditions.append({cond_name: cond_value})
    return conditions