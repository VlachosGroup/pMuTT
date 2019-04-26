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
__version__ = '1.2.6'

import re
import inspect
import itertools
from warnings import warn
import numpy as np
from matplotlib import pyplot as plt
from pMuTT.io.json import remove_class
from pMuTT import constants as c

class _pMuTTBase:
    """Generic parent class to all pMuTT objects"""

    def __init__(self):
        pass
    
    def __eq__(self, other):
        try:
            other_dict = other.to_dict()
        except AttributeError:
            # If other doesn't have to_dict method, is not equal
            return False
        return self.to_dict() == other_dict

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = dict(self.__dict__)
        obj_dict['class'] = str(self.__class__)
        return obj_dict

    @classmethod
    def from_dict(cls, json_obj):
        """Recreate an object from the JSON representation.

        Parameters
        ----------
            json_obj : dict
                JSON representation
        Returns
        -------
            Obj : Appropriate object
        """
        json_obj = remove_class(json_obj)
        return cls(**json_obj)


class _ModelBase(_pMuTTBase):
    """Generic parent class to all model type objects"""

    def __init__(self):
        pass

    def get_Cv(units, **kwargs):
        """Calculate the heat capacity (constant V)

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            kwargs : keyword arguments
                Parameters needed by ``get_CvoR``
        Returns
        -------
            Cv : float
                Heat capacity (constant V) in appropriate units
        """
        self.get_CvoR(**kwargs)*c.R(units)

    def get_Cp(units, **kwargs):
        """Calculate the heat capacity (constant P)

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            kwargs : keyword arguments
                Parameters needed by ``get_CpoR``
        Returns
        -------
            Cp : float
                Heat capacity (constant P) in appropriate units
        """
        self.get_CpoR(**kwargs)*c.R(units)

    def get_U(units, T=c.T0('K'), **kwargs):
        """Calculate the internal energy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float, optional
                Temperature in K. Default is 298.15 K
            kwargs : keyword arguments
                Parameters needed by ``get_UoRT``
        Returns
        -------
            U : float
                Internal energy in appropriate units
        """
        self.get_UoRT(**kwargs)*T*c.R('{}/K'.format(units))

    def get_H(units, T=c.T0('K'), **kwargs):
        """Calculate the enthalpy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float, optional
                Temperature in K. Default is 298.15 K
            kwargs : keyword arguments
                Parameters needed by ``get_HoRT``
        Returns
        -------
            H : float
                Enthalpy in appropriate units
        """
        self.get_HoRT(**kwargs)*T*c.R('{}/K'.format(units))

    def get_S(units, **kwargs):
        """Calculate the entropy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            kwargs : keyword arguments
                Parameters needed by ``get_SoR``
        Returns
        -------
            S : float
                Entropy in appropriate units
        """
        self.get_SoR(**kwargs)*c.R(units)

    def get_F(units, T=c.T0('K'), **kwargs):
        """Calculate the Helmholtz energy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float, optional
                Temperature in K. Default is 298.15 K
            kwargs : keyword arguments
                Parameters needed by ``get_FoRT``
        Returns
        -------
            F : float
                Hemholtz energy in appropriate units
        """
        self.get_FoRT(**kwargs)*T*c.R('{}/K'.format(units))

    def get_G(units, T=c.T0('K'), **kwargs):
        """Calculate the Gibbs energy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float, optional
                Temperature in K. Default is 298.15 K
            kwargs : keyword arguments
                Parameters needed by ``get_GoRT``
        Returns
        -------
            G : float
                Gibbs energy in appropriate units
        """
        self.get_GoRT(**kwargs)*T*c.R('{}/K'.format(units))

def plot_1D(obj, x_name, x_values, methods, nrows=None, ncols=None,
            figure=None, ax=None, **kwargs):
    """Make a 1D plot

    Parameters
    ----------
        obj : Any model object (e.g. :class:`~pMuTT.statmech.trans.IdealTrans`),
              any species object (e.g. :class:`~pMuTT.statmech.StatMech`) or
              :class:`~pMuTT.reaction.Reaction`)
            pMuTT object to evaluate
        x_name : str
            Name of variable to vary
        x_values : iterable object
            x values to use
        methods : tuple of str or str
            Methods to use
        nrows : int, optional
            Number of rows of the subplot grid. It not specified, defaults
            to length of ``methods``
        ncols : int, optional
            Number of columns of the subplot grid. If not specified,
            defaults to 1
        figure : `matplotlib.figure.Figure`_, optional
            Add plot to this figure. If not specified, one will be
            generated
        ax : (N,) list of `matplotlib.axes.Axes.axis`_, optional
            Adds plot to this axis. If not specified, one will be generated
        kwargs : keyword arguments
            Other variables to use in the calculation. Method specific
            arguments can be passed by having a key that corresponds to
            the method name
            
            e.g. kwargs = {'get_H': {'units': 'kcal/mol'},
                            'get_S': {'units': 'cal/mol/K'}}
    Returns
    -------
        figure : `matplotlib.figure.Figure`_
            Figure
        ax : (N,) list of `matplotlib.axes.Axes.axis`_
            Axes of the plots where N is the length of ``methods``

    .. _`matplotlib.figure.Figure`: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html
    .. _`matplotlib.axes.Axes.axis`: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axis.html
    """
    # Check if single method passed
    if not _is_iterable(methods):
        methods = (methods,)

    # If rows/columns not specified, assing default values
    if nrows is None:
        nrows = len(methods)
        ncols = 1

    # Create the subplots
    if ax is None:
        figure, ax = plt.subplots(nrows=nrows, ncols=ncols)
    # Force ax to be a list
    if nrows*ncols == 1:
        ax = [ax]
    # Evaluate obj for each method
    for i, method in enumerate(methods):
        # If adding data to existing plot, skip methods with None
        if method is None:
            continue
        method_kwargs = _get_specie_kwargs(method, **kwargs)
        fn = getattr(obj, method)
        y = np.zeros_like(x_values)
        for j, x in enumerate(x_values):
            method_kwargs[x_name] = x
            y[j] = _force_pass_arguments(fn, **method_kwargs)
        ax[i].plot(x_values, y)
        ax[i].set_xlabel(x_name)
        ax[i].set_ylabel(method.replace('get_', ''))
    return (figure, ax)

def plot_2D(obj, x1_name, x1_values, x2_name, x2_values, methods,
            nrows=None, ncols=None, **kwargs):
    """Make a 2D plot

    Parameters
    ----------
        obj : Any model object (e.g. :class:`~pMuTT.statmech.trans.IdealTrans`),
              any species object (e.g. :class:`~pMuTT.statmech.StatMech`) or
              :class:`~pMuTT.reaction.Reaction`)
            pMuTT object to evaluate
        x1_name : str
            Name of first variable to vary
        x1_values : iterable object
            x1 values to use
        x2_name : str
            Name of second variable to vary
        x2_values : iterable object
            x2 values to use
        methods : tuple of str or str
            Methods to use
        nrows : int, optional
            Number of rows of the subplot grid. It not specified, defaults
            to length of ``methods``
        ncols : int, optional
            Number of columns of the subplot grid. If not specified,
            defaults to 1
        kwargs : keyword arguments
            Other variables to use in the calculation. Method specific
            arguments can be passed by having a key that corresponds to
            the method name
            
            e.g. kwargs = {'get_H': {'units': 'kcal/mol'},
                            'get_S': {'units': 'cal/mol/K'}}
    Returns
    -------
        figure : `matplotlib.figure.Figure`_
            Figure
        ax : `matplotlib.axes.Axes.axis`_
            Axes of the plots.

    .. _`matplotlib.figure.Figure`: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html
    .. _`matplotlib.axes.Axes.axis`: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axis.html
    """
    # Check if single method passed
    if not _is_iterable(methods):
        methods = (methods,)

    # If rows/columns not specified, assing default values
    if nrows is None:
        nrows = len(methods)
        ncols = 1

    x2_mesh, x1_mesh = np.meshgrid(x2_values, x1_values)

    figure, ax = plt.subplots(nrows=nrows, ncols=ncols)
    # Force ax to be a list
    if nrows*ncols == 1:
        ax = [ax]
    c = []
    cbar = []
    for i, method in enumerate(methods):
        method_kwargs = _get_specie_kwargs(method, **kwargs)
        fn = getattr(obj, method)
        y = np.zeros(shape=(len(x1_values), len(x2_values)))
        for j, x1 in enumerate(x1_values):
            method_kwargs[x1_name] = x1
            for k, x2 in enumerate(x2_values):
                method_kwargs[x2_name] = x2
                y[j, k] = fn(**method_kwargs)
        # Create colormap
        c.append(ax[i].pcolormesh(x1_mesh, x2_mesh, y))
        # Set colorbar
        cbar.append(figure.colorbar(c[i], ax=ax[i]))
        cbar[i].ax.set_title(method.replace('get_', ''))
        # Set axis labels
        ax[i].set_xlabel(x1_name)
        ax[i].set_ylabel(x2_name)
    return (figure, ax, c, cbar)

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
            passed by having a key for the species' name mapping onto a
            dictionary with the arguments to pass.

            e.g. For the reaction: H2 + 0.5O2 = H2O, hydrogen pressure of 1 atm,
            oxygen pressure of 0.5 atm can be specified by passing the
            following:
            kwargs = {
                'T': 298.,
                'H2_kwargs': {'P': 1.},
                'O2_kwargs': {P': 0.5},
            }
    Returns
    -------
        specie_kwargs : dict
            Dictionary containing the specie-specific kwargs
    """
    specie_kwargs = kwargs.copy()
    # Remove any keys related to other species
    for key in kwargs.keys():
        if 'kwargs' in key:
            temp_kwargs = specie_kwargs.pop(key, {})
            if key == '{}_kwargs'.format(specie_name):
                specie_specific_kwargs = temp_kwargs
    # See if there was an entry for the specific species
    try:
        specie_kwargs.update(specie_specific_kwargs)
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