# -*- coding: utf-8 -*-
"""
PyMuTT.io_.excel

Read from/write to xlsx files of particular format.
"""

import numpy as np
import pandas as pd
import os
from ase.io import read
from PyMuTT import constants as c
from PyMuTT import parse_formula, get_molecular_weight
from PyMuTT.models.statmech import presets, StatMech
from PyMuTT.models.statmech.rot import get_geometry_from_atoms
from PyMuTT.models.statmech.rot import get_rot_temperatures_from_atoms

def read_excel(io, skiprows=[1], header=0, delimiter='.', **kwargs):
    """Reads an excel file and returns it as a list of dictionaries to
    initialize objects

    Parameters
    ----------
        io : str
            Name of the Excel spreadsheet
        skiprows : list, optional
            Rows to skip at the beginning (0-indexed).
            Default is [1] so comments can be put in that row
        header : int, optional
            Location to find header names (0-index). Default is 0
        delimiter : str, optional
            Delimiter to parse column names. Default is '.'
        **kwargs: keyword arguments
            Parameters used by `pandas.read_excel`_. Not required but some
            potentially useful parameters include:

            - sheet_name (str): Specify the name of the sheet you're reading
            - converters: Specify how to process certain columns
            - dtype (dict): Expected data type. Will be guessed if
              not specified
            - na_values (scalar, str, list-like, or dict): What strings to
              interpret as NaN
            - convert_float (bool): Converts integral floats to
              int (i.e. 1.0 --> 1)
    Returns
    -------
        excel_data : list of dict
            Can be used to initialize objects with the **kwargs syntax

    Notes
    -----
        Special rules exist for the following column headings

        - element
        - formula
        - atoms
        - statmech_model
        - vib_wavenumber
        - rot_temperatures
        - nasa.a_low
        - nasa.a_high

    .. _`pandas.read_excel`: https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_excel.html
    """
    input_data = pd.read_excel(io=io, skiprows=skiprows,
                               header=header, **kwargs)
    excel_path = os.path.dirname(io)
    thermos_out = []
    for row, row_data in input_data.iterrows():
        thermo_data = {}
        for col, cell_data in row_data.iteritems():
            # Special parsing instructions
            if pd.isnull(cell_data):
                continue
            elif 'element' in col:
                thermo_data = set_element(header=col, value=cell_data,
                                          output_structure=thermo_data,
                                          delimiter=delimiter)
            elif 'formula' in col:
                thermo_data = set_formula(formula=cell_data,
                                          output_structure=thermo_data)
            elif 'atoms' in col:
                thermo_data = set_atoms(path=cell_data, excel_path=excel_path,
                                        output_structure=thermo_data)
            elif 'statmech_model' in col:
                thermo_data = set_statmech_model(model=cell_data,
                                               output_structure=thermo_data)
            elif 'vib_wavenumber' in col:
                thermo_data = set_vib_wavenumbers(value=cell_data,
                                                 output_structure=thermo_data)
            elif 'rot_temperatures' in col:
                thermo_data = set_rot_temperatures(value=cell_data,
                                                   output_structure=thermo_data)
            elif 'nasa' in col:
                if 'a_low' in col:
                    thermo_data = set_nasa_a_low(header=col, value=cell_data,
                                                 output_structure=thermo_data)
                elif 'a_high' in col:
                    thermo_data = set_nasa_a_high(header=col, value=cell_data,
                                                  output_structure=thermo_data)
                else:
                    raise NotImplementedError('Does not support {}'
                                              .format(col))
            else:
                thermo_data[col] = cell_data
        thermos_out.append(thermo_data)
    return thermos_out


def set_element(header, value, output_structure, delimiter='.'):
    """Parses element header and assigns to output_structure['elements']

    Parameters
    ----------
        header : str
            String containing the element name. Element symbol should be at
            the end. e.g. 'element.O'
        value : int
            Amount found in formula
        output_structure : dict
            Structure to assign value. Will assign to
            output_structure['elements'][element]
        delimiter : str
            Delimiter for element. Element symbol should be at the end
    Returns
    -------
        output_structure : (dict)
            output_structure with new element added
    """
    element = header.split(delimiter)[-1]
    try:
        output_structure['elements'][element] = value
    except (NameError, KeyError):
        output_structure['elements'] = {element: value}
    return output_structure


def set_formula(formula, output_structure):
    """Parses stoichiometric formula unit and assigns to output_structure

    Parameters
    ----------
        formula : str
            Stoichiometric formula unit. e.g. H2O
            Note that an element cannot be specified multiple times.
            e.g. CH3OH is not supported
        output_structure : dict
            Structure to assign value. Will assign to
            output_structure['elements']
    Returns
    -------
        output_structure : dict
            output_structure with new elements added
    """
    elements = parse_formula(formula=formula)
    output_structure['elements'] = elements
    return output_structure


def set_atoms(path, output_structure, excel_path=None):
    """Reads the atoms object ans assigns to output_structure

    Parameters
    ----------
        path : str
            Location to import atoms object. If relative references used,
            the path should be relative to excel_path.
            See ase.read for supported formats
        excel_path : str
            Location where excel path is located
        output_structure : dict
            Structure to assign value. Will assign to output_structure['atoms']
    Returns
    -------
        output_structure: dict
            output_structure with atoms added
    """
    try:
        output_structure['atoms'] = read(path)
    except FileNotFoundError:
        try:
            output_structure['atoms'] = read(os.path.join(excel_path, path))
        except FileNotFoundError:
            print(path)
            raise FileNotFoundError('If using relative references for atoms '
                                    'files, use a path relative to the '
                                    'spreadsheet imported.', path)
    return output_structure


def set_statmech_model(model, output_structure):
    """Imports module and assigns the class to output_structure

    Parameters
    ----------
        model : str
            Thermodynamic model to import. Presets:

                - IdealGas
                - Harmonic
        output_structure : dict
            Structure to assign value. Will assign to
            output_structure['statmech_model']
    Returns
    -------
        output_structure: dict
            output_structure with new thermo model added
    """
    model = model.lower()
    output_structure['statmech_model'] = StatMech
    try:
        output_structure.update(presets[model])
    except KeyError:
        raise ValueError('Unsupported thermodynamic model, {}. See docstring '
                         'of presets in PyMuTT.models.statmech for supported '
                         'models.'.format(model))
    return output_structure


def set_vib_wavenumbers(value, output_structure):
    """Parses element header and assigns to output_structure['vib_wavenumber']

    Parameters
    ----------
        value : float
            Vibrational frequency in 1/cm
        output_structure : dict
            Structure to assign value. Will assign to
            output_structure['elements'][element]
    Returns
    -------
        output_structure: dict
            output_structure with new vibration added
    """
    try:
        output_structure['vib_wavenumbers'].append(value)
    except (NameError, KeyError):
        output_structure['vib_wavenumbers'] = [value]
    return output_structure


def set_rot_temperatures(value, output_structure):
    """Parses element header and assigns to output_structure['rot_temperatures']

    Parameters
    ----------
        value : float
            Vibrational frequency in 1/cm
        output_structure : dict
            Structure to assign value. Will assign to
            output_structure['elements'][element]
    Returns
    -------
        output_structure: dict
            output_structure with new vibration added
    """
    try:
        output_structure['rot_temperatures'].append(value)
    except (NameError, KeyError):
        output_structure['rot_temperatures'] = [value]
    return output_structure

def set_nasa_a_low(header, value, output_structure, delimiter='.'):
    """Parses a_low parameter for ``PyMuTT.models.empirical.nasa.Nasa`` object

    Parameters
    ----------
        header : str
            Name of the header. Used to determine coefficient. 
            Assumes zero index and header takes the format:

            nasa[delimiter]a_low[delimiter][index]
            e.g. nasa.a_low.0
        value : float
            a_low value
        output_structure : dict
            Structure to assign value. Will assign to output_structure['a_low']
        delimiter : str
            How to parse header to find the coefficient
    Returns
    -------
        output_structure : dict
            output_structure with a_low value added
    """
    i = int(header.split(delimiter)[-1])
    try:
        output_structure['a_low'][i] = value
    except KeyError:
        output_structure['a_low'] = np.zeros(7,)
        output_structure['a_low'][i] = value
    return output_structure


def set_nasa_a_high(header, value, output_structure, delimiter='.'):
    """Parses a_high parameter for ``PyMuTT.models.empirical.nasa.Nasa`` object

    Parameters
    ----------
        header : str
            Name of the header. Used to determine coefficient. 
            Assumes zero index and header takes the format:

            nasa[delimiter]a_high[delimiter][index]
            e.g. nasa.a_high.0
        value : float
            a_high value
        output_structure : dict
            Structure to assign value. Will assign to
            output_structure['a_high']
        delimiter : str
            How to parse header to find the coefficient
    Returns
    -------
        output_structure : dict
            output_structure with a_high value added
    """
    i = int(header.split(delimiter)[-1])
    try:
        output_structure['a_high'][i] = value
    except KeyError:
        output_structure['a_high'] = np.zeros(7,)
        output_structure['a_high'][i] = value
    return output_structure
