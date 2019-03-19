# -*- coding: utf-8 -*-
"""
pMuTT.io.vasp
Description: methods to read from vasp outcar files
"""

import os
import re


def get_vib_wavenumber_from_line(in_line):
    """Parses in_line for real frequencies

    Parameters
    ----------
    in_line: str
        Line containing frequency in OUTCAR
    Returns
    -------
    vib_wavenumber: float
        Vibrational wavenumber in cm-1
    """

    pattern = re.compile(r'(\d+\.?\d+) cm-1')  # pattern for frequency in cm-1
    m = pattern.search(in_line)
    try:
        vib_wavenumber = float(m[1])
        return vib_wavenumber
    except TypeError as e:
        print('No frequency in this line')
        raise e


def set_vib_wavenumbers_from_outcar(in_file, output_structure,
                                    min_frequency_cutoff=0.,
                                    include_imaginary=False):
    """Parses OUTCAR files for vibrational frequencies and assigns to
    output_structure['vib_wavenumber']. Imaginary frequencies are represented
    by negative numbers.

    Parameters
    ----------
        in_file: str
            OUTCAR file of frequency jobs
        output_structure: dict
            Structure to assign value. Will assign to
            output_structure['elements'][element]
        min_frequency_cutoff: float, optional
            Only frequencies less than min_frequency_cutoff (in cm-1)
            are read from OUTCAR. Default is 0 cm-1
        include_imaginary: bool, optional
            Whether imaginary frequencies should be included. Default is False.
    Returns
    -------
        output_structure: dict
            Output_structure with new vibration added
    Raises
    ------
        FileNotFoundError
            Raised if in_file does not exist
    """
    if not os.path.isfile(in_file):
        raise FileNotFoundError('Invalid OUTCAR filename: {}'.format(in_file))

    vib_wavenumbers = list()
    real_pattern = re.compile(r'f[ ]*=')  # pattern for real frequencies
    if include_imaginary:
        # pattern for imaginary frequencies
        imag_pattern = re.compile(r'f/i[ ]*=')

    with open(in_file, "r") as out_fp:
        outcar_contents = out_fp.readlines()
    for line in outcar_contents:
        if real_pattern.search(line) is not None:
            try:
                vib_wavenumber = get_vib_wavenumber_from_line(line)
            except TypeError:
                pass  # if no frequency in line continue
            else:
                if vib_wavenumber > min_frequency_cutoff:
                    vib_wavenumbers.append(vib_wavenumber)
        elif include_imaginary and imag_pattern.search(line) is not None:
            try:
                vib_wavenumber = get_vib_wavenumber_from_line(line)
            except TypeError:
                pass  # if no frequency in line continue
            else:
                vib_wavenumbers.append(-vib_wavenumber)

    if len(vib_wavenumbers) == 0:
        print('no frequencies found in file')
    output_structure['vib_wavenumbers'] = vib_wavenumbers
    return output_structure
