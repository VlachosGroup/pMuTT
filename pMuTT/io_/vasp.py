# -*- coding: utf-8 -*-
"""
pMuTT.io_.vasp
Description: methods to read from vasp outcar files
"""

import os
import re


def take_vib_wavenumber_from_line(in_line):
    """Parses in_line for real frequencies

    Parameters
    ----------
    in_line: str
            line containing frequency in OUTCAR
    Returns
    -------
    vib_wavenumber: float
            vibrational wavenumber in cm-1
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
                                    min_frequency_cutoff):
    """Parses OUTCAR files and assigns to output_structure['vib_wavenumber']

    Parameters
    ----------
        in_file: str
            OUTCAR file of frequency jobs
        output_structure: dict
            Structure to assign value. Will assign to
            output_structure['elements'][element]
        min_frequency_cutoff: float
            Frequencies > min_frequency_cutoff (cm-1)
            read from OUTCAR
    Returns
    -------
        output_structure: dict
            output_structure with new vibration added

    """
    if not os.path.isfile(in_file):
        raise FileNotFoundError('invalid outcar filename: {}'.format(in_file))
    real_vib_wavenumbers = list()
    real_pattern = re.compile(r'f[ ]*=')  # pattern for real frequencies
    try:
        with open(in_file, "r") as out_fp:
            outcar_contents = out_fp.readlines()
    except Exception as e:
        # Unknown error opening file
        raise e
    for line in outcar_contents:
        if real_pattern.search(line):
            try:
                vib_wavenumber = take_vib_wavenumber_from_line(line)
            except TypeError:
                pass  # if no frequency in line continue
            else:
                if vib_wavenumber > min_frequency_cutoff:
                    real_vib_wavenumbers.append(vib_wavenumber)
    if not real_vib_wavenumbers:
        print('no frequencies found in file')
    output_structure['vib_wavenumbers'] = real_vib_wavenumbers
    return output_structure
