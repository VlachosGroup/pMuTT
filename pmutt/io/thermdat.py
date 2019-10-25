# -*- coding: utf-8 -*-
"""
pmutt.io.thermdat

Read from/write to thermdat files.
"""

from datetime import datetime

import numpy as np

from pmutt import pmutt_list_to_dict
from pmutt.empirical.nasa import Nasa


def read_thermdat(filename, format='list', key='name'):
    """Directly read thermdat file that is in the Chemkin format

    Parameters
    ----------
        filename : str
            Input filename
        format : str, optional
            Format to output NASA polynomials. Supported options are:
            'list', 'tuple', 'dict'. Default is 'list'
        key : str, optional
            If `format` is 'dict', uses this attribute as the key for the
            output dictionary. Default is 'name'
    Returns
    -------
        Nasas : list, tuple or dict of :class:`~pmutt.empirical.nasa.Nasa`
    Raises
    ------
        FileNotFoundError
            If the file isn't found.
        IOError
            Invalid line number found.
    """
    
    species = []
    with open(filename, 'r') as f_ptr:
        for line in f_ptr:
            '''
            Lines to skip
            '''
            # Skip the header line
            if 'THERMO' in line:
                continue
            # Skip the end line
            if 'END' in line:
                continue
            # Skip blank lines
            if line == '\n':
                continue
            # Skip comment lines
            if line[0] == '!':
                continue
            # Skip header temperatures
            if _is_temperature_header(line):
                continue

            '''
            Parse lines
            '''
            line_num = _read_line_num(line)
            if line_num == 1:
                nasa_data = _read_line1(line)
            elif line_num == 2:
                nasa_data = _read_line2(line, nasa_data)
            elif line_num == 3:
                nasa_data = _read_line3(line, nasa_data)
            elif line_num == 4:
                nasa_data = _read_line4(line, nasa_data)
                species.append(Nasa(**nasa_data))
            else:
                err_msg = ('Invalid line number, {}, in thermdat file: {}'
                           ''.format(line_num, filename))
                raise IOError(err_msg)
    # Format the NASA polynomials in the required format
    if format == 'list':
        pass
    elif format == 'tuple':
        species = tuple(species)
    elif format == 'dict':
        species = pmutt_list_to_dict(species, key=key)
    else:
        err_msg = ('Unsupported format: {}. See pmutt.io.thermdat.read_thermdat'
                   ' docstring for supported formats.'.format(format))
        raise ValueError(err_msg)
    return species


def _get_fields(line, delimiter=' ', remove_fields=['', '\n']):
    """Gets the fields from a line delimited by delimiter and without entries
    in remove_fields

    Parameters
    ----------
        line : str
            Line to find the fields
        delimiter : str
            Text separating fields in string
        remove_fields : list of str
            Fields to delete
    Returns
    -------
        fields : list of str
            Fields of the line
    """
    for remove_field in remove_fields:
        line = line.replace(remove_field, '')
    all_fields = line.split(delimiter)
    fields = []
    for field in all_fields:
        if all([field != remove_field for remove_field in remove_fields]):
            # Add field if it does not match any of the remove_fields
            fields.append(field)
    return fields


def _is_temperature_header(line):
    """Determines if the line if the temperature header by seeing if the line
    only contains three numbers.

    Parameters
    ----------
        line : str
            Line to test
    Returns
    -------
        temperature_header : bool
            True if the line is the temperature header. False otherwise.
    """

    fields = _get_fields(line)
    n_num = 0
    for field in fields:
        # See if the field is a float
        try:
            float(field)
        except ValueError:
            # Contains text field
            return False
        else:
            n_num += 1
    # Temperature header contains 3 floats
    if n_num == 3:
        return True
    else:
        return False


def _read_line_num(line):
    """Reads the line number. Assumes the line number is the last character

    Parameters
    ----------
        line : str
            Line to be read
    Returns
    -------
        line_num : int
            Line number.
    """

    fields = _get_fields(line)
    return int(fields[-1])


def _read_line1(line):
    """Reads the first line of a thermdat specie

    Parameters
    ----------
        line : str
            Line 1 of thermdat specie
    Returns
    -------
        nasa_data : dict
            Nasa input fields
    """
    nasa_data = {}
    ref_pos = 24
    ref_offset = 5
    max_elements = 4
    phase_pos = 44

    # Store the name
    blank_pos = line.find(' ')
    nasa_data['name'] = line[:blank_pos]

    # Store the notes if any
    notes = line[blank_pos:ref_pos].strip()
    if len(notes) > 0:
        nasa_data['notes'] = notes

    # Store the elements
    nasa_data['elements'] = {}
    for i in range(max_elements):
        blank_pos = line.find(' ', ref_pos)
        # All the elements have been assigned
        if blank_pos == ref_pos:
            break

        element = line[ref_pos:blank_pos]
        ref_pos += ref_offset
        coeff = int(line[blank_pos:ref_pos])

        nasa_data['elements'][element] = coeff

    # Store the phase
    nasa_data['phase'] = line[phase_pos]

    # Store the temperatures
    fields = _get_fields(line[phase_pos+1:])
    nasa_data['T_low'] = float(fields[0])
    nasa_data['T_high'] = float(fields[1])
    nasa_data['T_mid'] = float(fields[2])
    return nasa_data


def _read_line2(line, nasa_data):
    """Reads the second line of a thermdat specie

    Parameters
    ----------
        line : str
            Line 2 of thermdat specie
        nasa_data : dict
            Pre-filled Nasa input fields
    Returns
        nasa_data : dict
            Nasa input fields
    """
    # Locations to find a values
    positions = [0, 15, 30, 45, 60]
    offset = 15

    nasa_data['a_high'] = np.zeros(7)

    for i, position in enumerate(positions):
        nasa_data['a_high'][i] = float(line[position:position+offset])
    return nasa_data


def _read_line3(line, nasa_data):
    """Reads the third line of a thermdat specie

    Parameters
    ----------
        line : str
            Line 3 of thermdat specie
        nasa_data : dict
            Pre-filled Nasa input fields
    Returns
    -------
        nasa_data : dict
            Nasa input fields
    """
    # Locations to find a values
    positions = [0, 15, 30, 45, 60]
    offset = 15

    nasa_data['a_low'] = np.zeros(7)

    j = 5  # Counter for a_high
    k = 0  # Counter for a_low
    for i, position in enumerate(positions):
        if i < 2:
            nasa_data['a_high'][j] = float(line[position:position+offset])
            j += 1
        else:
            nasa_data['a_low'][k] = float(line[position:position+offset])
            k += 1
    return nasa_data


def _read_line4(line, nasa_data):
    """Reads the third line of a thermdat specie

    Parameters
    ----------
        line : str
            Line 3 of thermdat specie
        nasa_data : dict
            Pre-filled Nasa input fields
    Returns
    -------
        nasa_data : dict
            Nasa input fields
    """
    # Locations to find a values
    positions = [0, 15, 30, 45]
    offset = 15

    j = 3
    for position in positions:
        nasa_data['a_low'][j] = float(line[position:position+offset])
        j += 1
    return nasa_data


def write_thermdat(nasa_species, filename=None, write_date=True, supp_data=None,
                   supp_txt=None, newline='\n'):
    """Writes thermdats in the Chemkin format

    Parameters
    ----------
        nasa_species : list or dict of :class:`~pmutt.empirical.nasa.Nasa`
            List of species to populate thermdat
        filename : str, optional
            Output file name. If not specified, returns thermdat as str
        supp_data : str, optional
            Additional thermdat entries to include. Must be in therndat format.
        supp_txt : str, optional
            Comment field to preceed nasa_species entries. Each line needs to
            begin with a ! so it is recognized as a comment.
        write_date : bool, optional
            Whether or not the date should be written. If False, writes the
            first 8 characters of ``notes`` attribute. Defaults to True
        newline : str, optional
            Newline character to use. Default is the Unix convention (\\n)
    Returns
    -------
        lines_out : str
            Thermdat lines as a string if ``filename`` is None
    """
    lines = []
    # Add header
    lines.append('THERMO ALL\n       100       500      1500\n')
    # Add supplementary data
    if supp_data is not None:
        if supp_data[-1] != '\n':
            supp_data += '\n'
        lines.append(supp_data)
    # Add supplementary text
    if supp_txt is not None:
        if supp_txt[-1] != '\n':
            supp_data += '\n'
        lines.append(supp_data)

    # Iterate over nasa_species using appropriate method
    if isinstance(nasa_species, dict):
        nasa_iter = nasa_species.values()
    else:
        nasa_iter = iter(nasa_species)

    for nasa_specie in nasa_iter:
        lines.append(_write_line1(nasa_specie, write_date))
        lines.append(_write_line2(nasa_specie))
        lines.append(_write_line3(nasa_specie))
        lines.append(_write_line4(nasa_specie))
    lines.append('END')

    # Write lines or return the string to the user
    lines_out = ''.join(lines)
    if filename is not None:
        with open(filename, 'w', newline=newline) as f_ptr:
            f_ptr.write(lines_out)
    else:
        return lines_out


def _write_line1(nasa_specie, write_date=True):
    """Writes the first line of the thermdat file, which contains information
    on the composition, phase, and temperature ranges

    Parameters
    ----------
        nasa_specie : :class:`~pmutt.empirical.nasa.Nasa`
            Nasa specie to take information from
        write_date : bool, optional
            Whether or not the date should be written. If False, writes the
            first 8 characters of ``notes`` attribute. Defaults to True
    Returns
    -------
        line : str
            Thermdat line
    """
    element_pos = [
        24,  # Element 1
        28,  # Element 1#
        29,  # Element 2
        33,  # Element 2#
        34,  # Element 3
        38,  # Element 3#
        39,  # Element 4
        43]  # Element 4#
    temperature_pos = [
        44,  # Phase
        45,  # T_low
        55,  # T_high
        65,  # T_mid
        79]  # Line num

    # Adjusts the position based on the number of elements
    line1_pos = [16]
    for element, val in nasa_specie.elements.items():
        if val > 0.:
            two_digit = len(str(val)) - 1
            line1_pos.append(element_pos.pop(0))
            line1_pos.append(element_pos.pop(0) - two_digit)
    line1_pos.extend(temperature_pos)

    # Creating a list of the text to insert
    if write_date:
        now = datetime.now()
        notes = now.strftime('%Y%m%d')
    elif (nasa_specie.notes is None) or (nasa_specie.notes == ''):
        notes = ''
    else:
        notes = nasa_specie.notes[:8]

    line1_fields = [nasa_specie.name, notes]
    for element, val in nasa_specie.elements.items():
        if val > 0.:
            line1_fields.extend([element, '%d' % val])
    line1_fields.extend([nasa_specie.phase,
                         '%.1f' % nasa_specie.T_low,
                         '%.1f' % nasa_specie.T_high,
                         '%.1f' % nasa_specie.T_mid])

    # Write the content with appropriate spacing
    line = ''
    for pos, field in zip(line1_pos, line1_fields):
        line += field
        line = _insert_space(pos, line)
    line += '1\n'
    return line

def _write_line2(nasa_specie):
    """Writes the second line of the thermdat file

    Parameters
    ----------
        nasa_specie : :class:`~pmutt.empirical.nasa.Nasa`
            Nasa specie to take information from
    Returns
    -------
        line : str
            Thermdat line
    """
    line = ('{: 2.8E}{: 2.8E}{: 2.8E}{: 2.8E}{: 2.8E}    2\n'
            ''.format(nasa_specie.a_high[0], nasa_specie.a_high[1],
                      nasa_specie.a_high[2], nasa_specie.a_high[3],
                      nasa_specie.a_high[4]))
    return line


def _write_line3(nasa_specie):
    """Writes the third line of the thermdat file

    Parameters
    ----------
        nasa_specie : :class:`~pmutt.empirical.nasa.Nasa`
            Nasa specie to take information from
    Returns
    -------
        line : str
            Thermdat line
    """
    line = ('{: 2.8E}{: 2.8E}{: 2.8E}{: 2.8E}{: 2.8E}    3\n'
            ''.format(nasa_specie.a_high[5], nasa_specie.a_high[6],
                      nasa_specie.a_low[0], nasa_specie.a_low[1],
                      nasa_specie.a_low[2]))
    return line


def _write_line4(nasa_specie):
    """Writes the fourth line of the thermdat file

    Parameters
    ----------
        nasa_specie : :class:`~pmutt.empirical.nasa.Nasa`
            Nasa specie to take information from
    Returns
    -------
        line : str
            Thermdat line
    """
    line = ('{: 2.8E}{: 2.8E}{: 2.8E}{: 2.8E}                   4\n'
            ''.format(nasa_specie.a_low[3], nasa_specie.a_low[4],
                      nasa_specie.a_low[5], nasa_specie.a_low[6]))
    return line


def _insert_space(end_index, string):
    """Inserts the number of spaces required given the string and the position
    of the next non-blank field.

    Parameters:
        end_index : int
            Expected string length
        string : str
            String to add spaces to
    Returns:
        string_with_blanks : str
            String with spaces padded on the end
    """
    string += ' ' * (end_index - len(string))
    return string
