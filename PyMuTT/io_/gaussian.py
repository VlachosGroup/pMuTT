import re


def read_pattern(filename, pattern, group=0, return_immediately=True):
    """Reads the pattern from the Gaussian log file.

    Parameters
    ----------
        filename : str
            Log file
        pattern : str
            Regular expression pattern
        group : int, optional
            Group to return. Default is 0
        return_immediately : bool, optional
            If True, returns after the first instance. If False, reads the
            whole file. Default is True
    Returns
    -------
        out_values : str or list of str
            Value(s) corresponding to pattern.
            - str if return_immediately is True
            - list if return_immediately is False
            If the pattern is not found, returns an empty list
    """
    out_values = []
    with open(filename, 'r') as f_ptr:
        for line in f_ptr:
            if re.search(pattern, line):
                if return_immediately:
                    return re.search(pattern, line).groups()[group]
                else:
                    data = re.search(pattern, line)
                    data_split = (data.groups()[0]).split()
                    for value in data_split:
                        out_values.append(value)
        else:
            return out_values


def read_zpe(filename):
    """Reads the zero-point energy from the Gaussian log file.

    Parameters
    ----------
        filename : str
            Log file
    Returns
    -------
        zero_point_energy : float
            Zero point energy in Hartree/particle
    """
    return float(read_pattern(filename=filename,
                              pattern='Zero-point correction=(.*?)\(',
                              group=0,
                              return_immediately=True))


def read_electronic_and_zpe(filename):
    """Reads the electronic energy and zero-point energy from the
    Gaussian log file.

    Parameters
    ----------
        filename : str
            Log file
    Returns
    -------
        electronic_and_zero_point_energy : float
            Electronic and zero point energy in Hartree/particle
    """
    return float(read_pattern(filename=filename,
                              pattern='Sum of electronic and zero-point '
                              'Energies=(.*)',
                              group=0,
                              return_immediately=True))


def read_frequencies(filename):
    """Reads the frequencies from the Gaussian log file.

    Parameters
    ----------
        filename : str
            Log file
    Returns
    -------
        frequencies : list of float
            Frequencies in 1/cm
    """
    # freq_out = []
    freq_patterns = read_pattern(filename=filename,
                                 pattern='Frequencies -- (.*)',
                                 group=0,
                                 return_immediately=False)
    return [float(freq) for freq in freq_patterns]
    # for freq_line in freq_patterns:
    #     print(freq_line)
    #     freq_out.append(float(freq) for freq in freq_line.split())
    # return freq_out


def read_rotational_temperatures(filename):
    """Reads the rotational temperatures from the Gaussian log file.

    Parameters
    ----------
        filename : str
            Log file
    Returns
    -------
        rotational_temperatures : list of float
            Rotational temperatures in K
    """
    rot_T_patterns = read_pattern(filename=filename,
                                  pattern='Rotational temperatures '
                                  '\(Kelvin\)(.*)',
                                  group=0,
                                  return_immediately=False)
    return [float(rot_T) for rot_T in rot_T_patterns]


def read_molecular_mass(filename):
    """Reads the molecular mass from the Gaussian log file.

    Parameters
    ----------
        filename : str
            Log file
    Returns
    -------
        molecular_mass : float
            Molecular mass in amu
    """
    data = read_pattern(filename=filename,
                        pattern='Molecular mass:(.*)',
                        group=0,
                        return_immediately=True).split()
    return float(data[0])


def read_rot_symmetry_num(filename):
    """Reads the rotational symmetry number from the Gaussian log file.

    Parameters
    ----------
        filename : str
            Log file
    Returns
    -------
        rot_symmetry_num : int
            Rotational symmetry number
    """
    return int(float(read_pattern(filename=filename,
                                  pattern='Rotational symmetry number(.*)',
                                  group=0,
                                  return_immediately=True)))
