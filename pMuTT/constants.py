# -*- coding: utf-8 -*-
"""
pMuTT.constants

Contains universal constants for catalysis research
"""

import numpy as np


def R(units):
    """Universal molar gas constant, R

    Parameters
    ----------
        units : str
            Units for R. Supported units

            =================   ===============================================
            Unit                Description
            =================   ===============================================
            J/mol/K             Joule per mole per kelvin
            kJ/mol/K            kiloJoule per mole per kelvin
            L kPa/mol/K         Liter kilopascal per mole per kelvin
            cm3 kPa/mol/K       Cubic centimeter kilopascal per mole per kelvin
            m3 Pa/mol/K         Cubic meter pascal per mole per kelvin
            cm3 MPa/mol/K       Cubic centimeter megapascal per mole per kelvin
            m3 bar/mol/K        Cubic meters bar per mole per kelvin
            L bar/mol/K         Liter bar per mole per kelvin
            L torr/mol/K        Liter torr per mole per kelvin
            cal/mol/K           Calorie per mole per kelvin
            kcal/mol/K          Kilocalorie per mole per kevin
            L atm/mol/K         Liter atmosphere per mole per kelvin
            cm3 atm/mol/K       Cubic centimeter atmosphere per mole per kelvin
            eV/K                Electron volt per molecule per kelvin
            Eh/K                Hartree per molecule per kelvin
            Ha/K                Hartree per molecule per kelvin
            =================   ===============================================

    Returns
    -------
        R : float
            Universal molar gas constant in appropriate units
    Raises
    ------
        KeyError
            If units is not supported.
    """
    R_dict = {
        'J/mol/K': 8.3144598,
        'kJ/mol/K': 8.3144598e-3,
        'L kPa/mol/K': 8.3144598,
        'cm3 kPa/mol/K': 8.3144598e3,
        'm3 Pa/mol/K': 8.3144598,
        'cm3 MPa/mol/K': 8.3144598,
        'm3 bar/mol/K': 8.3144598e-5,
        'L bar/mol/K': 8.3144598e-2,
        'L torr/mol/K': 62.363577,
        'cal/mol/K': 1.9872036,
        'kcal/mol/K': 1.9872036e-3,
        'L atm/mol/K': 0.082057338,
        'cm3 atm/mol/K': 82.057338,
        'eV/K': 8.6173303e-5,
        'Eh/K': 3.1668105e-06,
        'Ha/K': 3.1668105e-06,
    }
    try:
        return R_dict[units]
    except KeyError:
        raise KeyError('Invalid unit for R: {}. Use help(pMuTT.constants.R) '
                       'for accepted units.'.format(units))


def h(units, bar=False):
    """Planck's constant, h

    Parameters
    ----------
        units : str
            Units for h. Supported units

            =================   ===============================================
            Unit                Description
            =================   ===============================================
            J s                 Joule second
            kJ s                Kilojoule second
            eV s                Electron volt second
            Eh s                Hartree second
            Ha s                Hartree second
            =================   ===============================================

        bar : bool, optional
            If True, returns h/2*pi . Default is False
    Returns
    -------
        h : float
            Planck's constant in appropriate units
    Raises
    ------
        KeyError
            If units is not supported.
    """
    h_dict = {
        'J s': 6.626070040e-34,
        'kJ s': 6.626070040e-37,
        'eV s': 4.135667662e-15,
        'Eh s': 1.519829846E-16,
        'Ha s': 1.519829846E-16,

    }

    try:
        h_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}. Use help(pMuTT.constants.h) for '
                       'accepted units.'.format(units))

    if bar:
        return h_dict[units]/(2.*np.pi)
    else:
        return h_dict[units]


def kb(units):
    """Boltzmann constant

    Parameters
    ----------
        units : str
            Units for kb. Supported units

            =================   ===============================================
            Unit                Description
            =================   ===============================================
            J/K                 Joule per kelvin
            eV/K                Electron volt per kelvin
            cal/K               Calorie per kelvin
            Eh/K                Hartree per kelvin
            Ha/K                Hartree per kelvin
            =================   ===============================================

    Returns
    -------
        kb : float
            Boltzmann constant in appropriate units
    Raises
    ------
        KeyError
            If units is not supported.
    """
    kb_dict = {
        'J/K': 1.38064852e-23,
        'kJ/K': 1.38064852e-26,
        'eV/K': 8.6173303e-5,
        'cal/K': 3.2976230e-24,
        'kcal/K': 3.2976230e-27,
        'Eh/K': 3.1668105e-06,
        'Ha/K': 3.1668105e-06,
    }
    try:
        return kb_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}. Use help(pMuTT.constants.kb) for '
                       'accepted units.'.format(units))


def c(units):
    """Speed of light

    Parameters
    ----------
        units : str
            Units for c. Supported units

            =================   ===============================================
            Unit                Description
            =================   ===============================================
            m/s                 Meter per second
            cm/s                Centimeter per second
            =================   ===============================================

    Returns
    -------
        c : float
            Speed of light in appropriate units
    Raises
    ------
        KeyError
            If units is not supported.
    """
    c_dict = {
        'm/s': 299792458.,
        'cm/s': 299792458.e2,
    }
    try:
        return c_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}. Use help(pMuTT.constants.c) for '
                       'accepted units.'.format(units))


def m_e(units):
    """Mass of an electron

    Parameters
    ----------
        units : str
            Units for m_e. Supported units

            =================   ===============================================
            Unit                Description
            =================   ===============================================
            kg                  Kilograms
            g                   Grams
            amu                 Atomic mass units
            =================   ===============================================

    Returns
    -------
        m_e : float
            Mass of electron in appropriate units
    Raises
    ------
        KeyError
            If units is not supported.
    """
    m_e_dict = {
        'kg': 9.10938356e-31,
        'g': 9.10938356e-28,
        'amu': 5.48579909070e-4
    }
    try:
        return m_e_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}. Use help(pMuTT.constants.m_e) for '
                       'accepted units.'.format(units))


def m_p(units):
    """Mass of a proton

    Parameters
    ----------
        units : str
            Units for m_p. Supported units

            =================   ===============================================
            Unit                Description
            =================   ===============================================
            kg                  Kilograms
            g                   Grams
            amu                 Atomic mass units
            =================   ===============================================

    Returns
    -------
        m_p : float
            Mass of proton in appropriate units
    Raises
    ------
        KeyError
            If units is not supported.
    """
    m_p_dict = {
        'kg': 1.6726219e-27,
        'g': 1.6726219e-24,
        'amu': 1.007276466879
    }
    try:
        return m_p_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}. Use help(pMuTT.constants.m_p) for '
                       'accepted units.'.format(units))


def P0(units):
    """Reference pressure

    Parameters
    ----------
        units : str
            Units for P0. Supported units

            =================   ===============================================
            Unit                Description
            =================   ===============================================
            bar                 Bar
            atm                 Atmosphere
            Pa                  Pascal
            kPa                 Kilopascal
            MPa                 Megapascal
            psi                 Pound per square inch
            mmHg                Millimeter of mercury
            Torr                Torr
            =================   ===============================================

    Returns
    -------
        P0 : float
            Reference pressure in appropriate units
    Raises
    ------
        KeyError
            If units is not supported.
    """
    P0_dict = {
        'bar': 1.01325,
        'atm': 1.,
        'Pa': 101325.,
        'kPa': 101.325,
        'MPa': 0.101325,
        'psi': 14.6959,
        'mmHg': 760.,
        'Torr': 760.,
        None: 1
    }
    try:
        return P0_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}. Use help(pMuTT.constants.P0) for '
                       'accepted units.'.format(units))


def T0(units):
    """Room temperature

    Parameters
    ----------
        units : str
            Units for T0. Supported units

            =================   ===============================================
            Unit                Description
            =================   ===============================================
            K                   Kelvin
            C                   Celcius
            R                   Rankine
            F                   Fahrenheit
            =================   ===============================================

    Returns
    -------
        T0 : float
            Room temperature in appropriate units
    Raises
    ------
        KeyError
            If units is not supported.
    """
    T0_dict = {'K': 298.15,
               'C': 25.,
               'R': 533.07,
               'F': 73.4
               }
    try:
        return T0_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}. Use help(pMuTT.constants.T0) for '
                       'accepted units.'.format(units))


def V0(units):
    """Molar volume of ideal gas at standard temperature and pressure

    Parameters
    ----------
        units : str
            Units for V0. Supported units

            ========= =================
            Symbol    Unit
            ========= =================
            m3        Metre cubed
            cm3       Centimeter cubed
            mL        Milliliters
            L         Liters
            ========= =================
    """
    V0 = R('J/mol/K')*T0('K')/P0('Pa')
    return V0*convert_unit(from_='m3', to=units)


Na = 6.02214086e23  # Avogadro number
"""float : Avogadro's number in molecules/mol"""

e = 1.6021766208e-19  # Electronic charge
"""float : Electronic charge in Coulombs"""


def convert_unit(num=None, from_=None, to=None):
    """Converts units between two unit sets

    Parameters
    ----------
        num : float, optional
            Number to convert. I not specified, will return the appropriate
            conversion factor.
        from_ : str
            Units that num is currently in
        to : str
            Units you would like num to be in
    Returns
    -------
        conversion_num : float
            num in the appropriate units
    Raises
    ------
        ValueError
            If unit types are not consistent or not supported

    **Supported Units**

    *Energy*

    ========= =================
    Symbol    Unit
    ========= =================
    J         Joules
    kJ        KiloJoules
    eV        Electron Volts
    cal       Calories
    kcal      Kilocalories
    L atm     Liter atmospheres
    Eh        Hartree
    Ha        Hartree
    ========= =================

    *Energy/Amount*

    ===========      =====================
    Symbol           Unit
    ===========      =====================
    J/mol            Joules per mole
    kJ/mol           KiloJoules per mole
    cal/mol          Calories per mole
    kcal/mol         Kilocalories per mole
    eV/molecule      eV per molecule
    Eh/molecule      Hartree per molecule
    Ha/molecule      Hartree per molecule
    ===========      =====================

    *Time*

    ========= =================
    Symbol    Unit
    ========= =================
    s         Seconds
    min       Minutes
    hr        Hours
    ========= =================

    *Amount*

    ========= =================
    Symbol    Unit
    ========= =================
    molecule  Molecules
    mol       Moles
    ========= =================

    *Temperature*

    ========= =================
    Symbol    Unit
    ========= =================
    C         Celcius
    K         Kelvin
    ========= =================

    *Length*

    ========= =================
    Symbol    Unit
    ========= =================
    m         Meter
    cm        Centimeter
    nm        Nanometer
    A         Anstrom
    ========= =================

    *Area*

    ========= =================
    Symbol    Unit
    ========= =================
    m2        Meters squared
    cm2       Centimeter squared
    A2        Anstrom squared
    ========= =================

    *Volume*

    ========= =================
    Symbol    Unit
    ========= =================
    m3        Metre cubed
    cm3       Centimeter cubed
    mL        Milliliters
    L         Liters
    ========= =================

    *Mass*

    ========= =================
    Symbol    Unit
    ========= =================
    kg        Kilograms
    g         Grams
    amu       Atomic mass units
    ========= =================

    *Pressure*

    ========= =======================
    Symbol    Unit
    ========= =======================
    Pa        Pascals
    kPa       KiloPascals
    MPa       MegaPascals
    atm       Atmospheres
    bar       Bars
    mmHg      Millilmeters of Mercury
    psi       Pounds per square inch
    ========= =======================

    """

    type_dict = {
        'J': 'energy',
        'kJ': 'energy',
        'eV': 'energy',
        'cal': 'energy',
        'kcal': 'energy',
        'L atm': 'energy',
        'Eh': 'energy',
        'Ha': 'energy',
        'J/mol': 'energy/amount',
        'kJ/mol': 'energy/amount',
        'cal/mol': 'energy/amount',
        'kcal/mol': 'energy/amount',
        'eV/molecule': 'energy/amount',
        'Eh/molecule': 'energy/amount',
        'Ha/molecule': 'energy/amount',
        's': 'time',
        'min': 'time',
        'hr': 'time',
        'molecule': 'amount',
        'mol': 'amount',
        'C': 'temp',
        'K': 'temp',
        'm': 'length',
        'cm': 'length',
        'nm': 'length',
        'A': 'length',
        'm2': 'area',
        'cm2': 'area',
        'A2': 'area',
        'm3': 'volume',
        'cm3': 'volume',
        'mL': 'volume',
        'L': 'volume',
        'kg': 'mass',
        'g': 'mass',
        'amu': 'mass',
        'Pa': 'pressure',
        'kPa': 'pressure',
        'MPa': 'pressure',
        'atm': 'pressure',
        'bar': 'pressure',
        'mmHg': 'pressure',
        'torr': 'pressure',
        'psi': 'pressure'
    }

    unit_dict = {
        'J': 1.,
        'kJ': 1.e-3,
        'eV': 6.242e+18,
        'cal': 0.239006,
        'kcal': 0.000239006,
        'L atm': 101.33,
        'Eh': 4.359744650e-18,
        'Ha': 4.359744650e-18,
        'J/mol': 1.,
        'kJ/mol': 1.e-3,
        'cal/mol': 0.239006,
        'kcal/mol': 0.000239006,
        'eV/molecule': 6.242e+18/6.02214086e23,
        's': 1.,
        'min': 1./60.,
        'hr': 1./3600.,
        'mol': 1.,
        'molecule': 6.02214086e23,
        'C': 0.,
        'K': 273.15,
        'm': 1.,
        'cm': 100.,
        'nm': 1.e9,
        'A': 1.e10,
        'm2': 1.,
        'cm2': 1.e4,
        'A2': 1.e20,
        'm3': 1.,
        'cm3': 1.e6,
        'mL': 1.e6,
        'L': 1.e3,
        'kg': 1.,
        'g': 1.e3,
        'amu': 6.022e+26,
        'Pa': 1.,
        'kPa': 1.e-3,
        'MPa': 1.e-6,
        'atm': 9.86923e-6,
        'bar': 1.e-5,
        'mmHg': 0.00750062,
        'torr': 0.00750062,
        'psi': 0.000145038
    }

    # Check if the entry exists
    if type_dict.get(from_) is None:
        raise ValueError("%r not a supported unit. Use help(pMuTT.constants."
                         "convert_unit) for accepted units." % from_)
    if type_dict.get(to) is None:
        raise ValueError("%r not a supported unit. Use help(pMuTT.constants."
                         "convert_unit) for accepted units." % to)
    # Check that the unit types are the same
    from_type = type_dict[from_]
    to_type = type_dict[to]
    if from_type != to_type:
        raise ValueError("%r [Type %r] not compatible with %r [Type %r]. "
                         "Use help(pMuTT.constants.convert_unit) for "
                         "accepted units." % (from_, from_type, to, to_type))
    elif from_type == 'temp':
        if num is None:
            num = 0.
        return num + unit_dict[to] - unit_dict[from_]
    else:
        if num is None:
            num = 1.
        return num * unit_dict[to] / unit_dict[from_]


def energy_to_freq(energy):
    """Converts energy to frequency

    Parameters
    ----------
        energy : float
            Energy in J
    Returns
    -------
        freq : float
            Frequency in Hz
    """
    return energy/h('J s')


def energy_to_temp(energy):
    """Converts energy to temperature

    Parameters
    ----------
        energy : float
            Energy in J
    Returns
    -------
        temp : float
            Temperature in K
    """
    return energy/kb('J/K')


def energy_to_wavenumber(energy):
    """Converts energy to wavenumber

    Parameters
    ----------
        energy : float
            Energy in J
    Returns
    -------
        wavenumber : float
            Wavenumber in 1/cm
    """
    return energy/h('J s')/c('cm/s')


def freq_to_energy(freq):
    """Converts frequency to energy

    Parameters
    ----------
        freq : float
            Frequency in Hz
    Returns
    -------
        energy : float
            Energy in J
    """
    return freq*h('J s')


def freq_to_temp(freq):
    """Converts frequency to temperature

    Parameters
    ----------
        freq : float
            Frequency in Hz
    Returns
    -------
        temp : float
            Temperature in K
    """
    return freq*h('J s')/kb('J/K')


def freq_to_wavenumber(freq):
    """Converts frequency to wavenumber

    Parameters
    ----------
        freq : float
            Frequency in Hz
    Returns
    -------
        wavenumber : float
            Wavenumber in 1/cm
    """
    return freq/c('cm/s')


def inertia_to_temp(inertia):
    """Converts moment of inertia into rotational temperature

    Parameters
    ----------
        inertia : float
            Moment of inertia in kg*m2
    Returns
    -------
        rot_temperature : float
            Rotational temperature in K
    """
    return h('eV s', bar=True)**2/2./kb('eV/K')/inertia \
        * convert_unit(from_='eV', to='J')


def temp_to_energy(temp):
    """Converts temperature to energy

    Parameters
    ----------
        temp : float
            Temperature in K
    Returns
    -------
        energy : float
            Energy in J
    """
    return temp*kb('J/K')


def temp_to_freq(temp):
    """Converts temperature to frequency

    Parameters
    ----------
        temp : float
            Temperature in K
    Returns
    -------
        freq : float
            Frequency in Hz
    """
    return temp*kb('J/K')/h('J s')


def temp_to_wavenumber(temp):
    """Converts vibrational/rotational temperature to wavenumber

    Parameters
    ----------
        temp : float
            Temperature in K
    Returns
    -------
        wavenumber : float
            Wavenumber in 1/cm
    """
    return temp*kb('J/K')/c('cm/s')/h('J s')


def wavenumber_to_energy(wavenumber):
    """Converts wavenumbers (1/cm) to energies (eV)

    Parameters
    ----------
        wavenumber : float
            Wavenumber in 1/cm
    Returns
    -------
        energies : float
            Corresponding temperature in eV
    """
    return wavenumber*c('cm/s')*h('eV s')


def wavenumber_to_freq(wavenumber):
    """Converts wavenumber to frequency

    Parameters
    ----------
        wavenumber : float
            Wavenumber in 1/cm
    Returns
    -------
        freq : float
            Frequency in Hz
    """
    return wavenumber*c('cm/s')


def wavenumber_to_inertia(wavenumber):
    """Converts wavenumber (1/cm) to moment of inertia

    Parameters
    ----------
        wavenumber : float
            Wavenumber in 1/cm
    Returns
    -------
        mu : float
            Moment of inertia in kg*m2
    """
    return h('J s')/(8.*np.pi**2*wavenumber*c('cm/s'))


def wavenumber_to_temp(wavenumber):
    """Converts wavenumbers (1/cm) to temperatures (K)

    Parameters
    ----------
        wavenumber : float
            Wavenumber in 1/cm
    Returns
    -------
        temperature : float
            Corresponding temperature in K
    """
    return wavenumber*c('cm/s')*h('J s')/kb('J/K')


prefixes = {
    'Y': 1.e24,
    'Z': 1.e21,
    'E': 1.e18,
    'P': 1.e15,
    'T': 1.e12,
    'G': 1.e9,
    'M': 1.e6,
    'k': 1.e3,
    'm': 1.e-3,
    'mu': 1.e-9,
    'p': 1.e-12,
    'f': 1.e-15,
    'a': 1.e-18,
    'z': 1.e-21,
    'y': 1.e-24
}
"""dict : SI unit prefixes"""

atomic_weight = {
    1: 1.008,
    2: 4.002602,
    3: 6.938,
    4: 9.0121831,
    5: 10.806,
    6: 12.0116,
    7: 14.007,
    8: 15.999,
    9: 18.99840316,
    10: 20.1797,
    11: 22.98976928,
    12: 24.305,
    13: 26.9815385,
    14: 28.085,
    15: 30.973762,
    16: 32.06,
    17: 35.45,
    18: 39.948,
    19: 39.0983,
    20: 40.078,
    21: 44.955908,
    22: 47.867,
    23: 50.9415,
    24: 51.9961,
    25: 54.938044,
    26: 55.845,
    27: 58.933194,
    28: 58.6934,
    29: 63.546,
    30: 65.38,
    31: 69.723,
    32: 72.63,
    33: 74.921595,
    34: 78.971,
    35: 79.901,
    36: 83.798,
    37: 85.4678,
    38: 87.62,
    39: 88.90584,
    40: 91.224,
    41: 92.90637,
    42: 95.95,
    43: 98,
    44: 101.07,
    45: 102.9055,
    46: 106.42,
    47: 107.8682,
    48: 112.414,
    49: 114.818,
    50: 118.71,
    51: 121.76,
    52: 127.6,
    53: 126.90447,
    54: 131.293,
    55: 132.905452,
    56: 137.327,
    57: 138.90547,
    58: 140.116,
    59: 140.90766,
    60: 144.242,
    61: 145,
    62: 150.36,
    63: 151.964,
    64: 157.25,
    65: 158.92535,
    66: 162.5,
    67: 164.93033,
    68: 167.259,
    69: 168.93422,
    70: 173.054,
    71: 174.9668,
    72: 178.49,
    73: 180.94788,
    74: 183.84,
    75: 186.207,
    76: 190.23,
    77: 192.217,
    78: 195.084,
    79: 196.966569,
    80: 200.592,
    81: 204.382,
    82: 207.2,
    83: 208.9804,
    84: 209,
    85: 210,
    86: 222,
    87: 223,
    88: 226,
    89: 227,
    90: 232.0377,
    91: 231.03588,
    92: 238.02891,
    93: 237,
    94: 244,
    95: 243,
    96: 247,
    97: 247,
    98: 251,
    99: 252,
    100: 257,
    101: 258,
    102: 259,
    103: 262,
    104: 267,
    105: 268,
    106: 271,
    107: 272,
    108: 270,
    109: 276,
    110: 281,
    111: 280,
    112: 285,
    113: 284,
    114: 289,
    115: 288,
    116: 293,
    118: 294,
    'H': 1.008,
    'He': 4.002602,
    'Li': 6.938,
    'Be': 9.0121831,
    'B': 10.806,
    'C': 12.0116,
    'N': 14.007,
    'O': 15.999,
    'F': 18.99840316,
    'Ne': 20.1797,
    'Na': 22.98976928,
    'Mg': 24.305,
    'Al': 26.9815385,
    'Si': 28.085,
    'P': 30.973762,
    'S': 32.06,
    'Cl': 35.45,
    'Ar': 39.948,
    'K': 39.0983,
    'Ca': 40.078,
    'Sc': 44.955908,
    'Ti': 47.867,
    'V': 50.9415,
    'Cr': 51.9961,
    'Mn': 54.938044,
    'Fe': 55.845,
    'Co': 58.933194,
    'Ni': 58.6934,
    'Cu': 63.546,
    'Zn': 65.38,
    'Ga': 69.723,
    'Ge': 72.63,
    'As': 74.921595,
    'Se': 78.971,
    'Br': 79.901,
    'Kr': 83.798,
    'Rb': 85.4678,
    'Sr': 87.62,
    'Y': 88.90584,
    'Zr': 91.224,
    'Nb': 92.90637,
    'Mo': 95.95,
    'Tc': 98,
    'Ru': 101.07,
    'Rh': 102.9055,
    'Pd': 106.42,
    'Ag': 107.8682,
    'Cd': 112.414,
    'In': 114.818,
    'Sn': 118.71,
    'Sb': 121.76,
    'Te': 127.6,
    'I': 126.90447,
    'Xe': 131.293,
    'Cs': 132.905452,
    'Ba': 137.327,
    'La': 138.90547,
    'Ce': 140.116,
    'Pr': 140.90766,
    'Nd': 144.242,
    'Pm': 145,
    'Sm': 150.36,
    'Eu': 151.964,
    'Gd': 157.25,
    'Tb': 158.92535,
    'Dy': 162.5,
    'Ho': 164.93033,
    'Er': 167.259,
    'Tm': 168.93422,
    'Yb': 173.054,
    'Lu': 174.9668,
    'Hf': 178.49,
    'Ta': 180.94788,
    'W': 183.84,
    'Re': 186.207,
    'Os': 190.23,
    'Ir': 192.217,
    'Pt': 195.084,
    'Au': 196.966569,
    'Hg': 200.592,
    'Tl': 204.382,
    'Pb': 207.2,
    'Bi': 208.9804,
    'Po': 209,
    'At': 210,
    'Rn': 222,
    'Fr': 223,
    'Ra': 226,
    'Ac': 227,
    'Th': 232.0377,
    'Pa': 231.03588,
    'U': 238.02891,
    'Np': 237,
    'Pu': 244,
    'Am': 243,
    'Cm': 247,
    'Bk': 247,
    'Cf': 251,
    'Es': 252,
    'Fm': 257,
    'Md': 258,
    'No': 259,
    'Lr': 262,
    'Rf': 267,
    'Db': 268,
    'Sg': 271,
    'Bh': 272,
    'Hs': 270,
    'Mt': 276,
    'Ds': 281,
    'Rg': 280,
    'Cn': 285,
    'Uut': 284,
    'Fl': 289,
    'Uup': 288,
    'Lv': 293,
    'Uuo': 294,
}
"""dict : Atomic weight. The key can be the atomic number, the element symbol,
or the element name"""
