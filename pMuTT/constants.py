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

            =============  ===============================================  ============
            Unit           Description                                      Value
            =============  ===============================================  ============
            J/mol/K        Joule per mole per kelvin                        8.3144598
            kJ/mol/K       kiloJoule per mole per kelvin                    8.3144598e-3
            L kPa/mol/K    Liter kilopascal per mole per kelvin             8.3144598
            cm3 kPa/mol/K  Cubic centimeter kilopascal per mole per kelvin  8.3144598e3
            m3 Pa/mol/K    Cubic meter pascal per mole per kelvin           8.3144598
            cm3 MPa/mol/K  Cubic centimeter megapascal per mole per kelvin  8.3144598
            m3 bar/mol/K   Cubic meters bar per mole per kelvin             8.3144598e-5
            L bar/mol/K    Liter bar per mole per kelvin                    8.3144598e-2
            L torr/mol/K   Liter torr per mole per kelvin                   62.363577
            cal/mol/K      Calorie per mole per kelvin                      1.9872036
            kcal/mol/K     Kilocalorie per mole per kevin                   1.9872036e-3
            L atm/mol/K    Liter atmosphere per mole per kelvin             0.082057338
            cm3 atm/mol/K  Cubic centimeter atmosphere per mole per kelvin  82.057338
            eV/K           Electron volt per molecule per kelvin            8.6173303e-5
            Eh/K           Hartree per molecule per kelvin                  3.1668105e-6
            Ha/K           Hartree per molecule per kelvin                  3.1668105e-6
            =============  ===============================================  ============

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

            ====  ====================  =================
            Unit  Description           Value (bar=False)
            ====  ====================  =================
            J s   Joule second          6.626070040e-34
            kJ s  Kilojoule second      6.626070040e-37
            eV s  Electron volt second  4.135667662e-15
            Eh s  Hartree second        1.519829846E-16
            Ha s  Hartree second        1.519829846E-16
            ====  ====================  =================

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

            ======  =========================  ==============
            Unit    Description                Value
            ======  =========================  ==============
            J/K     Joule per kelvin           1.38064852e-23
            kJ/K    Kilojoule per kelvin       1.38064852e-26
            eV/K    Electron volt per kelvin   8.6173303e-5
            cal/K   Calorie per kelvin         3.2976230e-24
            kcal/K  Kilocalorie per kelvin     3.2976230e-27
            Eh/K    Hartree per kelvin         3.1668105e-6
            Ha/K    Hartree per kelvin         3.1668105e-6
            ======  =========================  ==============

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

            =====  =====================  ===========
            Unit   Description            Value
            =====  =====================  ===========
            m/s    Meter per second       299792458
            cm/s   Centimeter per second  299792458e2
            =====  =====================  ===========

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

            =====  =================  ================
            Unit   Description        Value
            =====  =================  ================
            kg     Kilograms          9.10938356e-31
            g      Grams              9.10938356e-28
            amu    Atomic mass units  5.48579909070e-4
            =====  =================  ================

    Returns
    -------
        m_e : float
            Mass of electron in appropriate units
    Raises
    ------
        KeyError
            If units is not supported.
    """
    try:
        return convert_unit(num=5.48579909070e-4, initial='amu', final=units)
    except KeyError:
        raise KeyError('Invalid unit: {}. Use help(pMuTT.constants.m_p) for '
                       'accepted units.'.format(units))

def m_p(units):
    """Mass of a proton

    Parameters
    ----------
        units : str
            Units for m_p. Supported units

            =====   =================  ==============
            Unit    Description        Value
            =====   =================  ==============
            kg      Kilograms          1.6726219e-27
            g       Grams              1.6726219e-24
            amu     Atomic mass units  1.007276466879
            =====   =================  ==============

    Returns
    -------
        m_p : float
            Mass of proton in appropriate units
    Raises
    ------
        KeyError
            If units is not supported.
    """
    try:
        return convert_unit(num=1.007276466879, initial='amu', final=units)
    except KeyError:
        raise KeyError('Invalid unit: {}. Use help(pMuTT.constants.m_p) for '
                       'accepted units.'.format(units))

def P0(units):
    """Reference pressure

    Parameters
    ----------
        units : str
            Units for P0. Supported units

            ====  =====================  ========
            Unit  Description            Value
            ====  =====================  ========
            bar   Bar                    1
            atm   Atmosphere             0.986923
            Pa    Pascal                 100000
            kPa   Kilopascal             100
            MPa   Megapascal             0.1
            psi   Pound per square inch  14.5038
            mmHg  Millimeter of mercury  750.062
            Torr  Torr                   750.062
            ====  =====================  ========

    Returns
    -------
        P0 : float
            Reference pressure in appropriate units
    Raises
    ------
        KeyError
            If units is not supported.
    """
    try:
        return convert_unit(num=1., initial='bar', final=units)
    except KeyError:
        raise KeyError('Invalid unit: {}. Use help(pMuTT.constants.P0) for '
                       'accepted units.'.format(units))


def T0(units):
    """Room temperature

    Parameters
    ----------
        units : str
            Units for T0. Supported units

            ====   ===========  ======
            Unit   Description  Value
            ====   ===========  ======
            K      Kelvin       298.15
            C      Celcius      25
            R      Rankine      536.67
            F      Fahrenheit   77
            ====   ===========  ======

    Returns
    -------
        T0 : float
            Room temperature in appropriate units
    Raises
    ------
        KeyError
            If units is not supported.
    """
    try:
        return convert_unit(num=298.15, initial='K', final=units)
    except KeyError:
        raise KeyError('Invalid unit: {}. Use help(pMuTT.constants.T0) for '
                       'accepted units.'.format(units))


def V0(units):
    """Molar volume of ideal gas at standard temperature and pressure

    Parameters
    ----------
        units : str
            Units for V0. Supported units

            ====  ================  ===================
            Unit  Description       Value
            ====  ================  ===================
            m3    Metre cubed       0.02446539540458919
            cm3   Centimeter cubed  24465.395404589188
            mL    Milliliters       24465.395404589188
            L     Liters            24.465395404589188
            ====  ================  ===================
    """
    V0 = R('J/mol/K')*T0('K')/P0('Pa')
    return V0*convert_unit(initial='m3', final=units)


Na = 6.02214086e23  # Avogadro number
"""float : Avogadro's number in molecules/mol"""

e = 1.6021766208e-19  # Electronic charge
"""float : Electronic charge in Coulombs"""


def convert_unit(num=None, initial=None, final=None):
    """Converts units between two unit sets

    Parameters
    ----------
        num : float, optional
            Number to convert. I not specified, will return the appropriate
            conversion factor.
        initial : str
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
    Unit      Description
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

    ===========      ============================
    Unit             Description
    ===========      ============================
    J/mol            Joules per mole
    kJ/mol           KiloJoules per mole
    cal/mol          Calories per mole
    kcal/mol         Kilocalories per mole
    eV/molecule      Electron volt per molecule
    Eh/molecule      Hartree per molecule
    Ha/molecule      Hartree per molecule
    ===========      ============================

    *Time*

    ========= =================
    Unit      Description
    ========= =================
    ps        Picosecond
    ns        Nanosecond
    ms        Millisecond
    s         Seconds
    min       Minutes
    hr        Hours
    day       Days
    yr        Year
    ========= =================

    *Amount*

    ========= =================
    Unit      Description
    ========= =================
    molecule  Molecules
    mol       Moles
    ========= =================

    *Temperature*

    ========= =================
    Unit      Description
    ========= =================
    C         Celcius
    K         Kelvin
    F         Fahrenheit
    R         Rankine
    ========= =================

    *Length*

    ========= =================
    Unit      Description
    ========= =================
    m         Meter
    cm        Centimeter
    nm        Nanometer
    A         Anstrom
    km        Kilometer
    inch      Inch
    ft        Foot
    mile      Mile
    ========= =================

    *Area*

    ========= ===================
    Unit      Description
    ========= ===================
    m2        Meters squared
    cm2       Centimeters squared
    A2        Anstroms squared
    km2       Kilometers squared
    inch2     Inches squared
    ft2       Feet squared
    ========= ===================

    *Volume*

    ========= =================
    Unit      Description
    ========= =================
    m3        Meter cubed
    cm3       Centimeter cubed
    mL        Milliliters
    L         Liters
    inch3     Inches cubed
    ft3       Feet cubed
    ========= =================

    *Mass*

    ========= =================
    Unit      Description
    ========= =================
    kg        Kilograms
    g         Grams
    amu       Atomic mass units
    lbs       Pounds
    ========= =================

    *Pressure*

    ========= =======================
    Unit      Description
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
        'ps': 'time',
        'ns': 'time',
        'ms': 'time',
        's': 'time',
        'min': 'time',
        'hr': 'time',
        'day': 'time',
        'molecule': 'amount',
        'mol': 'amount',
        'C': 'temp',
        'K': 'temp',
        'F': 'temp',
        'R': 'temp',
        'm': 'length',
        'cm': 'length',
        'nm': 'length',
        'km': 'length',
        'inch': 'length',
        'ft': 'length',
        'mile': 'length',
        'A': 'length',
        'm2': 'area',
        'cm2': 'area',
        'A2': 'area',
        'km2': 'area',
        'inch2': 'area',
        'ft2': 'area',
        'm3': 'volume',
        'cm3': 'volume',
        'mL': 'volume',
        'L': 'volume',
        'inch3': 'volume',
        'ft3': 'volume',
        'kg': 'mass',
        'g': 'mass',
        'amu': 'mass',
        'lbs': 'mass',
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
        'ps': 1.e12,
        'ns': 1.e9,
        'ms': 1.e3,
        's': 1.,
        'min': 1./60.,
        'hr': 1./3600.,
        'day': 1./3600./24.,
        'yr': 1./3600./24./365.25,
        'mol': 1.,
        'molecule': 6.02214086e23,
        'm': 1.,
        'cm': 100.,
        'nm': 1.e9,
        'km': 1.e-3,
        'inch': 39.3701,
        'ft': 3.28084,
        'mile': 1./1609.344,
        'A': 1.e10,
        'm2': 1.,
        'cm2': 1.e4,
        'A2': 1.e20,
        'km2': 1.e-6,
        'inch2': 1550.,
        'ft2': 10.7639,
        'm3': 1.,
        'cm3': 1.e6,
        'mL': 1.e6,
        'L': 1.e3,
        'inch3': 61023.7,
        'ft3': 35.3147,
        'kg': 1.,
        'g': 1.e3,
        'amu': 6.022e+26,
        'lbs': 2.20462,
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
    try:
        initial_type = type_dict[initial]
    except KeyError:
        raise ValueError("%r not a supported unit. Use help(pMuTT.constants."
                         "convert_unit) for accepted units." % initial)
    try:
        final_type = type_dict[final]
    except KeyError:
        raise ValueError("%r not a supported unit. Use help(pMuTT.constants."
                         "convert_unit) for accepted units." % final)

    # Check that the unit types are the same
    if initial_type != final_type:
        raise ValueError("%r [Type %r] not compatible with %r [Type %r]. "
                         "Use help(pMuTT.constants.convert_unit) for "
                         "accepted units." % (initial, initial_type, final,
                                              final_type))
    elif initial_type == 'temp':
        if num is None:
            num = 0.
        # Evaluating each combination
        if initial == final:
            result = num
        elif initial == 'C':
            if final == 'K':
                result = num + 273.15
            elif final == 'F':
                result = (num * 9./5.) + 32.
            elif final == 'R':
                result = (num + 273.15) * 1.8
        elif initial == 'K':
            if final == 'C':
                result = num - 273.15
            elif final == 'F':
                result = (num * 1.8) - 459.67
            elif final == 'R':
                result = num * 1.8
        elif initial == 'F':
            if final == 'C':
                result = (num - 32.) * 5./9.
            elif final == 'K':
                result = (num + 459.67)/1.8
            elif final == 'R':
                result = num + 459.67
        elif initial == 'R':
            if final == 'C':
                result = num/1.8 - 273.15
            elif final == 'K':
                result = num/1.8
            elif final == 'F':
                result = num - 459.67
    else:
        if num is None:
            num = 1.
        result = num * unit_dict[final] / unit_dict[initial]
    return result


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
        * convert_unit(initial='eV', final='J')


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


def debye_to_einstein(debye_temperature):
    """Converts Debye temperature to Einstein temperature

    Parameters
    ----------
        debye_temperature : float
            Debye temperature in K
    Returns
    -------
        einstein_temperature : float
            Einstein temperature in K
    """
    return (np.pi/6.)**(1./3.)*debye_temperature


def einstein_to_debye(einstein_temperature):
    """Converts Einstein temperature to Debye temperature

    Parameters
    ----------
        einstein_temperature : float
            Einstein temperature in K
    Returns
    -------
        debye_temperature : float
            Debye temperature in K
    """
    return einstein_temperature/(np.pi/6.)**(1./3.)

symmetry_dict = {
    'C1': 1,
    'Cs': 1,
    'C2': 2,
    'C2v': 2,
    'C3v': 3,
    'Cinfv': 1,
    'D2h': 4,
    'D3h': 6,
    'D5h': 10,
    'Dinfh': 2,
    'D3d': 6,
    'Td': 12,
    'Oh': 24
}
"""dict : Keys are point groups and the values are the symmetry numbers used for
rotational modes."""

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
    43: 98.,
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
    84: 209.,
    85: 210.,
    86: 222.,
    87: 223.,
    88: 226.,
    89: 227.,
    90: 232.0377,
    91: 231.03588,
    92: 238.02891,
    93: 237.,
    94: 244.,
    95: 243.,
    96: 247.,
    97: 247.,
    98: 251.,
    99: 252.,
    100: 257.,
    101: 258.,
    102: 259.,
    103: 262.,
    104: 267.,
    105: 268.,
    106: 271.,
    107: 272.,
    108: 270.,
    109: 276.,
    110: 281.,
    111: 280.,
    112: 285.,
    113: 284.,
    114: 289.,
    115: 288.,
    116: 293.,
    118: 294.,
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
    'Tc': 98.,
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
    'Pm': 145.,
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
    'Po': 209.,
    'At': 210.,
    'Rn': 222.,
    'Fr': 223.,
    'Ra': 226.,
    'Ac': 227.,
    'Th': 232.0377,
    'Pa': 231.03588,
    'U': 238.02891,
    'Np': 237.,
    'Pu': 244.,
    'Am': 243.,
    'Cm': 247.,
    'Bk': 247.,
    'Cf': 251.,
    'Es': 252.,
    'Fm': 257.,
    'Md': 258.,
    'No': 259.,
    'Lr': 262.,
    'Rf': 267.,
    'Db': 268.,
    'Sg': 271.,
    'Bh': 272.,
    'Hs': 270.,
    'Mt': 276.,
    'Ds': 281.,
    'Rg': 280.,
    'Cn': 285.,
    'Uut': 284.,
    'Fl': 289.,
    'Uup': 288.,
    'Lv': 293.,
    'Uuo': 294.,
}
"""dict : Atomic weight. The key can be the atomic number, the element symbol,
or the element name"""

symmetry_dict = {
   '10025-78-2': 3,
   '10025-87-3': 2,
   '10026-04-7': 12,
   '10026-11-6': 12,
   '10026-13-8': 6,
   '10031-24-0': 2,
   '10038-98-9': 12,
   '10405-27-3': 1,
   '10545-58-1': 2,
   '10545-99-0': 2,
   '1066-43-9': 3,
   '106-93-4': 2,
   '106-96-7': 1,
   '106-97-8': 2,
   '106-99-0': 2,
   '107-02-8': 1,
   '107-04-0': 1,
   '107-06-2': 2,
   '1070-74-2': 2,
   '107-12-0': 1,
   '107-22-2': 2,
   '107-31-3': 1,
   '1076-43-3': 12,
   '109-77-3': 2,
   '110-00-9': 2,
   '110-02-1': 2,
   '110-82-7': 6,
   '1111-89-3': 3,
   '115-10-6': 2,
   '115-11-7': 2,
   '116-14-3': 4,
   '123-91-1': 2,
   '127-18-4': 4,
   '13025-73-5': 3,
   '13450-92-5': 12,
   '13455-01-1': 3,
   '13465-71-9': 3,
   '13465-73-1': 3,
   '13465-74-2': 3,
   '13465-75-3': 2,
   '13465-76-4': 3,
   '13465-78-6': 3,
   '13465-84-4': 12,
   '13478-20-1': 3,
   '13536-94-2': 2,
   '13537-03-6': 3,
   '13537-07-0': 12,
   '13537-30-9': 3,
   '13537-33-2': 3,
   '13550-49-7': 3,
   '13569-43-2': 3,
   '13637-65-5': 3,
   '13637-87-1': 1,
   '13768-94-0': 2,
   '13769-76-1': 1,
   '13770-22-4': 1,
   '13777-22-5': 12,
   '13777-23-6': 12,
   '13777-25-8': 12,
   '13780-57-9': 4,
   '13818-89-8': 6,
   '13986-26-0': 12,
   '143-36-2': 3,
   '1449-65-6': 3,
   '1455-13-6': 1,
   '1482-88-8': 3,
   '14940-63-7': 1,
   '1495-50-7': 1,
   '151-56-4': 1,
   '1517-52-8': 6,
   '15586-97-7': 3,
   '156-59-2': 2,
   '156-60-5': 2,
   '15930-75-3': 3,
   '1630-77-9': 2,
   '1632-89-9': 1,
   '1632-99-1': 6,
   '1664-98-8': 2,
   '16921-96-3': 10,
   '17178-58-4': 1,
   '1735-17-7': 6,
   '1849-29-2': 1,
   '18820-63-8': 1,
   '2003-31-8': 1,
   '2003-32-9': 1,
   '2031-95-0': 3,
   '2036-39-7': 2,
   '20427-56-9': 12,
   '2206-26-0': 3,
   '2210-34-6': 1,
   '2314-97-8': 3,
   '2404-52-6': 3,
   '25604-70-0': 1,
   '25604-71-1': 1,
   '2614-35-9': 1,
   '26395-29-9': 1,
   '2699-79-8': 2,
   '287-23-0': 4,
   '2875-94-7': 2,
   '2875-95-8': 2,
   '2875-96-9': 2,
   '288-37-9': 2,
   '3017-23-0': 1,
   '3114-46-3': 1,
   '353-36-6': 1,
   '353-49-1': 1,
   '353-50-4': 2,
   '353-54-8': 3,
   '353-58-2': 1,
   '353-85-5': 3,
   '37230-84-5': 12,
   '39130-85-3': 3,
   '3982-91-0': 3,
   '4109-96-0': 2,
   '4122-13-8': 1,
   '420-32-6': 2,
   '460-12-8': 2,
   '460-19-5': 2,
   '463-49-0': 4,
   '463-71-8': 2,
   '50-00-0': 2,
   '503-17-3': 6,
   '503-28-6': 2,
   '504-64-3': 2,
   '506-78-5': 1,
   '506-82-1': 6,
   '507-16-4': 1,
   '507-25-5': 12,
   '544-97-8': 6,
   '557-99-3': 1,
   '558-13-4': 12,
   '558-20-3': 12,
   '558-21-4': 3,
   '56-23-5': 12,
   '593-53-3': 3,
   '593-61-3': 1,
   '593-63-5': 1,
   '593-74-8': 6,
   '593-75-9': 3,
   '593-95-3': 2,
   '594-15-0': 3,
   '594-18-3': 2,
   '594-73-0': 6,
   '6088-90-0': 1,
   '6088-91-1': 1,
   '6089-44-7': 1,
   '6111-63-3': 3,
   '624-61-3': 2,
   '624-65-7': 1,
   '624-74-8': 2,
   '64-18-6': 1,
   '64-19-7': 1,
   '6552-57-4': 2,
   '659-86-9': 1,
   '661-54-1': 3,
   '666-52-4': 2,
   '673-93-8': 3,
   '67-56-1': 1,
   '67-64-1': 2,
   '676-49-3': 3,
   '676-55-1': 2,
   '67-66-3': 3,
   '676-80-2': 3,
   '67-72-1': 6,
   '683-73-8': 4,
   '71-43-2': 12,
   '7446-70-0': 3,
   '74-82-8': 12,
   '74-83-9': 3,
   '74-84-0': 6,
   '74-85-1': 4,
   '74-86-2': 2,
   '74-87-3': 3,
   '74-88-4': 3,
   '74-89-5': 1,
   '74-90-8': 1,
   '74-95-3': 2,
   '74-96-4': 1,
   '74-97-5': 1,
   '74-98-6': 2,
   '74-99-7': 3,
   '75-00-3': 1,
   '75-05-8': 3,
   '75-07-0': 1,
   '75-09-2': 2,
   '75-19-4': 6,
   '75-21-8': 2,
   '75-25-2': 3,
   '75-35-4': 2,
   '75-44-5': 2,
   '75-46-7': 3,
   '7550-45-0': 12,
   '75-61-6': 2,
   '75-62-7': 3,
   '75-63-8': 3,
   '75-69-4': 3,
   '75-71-8': 2,
   '7572-29-4': 2,
   '75-72-9': 3,
   '75-73-0': 12,
   '76-16-4': 6,
   '7646-78-8': 12,
   '7647-18-9': 6,
   '7647-19-0': 6,
   '7664-41-7': 3,
   '7720-83-4': 12,
   '7727-18-6': 3,
   '7732-18-5': 2,
   '7758-95-4': 2,
   '7772-99-8': 2,
   '7782-65-2': 12,
   '7782-79-8': 1,
   '7783-41-7': 2,
   '7783-42-8': 1,
   '7783-46-2': 2,
   '7783-47-3': 2,
   '7783-54-2': 3,
   '7783-55-3': 3,
   '7783-61-1': 12,
   '7783-75-7': 24,
   '7783-77-9': 24,
   '7783-79-1': 24,
   '7783-80-4': 24,
   '7783-81-5': 24,
   '7783-82-6': 24,
   '7784-35-2': 3,
   '7784-36-3': 6,
   '7784-42-1': 3,
   '7784-45-4': 3,
   '7787-71-5': 2,
   '7789-20-0': 2,
   '7789-57-3': 3,
   '7789-59-5': 3,
   '7789-66-4': 12,
   '7789-67-5': 12,
   '7789-68-6': 12,
   '7790-91-2': 2,
   '7791-25-5': 2,
   '7803-51-2': 3,
   '7803-52-3': 3,
   '7803-62-5': 12,
   '78-93-3': 1,
   '79-20-9': 1,
   '79-28-7': 4,
   '79-35-6': 2,
   '811-98-3': 1,
   '819-01-2': 3,
   '917-96-4': 3,
   '920-42-3': 1,
   '992-94-9': 3
}
"""dict : Symmetry number. The key is the CAS-ID"""