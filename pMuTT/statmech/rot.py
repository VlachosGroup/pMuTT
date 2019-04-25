# -*- coding: utf-8 -*-

import itertools
from warnings import warn
import numpy as np
from pMuTT import _pMuTTBase
from pMuTT import constants as c
from pMuTT.io.json import remove_class

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

class RigidRotor(_pMuTTBase):
    """Rotational mode using the rigid rotor assumption. Equations sourced from:
    
    * Sandler, S. I. An Introduction to Applied Statistical Thermodynamics;
      John Wiley & Sons, 2010.

    Attributes
    ----------
        symmetrynumber : float or str
            Symmetry number of molecule. If a string is inputted, it represents 
            the point group. Supported options are shown below.

            ===========    ===============
            Point group    symmetry number
            ===========    ===============
            C1             1
            Cs             1
            C2             2
            C2v            2
            C3v            3
            Cinfv          1
            D2h            4
            D3h            6
            D5h            10
            Dinfh          2
            D3d            6
            Td             12
            Oh             24
            ===========    ===============

            See DOI for more details: 10.1007/s00214-007-0328-0
        rot_temperatures : list of float, optional
            Rotational temperatures in K
        geometry : str
            Geometry of molecule. Accepted options are:

            - monatomic
            - linear
            - nonlinear
        atoms : ase.Atoms object, optional
            An atoms object can be used to calculate rot_temperatures and
            guess geometry
        degree_tol : float
            Degree tolerance to estimate geometry. Only required if estimating
            geometry or rot_temperatures
    """

    def __init__(self, symmetrynumber, rot_temperatures=None, geometry=None,
                 atoms=None, degree_tol=5.):
        if isinstance(symmetrynumber, str):
            try:
                symmetrynumber = symmetry_dict[symmetrynumber]
            except KeyError:
                raise ValueError('Point group, {}, not supported. '
                                 'See :class:`~pMuTT.statmech.rot.RigidRotor '
                                 'for supported '
                                 'options.'.format(symmetrynumber))

        self.symmetrynumber = symmetrynumber
        if rot_temperatures is None and atoms is not None:
            self.rot_temperatures = get_rot_temperatures_from_atoms(
                    atoms=atoms, degree_tol=degree_tol)
        else:
            self.rot_temperatures = rot_temperatures

        if geometry is None and atoms is not None:
            self.geometry = get_geometry_from_atoms(
                    atoms=atoms, degree_tol=degree_tol)
        else:
            self.geometry = geometry

    def get_q(self, T):
        """Calculates the partition function

        :math:`q^{rot}=0` if monatomic

        :math:`q^{rot}=\\frac{T}{\\sigma}\\prod_i\\frac{1}{\\Theta_R}` if linear

        :math:`q^{rot}=\\frac{\\sqrt{\\pi}}{\\sigma}(T^3\\prod_i\\frac{1}{
        \\Theta_{Ri}})^\\frac{1}{2}` if nonlinear

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            q_rot : float
                Rotational partition function
        """
        if self.geometry == 'monatomic':
            return 0.
        elif self.geometry == 'linear':
            return T/self.symmetrynumber/np.prod(self.rot_temperatures)
        elif self.geometry == 'nonlinear':
            return np.sqrt(np.pi)/self.symmetrynumber \
                * (T**3/np.prod(self.rot_temperatures))**0.5
        else:
            raise ValueError(
                'Geometry, {}, not supported.'.format(self.geometry))

    def get_CvoR(self):
        """Calculates the dimensionless heat capacity at constant volume

        :math:`\\frac{C_V^{rot}}{R}=0` if monatomic

        :math:`\\frac{C_V^{rot}}{R}=1` if linear

        :math:`\\frac{C_V^{rot}}{R}=1.5` if nonlinear

        Returns
        -------
            CvoR_rot : float
                Rotational dimensionless heat capacity at constant volume
        """
        if self.geometry == 'monatomic':
            return 0.
        elif self.geometry == 'linear':
            return 1.
        elif self.geometry == 'nonlinear':
            return 1.5
        else:
            raise ValueError(
                'Geometry, {}, not supported.'.format(self.geometry))

    def get_CpoR(self):
        """Calculates the dimensionless heat capacity at constant pressure

        :math:`\\frac{C_P^{rot}}{R}=\\frac{C_V^{rot}}{R}`

        Returns
        -------
            CpoR_rot : float
                Rotational dimensionless heat capacity at constant pressure
        """
        return self.get_CvoR()

    def get_UoRT(self):
        """Calculates the dimensionless internal energy

        :math:`\\frac{U^{rot}}{RT}=0` if monatomic

        :math:`\\frac{U^{rot}}{RT}=1` if linear

        :math:`\\frac{U^{rot}}{RT}=1.5` if nonlinear

        Returns
        -------
            UoRT_rot : float
                Rotational dimensionless internal energy
        """
        if self.geometry == 'monatomic':
            return 0.
        elif self.geometry == 'linear':
            return 1.
        elif self.geometry == 'nonlinear':
            return 1.5
        else:
            raise ValueError(
                'Geometry, {}, not supported.'.format(self.geometry))

    def get_HoRT(self):
        """Calculates the dimensionless enthalpy

        :math:`\\frac{U^{rot}}{RT}=\\frac{H^{rot}}{RT}`

        Returns
        -------
            HoRT_rot : float
                Rotational dimensionless enthalpy
        """
        return self.get_UoRT()

    def get_SoR(self, T):
        """Calculates the dimensionless entropy

        :math:`\\frac{S^{rot}}{R}=0` if monatomic

        :math:`\\frac{S}{R}=\\log(\\frac{T}{\\sigma}\\prod_i \\frac{1}{
        \\Theta_{R,i}})+1` if linear

        :math:`\\frac{S^{rot}}{R}=\\log\\bigg(\\frac{\\sqrt\\pi }{\\sigma}
        \\prod_i{(\\frac{T^3}{\\Theta_{R,i}})^\\frac{1}{2}}\\bigg)+1.5`
        if nonlinear

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            SoR_rot : float
                Rotational dimensionless entropy
        """
        if self.geometry == 'monatomic':
            return 0.
        elif self.geometry == 'linear':
            return np.log(T/self.symmetrynumber
                          / np.prod(self.rot_temperatures)) + 1.
        elif self.geometry == 'nonlinear':
            return np.log(np.sqrt(np.pi)/self.symmetrynumber
                          * (T**3/np.prod(self.rot_temperatures))**0.5) + 1.5
        else:
            raise ValueError(
                'Geometry, {}, not supported.'.format(self.geometry))

    def get_FoRT(self, T):
        """Calculates the dimensionless Helmholtz energy

        :math:`\\frac{A^{rot}}{RT}=\\frac{U^{rot}}{RT}-\\frac{S^{rot}}{R}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            FoRT_rot : float
                Rotational dimensionless Helmholtz energy
        """
        return self.get_UoRT()-self.get_SoR(T=T)

    def get_GoRT(self, T):
        """Calculates the dimensionless Gibbs energy

        :math:`\\frac{G^{rot}}{RT}=\\frac{H^{rot}}{RT}-\\frac{S^{rot}}{R}`

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            GoRT_rot : float
                Rotational dimensionless Gibbs energy
        """
        return self.get_HoRT()-self.get_SoR(T=T)

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        return {'class': str(self.__class__),
                'symmetrynumber': self.symmetrynumber,
                'geometry': self.geometry,
                'rot_temperatures': list(self.rot_temperatures)}

def get_rot_temperatures_from_atoms(atoms, geometry=None, degree_tol=5.):
    """Calculate the rotational temperatures from ase.Atoms object

    Parameters
    ----------
        atoms : ase.Atoms object
            Atoms object
        geometry : str, optional
            Geometry of molecule. If not specified, it will be guessed from
            Atoms object.
        degree_tol : float, optional
            Degree tolerance in degrees. Default is 5 degrees
    Returns
    -------
        rot_temperatures : list of float
            Rotational temperatures
    """
    if geometry is None:
        geometry = get_geometry_from_atoms(atoms=atoms, degree_tol=degree_tol)

    rot_temperatures = []
    for moment in atoms.get_moments_of_inertia():
        if np.isclose(0., moment):
            continue
        moment_SI = moment*c.convert_unit(initial='amu', final='kg') \
            * c.convert_unit(initial='A2', final='m2')
        rot_temperatures.append(c.inertia_to_temp(moment_SI))

    if geometry == 'monatomic':
        # Expecting all modes to be 0
        assert np.isclose(np.sum(rot_temperatures), 0.)
        return [0.]
    elif geometry == 'linear':
        # Expecting one mode to be 0 and the other modes to be identical
        if not np.isclose(rot_temperatures[0], rot_temperatures[1]):
            warn('Expected rot_temperatures for linear specie, {}, to be '
                 'similar. Values found were:{}'
                 .format(atoms, rot_temperatures))
        return [max(rot_temperatures)]
    elif geometry == 'nonlinear':
        # Expecting 3 modes. May or may not be equal
        return rot_temperatures
    else:
        raise ValueError(
            'Geometry, {}, not supported.'.format(geometry))


def get_geometry_from_atoms(atoms, degree_tol=5.):
    """Estimate the geometry using the ase.Atoms object

    Parameters
    ----------
        atoms : ase.Atoms object
            Atoms object
        degree_tol : float, optional
            Degree tolerance in degrees. Default is 5 degrees
    Returns
    -------
        geometry : str
            Geometry
    """
    if len(atoms) == 1:
        return 'monatomic'
    elif len(atoms) == 2:
        return 'linear'
    else:
        for i, j, k in itertools.combinations(range(len(atoms)), 3):
            angle = atoms.get_angle(i, j, k)
            # If the angles are NOT close to 0 and 180 degrees
            if not np.isclose(angle, 0., atol=degree_tol) \
               and not np.isclose(angle, 180., atol=degree_tol):
                return 'nonlinear'
        else:
            return 'linear'
