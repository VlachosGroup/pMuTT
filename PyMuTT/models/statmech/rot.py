# -*- coding: utf-8 -*-

import itertools
import numpy as np
from PyMuTT import constants as c

class IdealRot:
    """Rotational mode using the ideal gas assumption

    Attributes
    ----------
        symmetrynumber : float
            Symmetry number of molecule. 

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
    """

    def __init__(self, symmetrynumber, rot_temperatures=None, geometry=None):
        self.symmetrynumber = symmetrynumber
        self.rot_temperatures = rot_temperatures
        self.geometry = geometry
                    
    @classmethod
    def from_atoms(cls, symmetrynumber, atoms):
        """Calculates the rotational temperatures and geometry from atoms
        object
        
        Parameters
        ----------
        symmetrynumber : int
            Symmetry number.
        atoms : ase.Atoms object, optional
            An atoms object can be used to calculate rot_temperatures and
            guess geometry
        """
        symmetrynumber = symmetrynumber
        geometry = get_geometry_from_atoms(atoms=atoms)
        rot_temperatures = get_rot_temperatures_from_atoms(atoms=atoms, 
            geometry=geometry)

        idealrot = cls(symmetrynumber=symmetrynumber, 
            rot_temperatures=rot_temperatures, geometry=geometry)
        return idealrot

    def get_q(self, T):
        """Calculates the partition function

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
                *(T**3/np.prod(self.rot_temperatures))**0.5
        else:
            raise ValueError(
                'Geometry, {}, not supported.'.format(self.geometry))
        
    def get_CvoR(self):
        """Calculates the dimensionless heat capacity at constant volume

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

        Returns
        -------
            CpoR_rot : float
                Rotational dimensionless heat capacity at constant pressure
        """
        return self.get_CvoR()
    
    def get_UoRT(self):
        """Calculates the imensionless internal energy

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

        Returns
        -------
            HoRT_rot : float
                Rotational dimensionless enthalpy
        """
        return self.get_UoRT()

    def get_SoR(self, T):
        """Calculates the dimensionless entropy

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
            return np.log(T/self.symmetrynumber \
                /np.prod(self.rot_temperatures)) + 1.
        elif self.geometry == 'nonlinear':
            return np.log(np.sqrt(np.pi)/self.symmetrynumber* \
                (T**3/np.prod(self.rot_temperatures))**0.5) + 1.5
        else:
            raise ValueError(
                'Geometry, {}, not supported.'.format(self.geometry))

    def get_AoRT(self, T):
        """Calculates the dimensionless Helmholtz energy

        Parameters
        ----------
            T : float
                Temperature in K
        Returns
        -------
            AoRT_rot : float
                Rotational dimensionless Helmholtz energy
        """
        return self.get_UoRT()-self.get_SoR(T=T)

    def get_GoRT(self, T):
        """Calculates the dimensionless Gibbs energy

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

def get_rot_temperatures_from_atoms(atoms, geometry):
    """Calculate the rotational temperatures from ase.Atoms object

    Parameters
    ----------
        atoms : ase.Atoms object
            Atoms object
    Returns
    -------
        rot_temperatures : list of float
            Rotational temperatures
    """
    rot_temperatures = []
    for moment in atoms.get_moments_of_inertia():
        if np.isclose(0., moment):
            continue
        moment_SI = moment*c.convert_unit(from_='amu', to='kg') \
            *c.convert_unit(from_='A2', to='m2')
        rot_temperatures.append(
            c.h('eV s', bar=True)**2/2/c.kb('eV/K')/moment_SI \
                *c.convert_unit(from_='eV', to='J'))

    if geometry == 'monatomic':
        assert np.isclose(np.sum(rot_temperatures), 0.)
        return [0.]
    elif geometry == 'linear':
        assert np.isclose(rot_temperatures[0], rot_temperatures[1])
        return [max(rot_temperatures)]
    elif geometry == 'nonlinear':
        return rot_temperatures
    else:
        raise ValueError(
            'Geometry, {}, not supported.'.format(geometry))

def get_geometry_from_atoms(atoms, degree_tol=10.):
    """Estimate the geometry using the ase.Atoms object

    Parameters
    ----------
        atoms : ase.Atoms object
            Atoms object
        degree_tol : float
            Degree tolerance in degrees
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
            if np.isclose(angle, 0., atol=degree_tol) \
                or np.isclose(angle, 180., atol=degree_tol):
                return 'linear'
        else:
            return 'nonlinear'
