# -*- coding: utf-8 -*-
'''
         -----------------------------------------------------
          Calculate thermodynamic data (S298, H298, and Cp(T)
          from ab initio DFT data (energies and frequencies)
              providng input thermodynamics files for
             KMC (Zacros) and MKM (Chemkin and Matlab)

                     Vlachos Research Group
                Chemical and Biomolecular Egineering
                      University of Delaware

                         Jonathan A Lym
                     Gerhard R Wittreich, P.E.
          -----------------------------------------------------

Created on Fri Mar 31 16:20:43 2017

@author: wittregr

Adopted from Matlab code written and modified by:

                     Vassili Vorotnikov
                            and
                         Geun Ho Gu

 This program contains the class objects used to read energy, vibration and
 molecular configuration data and determine the standard entropy and enthalpy
 and heat capacities at various temperatures.

'''
import numpy as _np
import os as _os
import ase.io as _ase
import ase.thermochemistry as _thermo
import re as _re
import scipy.interpolate as _sp
import datetime as _datetime
from GRW_constants import constant as _c


class Particle(object):

    Cp_Range = _np.linspace(100, 1500, 15)
    VibScalingFactor = 1         # Vibrational Scaling Factor

    def __init__(self, data, dict, Base_path, Tstp=298.15, Pstp=100000.):

        '''
        Fill object with species name and associated thermodynamic data
        '''
        self.name = str(data[dict['name']])         # Species name
        if 'numh' in dict:
            if data[dict['numc']] != '':
                self.carbon = int(data[dict['numc']])       # No. of C atoms
            else:
                self.carbon = int(0)
            if data[dict['numh']] != '':
                self.hydrogen = int(data[dict['numh']])     # No. of H atoms
            else:
                self.hydrogen = int(0)
            if data[dict['numo']] != '':
                self.oxygen = int(data[dict['numo']])       # No. of O atoms
            else:
                self.oxygen = int(0)
            if data[dict['numn']] != '':
                self.nitrogen = int(data[dict['numn']])     # No. of N atoms
            else:
                self.nitrogen = int(0)
        if 'mw' in dict:
            self.MW = float(data[dict['mw']])/_c.NA/1000
        else:
            self.MW = (self.carbon * _c.MW_carbon +
                       self.hydrogen * _c.MW_hydorgen +
                       self.oxygen * _c.MW_oxygen +
                       self.nitrogen * _c.MW_nitrogen)/_c.NA/1000
        if not hasattr(self, 'Inertia') and 'totengpath' in dict:
            self.totengpath = _os.path.join(*_re.split(r'\\|/',
                                            str(data[dict['totengpath']]).
                                            strip('.').strip('\\')))
        if 'etotal' in dict:
            self.etotal = float(data[dict['etotal']])   # Total energy
        if 'hf_tstp' in dict:
            self.hf_Tstp = float(data[dict['hf_tstp']])
        if 'edisp' in dict and data[dict['edisp']] != '':
            self.edisp = float(data[dict['edisp']])     # Dispersion energy
        else:
            self.edisp = float(0)
        if 'spin' in dict:
            self.spin = float(data[dict['spin']])   # Spin
        else:
            self.spin = float(0)
        self.numvibfreq = int(data[dict['numvibfreq']])  # No. vib frequencies
        if not hasattr(self, 'phase'):
            self.phase = None                       # Phase (G=gas, S=surface)
        self.vibfreq = []                           # Vibration frequencies
        for x in range(0, self.numvibfreq):
            self.vibfreq.append(Particle.VibScalingFactor *
                                float(data[dict['vibfreq'] + x]))
        self.vibfreq = _np.array(self.vibfreq)
        self.Base_path = Base_path
        self.Tstp = Tstp
        self.Pstp = Pstp
        self.ThermoProperties()

    def ThermoProperties(self):
        '''
        Calculate all thermodynamic properties from input data
        '''

        '''
        Get rotational data from VASP CONTCAR for gas species using
        Atomic Simulation Environment (ASE) libraries for python
        '''
        if self.phase == 'G':
            if hasattr(self, 'Inertia'):
                self.I3 = self.Inertia
            else:
                filepath = _os.path.join(self.Base_path,
                                         self.totengpath,
                                         'CONTCAR')
                self.VASP = _ase.read(filepath)
                self.I3 = self.VASP.get_moments_of_inertia() *\
                    _c.A2_to_m2*_c.amu_to_kg
                self.MW = sum(self.VASP.get_masses())/_c.NA/1000.
            self.T_I = _c.h1**2/(8*_np.pi**2*_c.kb1)
        '''
        Calulcate common frequency data for vibrational components
        '''
        self.nu = self.vibfreq * 100 * _c.c2
        self.theta = _c.h1 * self.nu / _c.kb1
        '''
        Initialize ASE thermochemistry modules
        '''
        self.vib_energy = self.nu * _c.h3  # Energies in eV
        if self.phase == "G":
            self.ASE = _thermo.IdealGasThermo(vib_energies=self.vib_energy,
                                              geometry=self.Geometry,
                                              potentialenergy=self.etotal,
                                              spin=self.spin,
                                              atoms=self.VASP,
                                              symmetrynumber=self.sigma)

        else:
            self.ASE = _thermo.HarmonicThermo(vib_energies=self.vib_energy,
                                              potentialenergy=self.etotal)

        '''
        Call Entropy method to calculate standard state entropy
        '''
        self.Calc_Entropy()

        '''
        Call Heat Capacity method to calculate heat capacities at the
        temperature range specified by Cp_Range
        '''
        self.Calc_HeatCapacities()

        '''
        Call Enthalpy method to calculate standard state enthalpy
        '''
        if not hasattr(self, 'hf_Tstp'):
            self.Calc_Enthalpy()

    def Calc_HeatCapacities(self):
        '''
        Calculate the vibrational, rotational and translational components and
        total heat capacity of a species for a range of temperatures
        '''
        '''
        Calculate vibrational contribution to heat capacity for temperature
        range specified in Cp_Range for linear and non-linear species
        '''
        vib = _np.zeros([len(self.theta), len(self.Cp_Range)])
        for x in range(0, len(self.Cp_Range)):
            vib[:, x] = _np.divide(self.theta, self.Cp_Range[x])
        self.Cp_vib = sum(_np.divide(_np.multiply(_np.power(vib, 2),
                                                  _np.exp(-vib)),
                                     _np.power(1-_np.exp(-vib), 2)))
        if self.phase == 'G':
            '''
            Translational and rotational Gas phase calculation
            '''
            self.Cp_trans = 3./2
            self.Cv2Cp = 1.
            if self.islinear == 0:
                '''
                Non-Linear species
                '''
                self.Cp_rot = 3./2
            else:
                '''
                Linear species
                '''
                self.Cp_rot = 1.
        else:
            '''
            Surface species
            '''
            self.Cp_rot = 0.
            self.Cp_trans = 0.
            self.Cv2Cp = 0.
        '''
        Sum all contribution to heat capacity for total heat capapcity
        '''
        self.Cp = _c.R1*(self.Cp_trans + self.Cp_rot +
                         self.Cp_vib + self.Cv2Cp)

    def Calc_Enthalpy(self):

        if self.phase == 'G':
            self.dfth = self.ASE.get_enthalpy(
                temperature=self.Tstp,
                verbose=False)/_c.R3*_c.R4
        else:
            self.dfth = self.ASE.get_internal_energy(
                temperature=self.Tstp,
                verbose=False)/_c.R3*_c.R4

    def Calc_Entropy(self):

        if self.phase == "G":
            self.S_Tstp = self.ASE.get_entropy(
                temperature=self.Tstp,
                pressure=self.Pstp,
                verbose=False)/_c.R3*_c.R1
        else:
            self.S_Tstp = self.ASE.get_entropy(
                temperature=self.Tstp,
                verbose=False)/_c.R3*_c.R1


class Reference(Particle):
    '''
    SubClass object to add specific fields for reference species
    '''
    def __init__(self, data, dict, Base_path, Tstp=298.15):
        if 'sigma' in dict and data[dict['sigma']] != '':
            self.sigma = int(data[dict['sigma']])            # Sigma
        else:
            self.sigma = int(0)
        if 'islinear' in dict and data[dict['islinear']] != '':
            self.islinear = int(data[dict['islinear']])   # Is molecule linear?
            Geometry = ['nonlinear', 'linear']
            self.Geometry = Geometry[self.islinear]
        else:
            self.linear = int(-1)
        if 'hf298nist' in dict:
            self.hf298nist = float(data[dict['hf298nist']])  # NIST Std enthpy
        if 'inertia' in dict:
            if data[dict['inertia']] != '':
                self.Inertia = float(data[dict['inertia']])
            else:
                self.Inertia = float(0)
        self.phase = str.upper(data[dict['phase']])      # Phase
        if 'a_st' in dict and data[dict['a_st']] != '':
            self.A_st = float(data[dict['a_st']])
        super(Reference, self).__init__(data, dict, Base_path, Tstp)

    @staticmethod
    def BasisSet(RefSpecies):
        A = []
        b_nist = []
        b_dfth = []
        for x in range(0, len(RefSpecies)):
            A.append([RefSpecies[x].carbon, RefSpecies[x].hydrogen,
                      RefSpecies[x].oxygen, RefSpecies[x].nitrogen])
            b_nist.append([RefSpecies[x].hf298nist])
            b_dfth.append([RefSpecies[x].dfth])
        ref = _np.linalg.lstsq(A, _np.array(b_nist) - _np.array(b_dfth),
                               rcond=None)[0]
        return(ref)


class Target(Particle):
    '''
    SubClass object to add specific fields for target surface species
    '''
    def __init__(self, data, dict, Base_path, Tstp=298.15):
        self.surface = str(data[dict['surface']])          # Surface
        self.functional = str(data[dict['functional']])    # Functional
        self.kpoints = str(data[dict['kpoints']])          # k-Points
        self.vibfreqpath = str(data[dict['vibfreqpath']])  # Unused
        self.phase = None                                  # Phase
        super(Target, self).__init__(data, dict, Base_path, Tstp)

    @staticmethod
    def ReferenceDFT(Species, Surface, Basis):
        for x in range(0, len(Species)):
            Molecule = _np.array([Species[x].carbon, Species[x].hydrogen,
                                  Species[x].oxygen, Species[x].nitrogen])
            if Species[x].phase == 'G':
                Species[x].hf_Tstp = (Species[x].dfth +
                                      _np.dot(Molecule, Basis))[0]
                if hasattr(Species[x], 'edisp'):
                    Species[x].convedisp = (Species[x].edisp *
                                            _c.ev_atom_2_kcal_mol)
            else:
                Slab = next((y for y in Surface if y.name ==
                             Species[x].surface),
                            None)
                if Slab is None:
                    print('Error')
                else:
                    Species[x].hf_Tstp = (Species[x].dfth +
                                          _np.dot(Molecule, Basis) -
                                          Slab.etotal *
                                          _c.ev_atom_2_kcal_mol)[0]
                    if hasattr(Species[x], 'edisp'):
                        Species[x].convedisp = (Species[x].edisp *
                                                _c.ev_atom_2_kcal_mol -
                                                Slab.edisp *
                                                _c.ev_atom_2_kcal_mol)
        return(Species)

    @staticmethod
    def CreateThermdat(Species, Base_path, Output):
        '''
        Calculate the seven coefficients for the NASA polynomials for two
        temperature ranges and output the results in a Chemkin format thermdat
        file
        '''
        T_mid = 500
        Tstp = Species[0].Tstp

        def HS_NASA(T, a):
            '''
            7-coefficient NASA polynomials for enthalpy and entropy
            '''
            Enthalpy = a[0] + a[1]*T/2 + a[2]*T**2/3 + \
                a[3]*T**3/4 + a[4]*T**4/5
            Entropy = a[0]*_np.log(T) + a[1]*T + a[2]*T**2/2 + \
                a[3]*T**3/3 + a[4]*T**4/4
            return[Enthalpy, Entropy]

        for x in range(0, len(Species)):
            T_rng_low = _np.linspace(min(Species[x].Cp_Range), T_mid, 1600)
            T_rng_high = _np.linspace(T_mid, max(Species[x].Cp_Range), 4000)
            T_func = _sp.InterpolatedUnivariateSpline(Species[x].Cp_Range,
                                                      Species[x].Cp/_c.R1, k=4)
            '''
            Fit coefficients A1-A5 to heat capacity data
            '''
            Species[x].a_low = _np.polyfit(T_rng_low,
                                           T_func(T_rng_low), 4)[::-1]
            Species[x].a_high = _np.polyfit(T_rng_high,
                                            T_func(T_rng_high), 4)[::-1]
            '''
            Correct A1 high temperature range coefficient to eliminate
            discontinuity between high and low temperature range polynomials
            '''
            Species[x].a_high[0] = Species[x].a_high[0] + \
                (_np.polyval(Species[x].a_low[::-1], T_mid) -
                 _np.polyval(Species[x].a_high[::-1], T_mid))

            '''
            Determine A6 coefficient for enthalpy calculations
            '''
            a6_high = (Species[x].hf_Tstp/_c.R1/Tstp*1000 -
                       HS_NASA(Tstp, Species[x].a_high)[0])*Tstp
            a6_low = (Species[x].hf_Tstp/_c.R1/Tstp*1000 -
                      HS_NASA(Tstp, Species[x].a_low)[0])*Tstp
            '''
            Correct A6 high temperature range coefficient to eliminate
            discontinuity between high and low temperature range polynomials
            '''
            a6_high_delta = (HS_NASA(T_mid, Species[x].a_low)[0] +
                             a6_low/T_mid) - \
                            (HS_NASA(T_mid,
                             Species[x].a_high)[0] + a6_high/T_mid)
            a6_high = a6_high + a6_high_delta * T_mid
            Species[x].a_high = _np.append(Species[x].a_high, a6_high)
            Species[x].a_low = _np.append(Species[x].a_low, a6_low)

            '''
            Determine A7 coefficient for entropy calculations
            '''
            a7_high = Species[x].S_Tstp/_c.R1 - \
                HS_NASA(Tstp, Species[x].a_high)[1]
            a7_low = Species[x].S_Tstp/_c.R1 - \
                HS_NASA(Tstp, Species[x].a_low)[1]

            '''
            Correct A7 high temperature range coefficient to eliminate
            discontinuity between high and low temperature range polynomials
            '''
            a7_high_delta = (HS_NASA(T_mid, Species[x].a_low)[1] +
                             a7_low) - (HS_NASA(T_mid,
                                        Species[x].a_high)[1] + a7_high)
            a7_high = a7_high + a7_high_delta
            Species[x].a_high = _np.append(Species[x].a_high, a7_high)
            Species[x].a_low = _np.append(Species[x].a_low, a7_low)

        '''
        Write the species name, seven NASA coefficients for both a high and
        a low temperature range and other data in the Chemkin thermdat
        file format
        '''
        if _os.path.isdir(_os.path.join(Base_path, Output)) is False:
            _os.mkdir(_os.path.join(Base_path, Output))
        filepath = _os.path.join(Base_path, Output, 'thermdat')
        fid = open(filepath, 'w')
        fid.truncate()
        '''
        Write thermdat file header
        '''
        fid.write('THERMO ALL\n')
        fid.write('%10.0f%10.0f%10.0f\n' % (min(Species[x].Cp_Range),
                                            T_mid,
                                            max(Species[x].Cp_Range)))
        for s in range(0, _np.size(Species)):
            '''
            Write header line for each species on line 1
            '''
            fid.write('%-16s' % (Species[s].name))
            fid.write('%-8s' % (_datetime.date.today().strftime('%Y%m%d')))
            fid.write('%1s%4i' % ('C', Species[s].carbon))
            fid.write('%1s%4i' % ('H', Species[s].hydrogen))
            fid.write('%1s%4i' % ('O', Species[s].oxygen))
            fid.write('%1s%4i' % ('N', Species[s].nitrogen))
            if Species[s].name.find('(S)'):
                fid.write('S')
            else:
                fid.write('G')
            fid.write('%10.1f%10.1f%8.1f' % (min(Species[x].Cp_Range),
                                             max(Species[x].Cp_Range),
                                             T_mid))
            fid.write('%6s%1i\n' % ('', 1))
            '''
            Write first five NASA coefficients for
            low temperature range on line 2
            '''
            for x in range(0, 5):
                fid.write('%15E' % (Species[s].a_high[x]))
            fid.write('%4s%1i\n' % ('', 2))
            '''
            Write final two NASA coefficients for
            low temperature range on line 2
            '''
            for x in range(0, 2):
                fid.write('%15E' % (Species[s].a_high[x+5]))
            '''
            Write first three NASA coeficients for
            high temperature range on line 3
            '''
            for x in range(0, 3):
                fid.write('%15E' % (Species[s].a_low[x]))
            fid.write('%4s%1i\n' % ('', 3))
            '''
            Write final four NASA coefficients for
            high temperature range on line 4
            '''
            for x in range(0, 4):
                fid.write('%15E' % (Species[s].a_low[x+3]))
            fid.write('%19s%1i\n' % ('', 4))
        '''
        Write file footer and close the file
        '''
        fid.write('END\n')
        fid.close()

        return(Species)


class Surface:
    '''
    Class object to populate slab energies for surfaces
    '''
    def __init__(self, data, dict):
        self.name = str(data[dict['name']])          # Surface name
        self.etotal = float(data[dict['etotal']])    # Total energy-DFT
        self.edisp = float(data[dict['edisp']])      # Dispersion energy-DFT


def DFTFileRead(filepath):
    fid = open(filepath, 'r')
    file = fid.read()
    lines = file.splitlines()
    dict_array = lines[2].lower().split('\t')
    dict = {}
    for x in range(0, len(dict_array)):
        dict[dict_array[x]] = x
    return(lines, dict)
