# -*- coding: utf-8 -*-

import numpy as np
from ase.build import molecule
from pmutt.statmech import StatMech, presets
from pmutt.empirical.nasa import Nasa
from pmutt.empirical.shomate import Shomate
from pmutt.empirical.references import Reference, References
from pmutt.reaction import Reaction, Reactions

O2_nasa = Nasa(name='O2', T_low=200., T_mid=1000., T_high=3500.,
               elements={'O': 2}, phase='G',
               a_low=np.array([3.78245636E+00, -2.99673416E-03, 9.84730201E-06,
                               -9.68129509E-09, 3.24372837E-12, -1.06394356E+03,
                               3.65767573E+00]),
               a_high=np.array([3.28253784E+00, 1.48308754E-03, -7.57966669E-07,
                                2.09470555E-10, -2.16717794E-14,
                                -1.08845772E+03, 5.45323129E+00]))

H2_shomate = Shomate(name='H2', T_low=298., T_high=1000., phase='G',
                     elements={'H':2}, 
                     a=np.array([33.066178, -11.363417, 11.432816, -2.772874,
                                 -0.158558, -9.980797, 172.707974, 0.]))

H2_ref = Reference(name='H2', elements={'H': 2}, potentialenergy=-6.7598,
                   atoms=molecule('H2'), symmetrynumber=2, spin=0.,
                   vib_wavenumbers=[4306.1793], T_ref=298., HoRT_ref=0.,
                   **presets['idealgas'])

H2O_ref = Reference(name='H2O', elements={'H': 2, 'O': 1},
                    potentialenergy=-14.2209, symmetrynumber=2, spin=0.,
                    atoms=molecule('H2O'),
                    vib_wavenumbers=[3825.434, 3710.2642, 1582.432],
                    T_ref=298., HoRT_ref=-97.60604334, **presets['idealgas'])

refs = References(references=[H2_ref, H2O_ref])

H2O_statmech = StatMech(name='H2O', phase='G',
                        atoms=molecule('H2O'), potentialenergy=-14.2209,
                        symmetrynumber=2, spin=0, references=refs,
                        vib_wavenumbers=[3825.434, 3710.264, 1582.432],
                        **presets['idealgas'])

H2O_TS_statmech = StatMech(name='H2O', elements={'H': 2, 'O': 1}, phase='G',
                        atoms=molecule('H2O'), potentialenergy=0.,
                        symmetrynumber=2, spin=0, references=refs,
                        vib_wavenumbers=[3825.434, 3710.264, 1582.432],
                        **presets['idealgas'])

species = {
    'H2': H2_shomate,
    'H2O': H2O_statmech,
    'O2': O2_nasa,
    'H2O_TS': H2O_TS_statmech,
}

rxn = Reaction.from_string('H2 + 0.5O2 = H2O_TS = H2O', species)