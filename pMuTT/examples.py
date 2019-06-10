# -*- coding: utf-8 -*-

import numpy as np
from pMuTT.statmech import StatMech, presets
from pMuTT.empirical.nasa import Nasa
from pMuTT.empirical.shomate import Shomate
from pMuTT.reaction import Reaction, Reactions

O2_nasa = Nasa(name='O2', T_low=200., T_mid=1000., T_high=3500.,
               elements={'O': 2},
               a_low=np.array([3.78245636E+00, -2.99673416E-03, 9.84730201E-06,
                               -9.68129509E-09, 3.24372837E-12, -1.06394356E+03,
                               3.65767573E+00]),
               a_high=np.array([3.28253784E+00, 1.48308754E-03, -7.57966669E-07,
                                2.09470555E-10, -2.16717794E-14,
                                -1.08845772E+03, 5.45323129E+00]))

H2_shomate = Shomate(name='H2', T_low=298., T_high=1000., 
                     a=np.array([33.066178, -11.363417, 11.432816, -2.772874,
                                 -0.158558, -9.980797, 172.707974, 0.]))



# H2O_statmech = StatMech(**presets['idealgas'])