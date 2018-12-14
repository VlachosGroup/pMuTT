import numpy as np
import matplotlib.pyplot as plt
from pMuTT import constants as c
from pMuTT.empirical.nasa import Nasa

# Gas phase heat capacity data (in J/mol/K) for CH3OH from NIST
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Units=SI&Mask=1#Thermo-Gas
T = np.array([50, 100, 150, 200, 273.15, 298.15, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2250, 2500, 2750, 3000])
Cp = np.array([34.00, 36.95, 38.64, 39.71, 42.59, 44.06, 44.17, 51.63, 59.7, 67.19, 73.86, 79.76, 84.95, 89.54, 93.57, 97.12, 100.24, 102.98, 105.4, 110.2, 113.8, 116.5, 118.6, 120, 121])
CpoR = Cp/c.R('J/mol/K')

#Enthalpy of Formation for CH3OH (in kJ/mol) from NIST
T_ref = c.T0('K')
H_ref = -205.
HoRT_ref = H_ref/c.R('kJ/mol/K')/T_ref
#Standard molar entropy (in J/mol/K) from Wikipedia, https://en.wikipedia.org/wiki/Methanol_(data_page)
S_ref = 239.9
SoR_ref = S_ref/c.R('J/mol/K')

#Units to plot the figure
Cp_units = 'J/mol/K'
H_units = 'kJ/mol'
S_units = 'J/mol/K'
G_units = 'kJ/mol'

#Input the experimental data and fitting to a NASA polynomial
CH3OH_nasa = Nasa.from_data(name ='CH3OH', T=T, CpoR=CpoR, T_ref=T_ref, HoRT_ref=HoRT_ref, SoR_ref=SoR_ref)

#Compare the Nasa polynomial to the input data
fig, axes = CH3OH_nasa.plot_empirical(Cp_units=Cp_units, H_units=H_units, S_units=S_units, G_units=G_units)
axes[0].plot(T, Cp, 'ko')
axes[1].plot(T_ref, H_ref, 'ko')
axes[2].plot(T_ref, S_ref, 'ko')
axes[3].plot(T_ref, H_ref - T_ref * S_ref * c.convert_unit(from_='J', to='kJ'), 'ko')
plt.show()