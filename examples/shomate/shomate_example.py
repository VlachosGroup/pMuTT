import numpy as np
from matplotlib import pyplot as plt
from pMuTT.empirical.shomate import Shomate

T_low = 298.
T_high = 1000.
T_ref = 298.
T = np.linspace(T_low, T_high)

'''Initializing a Shomate polynomial directly'''
shomate_H2 = Shomate(name='H2', T_low=T_low, T_high=T_high,
                     a=np.array([33.066178, -11.363417, 11.432816, -2.772874,
                                 -0.158558, -9.980797, 172.707974, 0.]))
fig1, ax1 = shomate_H2.plot_empirical()
ax1[0].set_title('Directly inputting polynomial')

'''Using shomate to calculate thermodynamic properties'''
CpoR = shomate_H2.get_CpoR(T=T)
HoRT_ref = shomate_H2.get_HoRT(T=T_ref)
SoR_ref = shomate_H2.get_SoR(T=T_ref)

'''Initializing a Shomate using data'''
shomate_H2_data = Shomate.from_data(name='H2', T=T, CpoR=CpoR, T_ref=T_ref,
                                    HoRT_ref=HoRT_ref, SoR_ref=SoR_ref)
fig2, ax2 = shomate_H2_data.plot_empirical()
ax2[0].set_title('Fitting to data')
plt.show()
