import os
import numpy as np
from matplotlib import pyplot as plt
from pMuTT.io_.thermdat import read_thermdat
from pMuTT import pMuTT_list_to_dict
from pMuTT import reaction as rxn

'''
Read the thermdat
'''
base_path = os.path.dirname(__file__)
nasa_list = read_thermdat(os.path.join(base_path, 'thermdat'))
nasa_dict = pMuTT_list_to_dict(nasa_list, key='name')

'''
Create the reaction object
'''
reaction_str = 'H2+0.5O2=H2O'
rxn_nasa = rxn.Reaction.from_string(reaction_str=reaction_str,
                                    species=nasa_dict)

'''
Check that the equation is balanced
'''
rxn_nasa.check_element_balance()

'''
Get thermodynamics from Reaction
'''
T = np.linspace(200., 3500.)
HoRT_rxn = rxn_nasa.get_delta_HoRT(T=T)
SoR_rxn = rxn_nasa.get_delta_SoR(T=T)
GoRT_rxn = rxn_nasa.get_delta_GoRT(T=T)

'''
Plot the data
'''
f, ax = plt.subplots(3, sharex=True)
ax[0].set_title(rxn_nasa.to_str())
ax[0].plot(T, HoRT_rxn, 'r-')
ax[0].set_ylabel('H/RT')

ax[1].plot(T, SoR_rxn, 'g-')
ax[1].set_ylabel('S/R')

ax[2].plot(T, GoRT_rxn, 'b-')
ax[2].set_ylabel('G/RT')
ax[2].set_xlabel('Temperature (K)')
plt.show()
