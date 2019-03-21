import os
import json

import numpy as np
from matplotlib import pyplot as plt

from pMuTT import pMuTT_list_to_dict
from pMuTT.reaction import Reaction
from pMuTT.io_.thermdat import read_thermdat
from pMuTT.io_.excel import read_excel
from pMuTT.io_.json import pMuTTEncoder

os.chdir(os.path.dirname(__file__))

'''Read species'''
nasa_species = read_thermdat('thermdat')
nasa_dict = pMuTT_list_to_dict(nasa_species)

'''Read reactions'''
rxns_data = read_excel('reactions.xlsx')
rxns = [Reaction.from_string(species=nasa_dict, **rxn_data) for rxn_data in rxns_data]

'''Save reactions as json'''
with open('reactions.json', 'w') as f_ptr:
    json.dump(rxns, f_ptr, cls=pMuTTEncoder, indent=True)

'''
Get thermodynamics from Reaction
'''
T = np.linspace(200., 3500.)
for rxn in rxns:
    HoRT_rxn = rxn.get_delta_HoRT(T=T)
    SoR_rxn = rxn.get_delta_SoR(T=T)
    GoRT_rxn = rxn.get_delta_GoRT(T=T)

    '''
    Plot the data
    '''
    f, ax = plt.subplots(3, sharex=True)
    ax[0].set_title(rxn.to_str())
    ax[0].plot(T, HoRT_rxn, 'r-')
    ax[0].set_ylabel('H/RT')

    ax[1].plot(T, SoR_rxn, 'g-')
    ax[1].set_ylabel('S/R')

    ax[2].plot(T, GoRT_rxn, 'b-')
    ax[2].set_ylabel('G/RT')
    ax[2].set_xlabel('Temperature (K)')
plt.show()