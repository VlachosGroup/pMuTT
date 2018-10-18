import os
import json
import numpy as np
from pprint import pprint
from matplotlib import pyplot as plt
from pMuTT.io_.thermdat import read_thermdat
from pMuTT.io_.jsonio import pMuTTEncoder, json_to_pMuTT
from pMuTT.models.empirical.nasa import Nasa

base_path = os.path.dirname(__file__)
thermdat_path = os.path.join(base_path, 'thermdat')
json_path = os.path.join(base_path, 'nasa_obj.json')

'''
Loading thermdat file
'''
#Thermdat file from http://combustion.berkeley.edu/gri_mech/version30/files30/thermo30.dat
species = read_thermdat(thermdat_path)

'''
Saving and Loading from JSON
'''
with open(json_path, 'w') as f_ptr:
    json.dump(species, f_ptr, cls=pMuTTEncoder, indent=True)

with open(json_path, 'r') as f_ptr:
    species_json = json.load(f_ptr, object_hook=json_to_pMuTT)

'''
Plotting to ensure the object was written and read correctly
'''
for specie, specie_json in zip(species, species_json):
    T = np.linspace(specie.T_low, specie.T_high)
    plt.figure()
    f, ax = plt.subplots(4, sharex=True)    
    # CpoR Plot
    ax[0].plot(T, specie.get_CpoR(T=T), 'r-', label='Thermdat')
    ax[0].plot(T, specie_json.get_CpoR(T=T), 'b*', label='JSON')
    ax[0].legend()
    ax[0].set_ylabel('Cp/R')

    # HoRT Plot
    ax[1].plot(T, specie.get_HoRT(T=T), 'r-')
    ax[1].plot(T, specie_json.get_HoRT(T=T), 'b*')
    ax[1].set_ylabel('H/RT')

    # SoR Plot
    ax[2].plot(T, specie.get_SoR(T=T), 'r-')
    ax[2].plot(T, specie_json.get_SoR(T=T), 'b*')
    ax[2].set_ylabel('S/R')

    # HoRT Plot
    ax[3].plot(T, specie.get_GoRT(T=T), 'r-')
    ax[3].plot(T, specie_json.get_GoRT(T=T), 'b*')
    ax[3].set_ylabel('G/RT')
    ax[3].set_xlabel('T (K)')
    f.savefig(os.path.join(base_path, '{}.png'.format(specie.name)))
