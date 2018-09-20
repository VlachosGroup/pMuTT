import os
import json
import numpy as np
from pprint import pprint
from matplotlib import pyplot as plt
from PyMuTT.io_.thermdat import read_thermdat
from PyMuTT.io_.jsonio import PyMuTTEncoder, json_to_PyMuTT
from PyMuTT.models.empirical.nasa import Nasa

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
    json.dump(species, f_ptr, cls=PyMuTTEncoder, indent=True)

with open(json_path, 'r') as f_ptr:
    species_json = json.load(f_ptr, object_hook=json_to_PyMuTT)

'''
Plotting to ensure the object was written and read correctly
'''
for specie, specie_json in zip(species, species_json):
    Ts = np.linspace(specie.T_low, specie.T_high)
    plt.figure()
    f, ax = plt.subplots(4, sharex=True)    
    # CpoR Plot
    ax[0].plot(Ts, specie.get_CpoR(Ts=Ts), 'r-', label='Thermdat')
    ax[0].plot(Ts, specie_json.get_CpoR(Ts=Ts), 'b*', label='JSON')
    ax[0].legend()
    ax[0].set_ylabel('Cp/R')

    # HoRT Plot
    ax[1].plot(Ts, specie.get_HoRT(Ts=Ts), 'r-')
    ax[1].plot(Ts, specie_json.get_HoRT(Ts=Ts), 'b*')
    ax[1].set_ylabel('H/RT')

    # SoR Plot
    ax[2].plot(Ts, specie.get_SoR(Ts=Ts), 'r-')
    ax[2].plot(Ts, specie_json.get_SoR(Ts=Ts), 'b*')
    ax[2].set_ylabel('S/R')

    # HoRT Plot
    ax[3].plot(Ts, specie.get_GoRT(Ts=Ts), 'r-')
    ax[3].plot(Ts, specie_json.get_GoRT(Ts=Ts), 'b*')
    ax[3].set_ylabel('G/RT')
    ax[3].set_xlabel('T (K)')
    f.savefig(os.path.join(base_path, '{}.png'.format(specie.name)))
