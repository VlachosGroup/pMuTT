from pMuTT.reaction import Reaction
from pMuTT.reaction.bep import BEP
from pMuTT.statmech import StatMech, presets
from pMuTT import constants as c
from pMuTT.io_.jsonio import pMuTTEncoder, json_to_pMuTT
import numpy as np
import os
import json
from pprint import pprint

os.chdir(os.path.dirname(__file__))
T = c.T0('K')
# Factor to convert potential energy to dimensionless number
dim_factor = c.R('eV/K')*T

species = {
    'H2': StatMech(name='H2', potentialenergy=2.*dim_factor, 
                   **presets['electronic']),
    'O2': StatMech(name='O2', potentialenergy=4.*dim_factor, 
                   **presets['electronic']),
    'H2O': StatMech(name='H2O', potentialenergy=3.*dim_factor, 
                    **presets['electronic']),
    'H2O_TS': StatMech(name='H2O', potentialenergy=5.*dim_factor, 
                       **presets['electronic']),
}

'''Initializing the reaction with the BEP relationship'''
rxn = Reaction.from_string(reaction_str='H2+0.5O2=H2O_TS=H2O', species=species,
                           descriptor='delta', slope=10., intercept=20.,
                           bep=BEP)
print(rxn)
print('Delta H/RT = {}'.format(rxn.get_delta_HoRT(T=T)))
print('BEP Forward E/RT = {}'.format(
        rxn.get_EoRT_act(T=T, rev=False, method='bep')))
print('BEP Reverse E/RT = {}'.format(
        rxn.get_EoRT_act(T=T, rev=True, method='bep')))
print('TS Forward E/RT = {}'.format(
        rxn.get_EoRT_act(T=T, rev=True, method='ts')))
print('"Any" E/RT (Expecting the same answer as TS)= {}'.format(
        rxn.get_EoRT_act(T=T, rev=True, method='any')))

'''Writing and reading from JSON file'''
with open('bep.json', 'w') as f_ptr:
    json.dump(rxn, f_ptr, cls=pMuTTEncoder, indent=True)

with open('bep.json', 'r') as f_ptr:
    rxn_io = json.load(f_ptr, object_hook=json_to_pMuTT)

print('-')
print('Delta H/RT from JSON file = {}'.format(rxn_io.get_delta_HoRT(T=T)))
print('BEP Forward E/RT from JSON file = {}'.format(
        rxn.get_EoRT_act(T=T, rev=False, method='bep')))