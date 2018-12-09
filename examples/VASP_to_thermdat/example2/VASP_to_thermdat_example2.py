import os
import matplotlib.pyplot as plt
from pMuTT import constants as c
from pMuTT.io_.excel import read_excel
from pMuTT.io_.thermdat import write_thermdat
from pMuTT.models.empirical.nasa import Nasa
from pMuTT.models.empirical.references import Reference, References

'''
User inputs
'''
base_path = os.path.dirname(__file__)

#Reference information
refs_path = '{}/references.xlsx'.format(base_path)

#Surface information
surfaces_path = '{}/surfaces.xlsx'.format(base_path)

#Input information
species_path = '{}/input_data.xlsx'.format(base_path)
T_low = 100.
T_high = 1500. #K

#Output information
thermdat_path = '{}/thermdat'.format(base_path)

#Miscellaneous options
save_plot = True

'''
Processing References
'''
#Import from excel
refs_data = read_excel(io=refs_path)
refs = References(references=[Reference(**ref_data) for ref_data in refs_data])

'''
Processing Surfaces
'''
surfaces_data = read_excel(io=surfaces_path)

'''
Processing Input Species
'''
#Import from excel
species_data = read_excel(io=species_path)
#Adjust potential energy to subtract contribution from surface
for specie_data in species_data:
    #Find appropriate surface
    for surface_data in surfaces_data:
        if surface_data['name'] in specie_data['notes']:
            specie_data['potentialenergy'] -= surface_data['potentialenergy']
            break
species = [Nasa.from_statmech(references=refs, T_low=T_low, T_high=T_high, \
                              **specie_data) for specie_data in species_data]

'''
Printing Out Results
'''
write_thermdat(nasa_species=species, filename=thermdat_path)
if save_plot:
    for nasa in species:
        nasa.plot_statmech_and_empirical(Cp_units='J/mol/K', H_units='kJ/mol', 
                                         S_units='J/mol/K', G_units='kJ/mol')
        thermo_plot_path = os.path.join( \
                *[base_path, 'thermo_plots', '{}.png'.format(nasa.name)])
        plt.savefig(thermo_plot_path)