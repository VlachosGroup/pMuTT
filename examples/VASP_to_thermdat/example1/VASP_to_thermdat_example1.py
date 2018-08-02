import matplotlib.pyplot as plt
from pprint import pprint

from PyMuTT import constants as c
from PyMuTT.io_.excel import read_excel
from PyMuTT.io_.thermdat import write_thermdat
from PyMuTT.models.empirical import BaseThermo
from PyMuTT.models.empirical.nasa import Nasa
from PyMuTT.models.empirical.references import References

'''
User inputs
'''
#Reference information
refs_path = './references.xlsx'

#Input information
species_path = './thermdat_input.xlsx'
T_low = 200.
T_high = 1100. #K

#Output information
thermdat_path = './thermdat'

#Miscellaneous options
show_plot = True

'''
Processing References
'''
#Import from excel
refs_input = read_excel(io=refs_path)
refs = References([BaseThermo(**ref_input) for ref_input in refs_input])
print('Reference Input:')
pprint(refs_input)

'''
Processing Input Species
'''
#Import from excel
species_data = read_excel(io=species_path)
species = [Nasa(references=refs, T_low=T_low, T_high=T_high, T_ref=c.T0('K'), **specie_data) for specie_data in species_data]
print('Species Input:')
pprint(species_data)

'''
Printing Out Results
'''
write_thermdat(nasa_species=species, filename=thermdat_path)
if show_plot:
	for specie in species:
		specie.plot_thermo_model_and_empirical(Cp_units='J/mol/K', H_units='kJ/mol', S_units='J/mol/K', G_units='kJ/mol')
	plt.show()

