import matplotlib.pyplot as plt
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

#Surface information
surfaces_path = './surfaces.xlsx'

#Input information
species_path = './input_data.xlsx'
T_low = 100.
T_high = 1500. #K

#Output information
thermdat_path = './thermdat'

#Miscellaneous options
show_plot = True

'''
Processing References
'''
#Import from excel
refs_data = read_excel(io=refs_path)
refs = References([BaseThermo(**ref_data) for ref_data in refs_data])

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
species = [Nasa(references=refs, T_low=T_low, T_high=T_high, T_ref=c.T0('K'), **specie_data) for specie_data in species_data]

'''
Printing Out Results
'''
write_thermdat(nasa_species=species, filename=thermdat_path)
if show_plot:
	for nasa in species:
		nasa.plot_thermo_model_and_empirical(Cp_units='J/mol/K', H_units='kJ/mol', S_units='J/mol/K', G_units='kJ/mol')
	plt.show()