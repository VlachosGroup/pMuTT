import matplotlib.pyplot as plt
from Thermochemistry import constants as c
from Thermochemistry.io_.excel import read_excel
from Thermochemistry.io_.thermdat import write_thermdat
from Thermochemistry.models.empirical import BaseThermo
from Thermochemistry.models.empirical.thermdat import Thermdat
from Thermochemistry.models.empirical.references import References

'''
User inputs
'''
#Reference information
refs_path = './references.xlsx'

#Input information
thermdats_in_path = './thermdat_input.xlsx'
thermdats_out_path = './thermdat'
T_low = 500.
T_high = 1100. #K

#Miscellaneous options
show_plot = True

'''
Processing References
'''
#Import from excel
refs_input = read_excel(io=refs_path)
refs = References([BaseThermo(**ref_input) for ref_input in refs_input])
refs.calc_offset()

'''
Processing Input Species
'''
#Import from excel
thermdats_data = read_excel(io=thermdats_in_path)
thermdats = [Thermdat(**thermdat_data) for thermdat_data in thermdats_data]
for thermdat in thermdats:
	thermdat.T_ref = c.T0('K')
	thermdat.HoRT_ref = thermdat.thermo_model.get_HoRT(Ts=c.T0('K')) + refs.get_specie_offset(thermdat.elements)
	thermdat.calc_nasa(T_low=T_low, T_high=T_high, T_ref=c.T0('K'))

'''
Printing Out Results
'''
write_thermdat(thermdats=thermdats, filename=thermdats_out_path)
if show_plot:
	for thermdat in thermdats:
		thermdat.plot_thermo_model_and_empirical(Cp_units='J/mol/K', H_units='kJ/mol', S_units='J/mol/K', G_units='kJ/mol')
	plt.show()