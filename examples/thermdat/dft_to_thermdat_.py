from py_box3.thermo import dft_to_thermdat

parameters = {
	'input_path': '.\\thermdat_input_Fe3O4.csv', #csv file where the DFT information is located
	'ref_path': '.\\thermdat_ref.csv', #csv file that has the experimental gas-phase reference species
	'T_low': 500., #Lower temperature to fit in K
	'T_high': 1100., #Higher temperature to fit in K
	'write_files': True, #Whether a thermdat file should be written
	'out_path': 'thermdat_Fe3O4_DFT', #The name of the thermdat file
	'add_gas_species': True, #Whether there is a separate thermdat file with gas-phase species that was not included in input_path or ref_path (useful for species like He or N)
	'gas_path': 'thermdat_N2' #Name of the file that contains the gas species
}

dft_to_thermdat(**parameters)