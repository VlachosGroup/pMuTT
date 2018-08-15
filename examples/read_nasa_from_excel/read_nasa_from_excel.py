import os
from pprint import pprint
from PyMuTT.io_.excel import read_excel
from PyMuTT.models.empirical.nasa import Nasa

base_path = os.path.dirname(__file__)
#Thermdat file from http://combustion.berkeley.edu/gri_mech/version30/files30/thermo30.dat
species_data =	read_excel('{}/input_data.xlsx'.format(base_path))
species = [Nasa(**specie_data) for specie_data in species_data]

#Printing information related to each specie
pprint(species_data)
# for specie in species:
# 	print('Name: {}'.format(specie.name))
# 	for key, val in specie.__dict__.items():
# 		if key != 'name':
# 			print('\t{}\t{}'.format(key, val))
	