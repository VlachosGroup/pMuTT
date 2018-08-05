from PyMuTT.io_.thermdat import read_thermdat
from PyMuTT.models.empirical.nasa import Nasa

#Thermdat file from http://combustion.berkeley.edu/gri_mech/version30/files30/thermo30.dat
species = read_thermdat('thermdat')

#Printing information related to each specie
for specie in species:
	print('Name: {}'.format(specie.name))
	for key, val in specie.__dict__.items():
		if key != 'name':
			print('\t{}\t{}'.format(key, val))
	