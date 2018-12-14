import os
from pprint import pprint
from pMuTT.io_.excel import read_excel
from pMuTT.empirical.nasa import Nasa

base_path = os.path.dirname(__file__)
#Thermdat file from http://combustion.berkeley.edu/gri_mech/version30/files30/thermo30.dat
species_data = read_excel('{}/input_data.xlsx'.format(base_path))
species = [Nasa(**specie_data) for specie_data in species_data]

#Printing information related to each specie
pprint(species)