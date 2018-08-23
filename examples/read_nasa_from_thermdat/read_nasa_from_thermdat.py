import os
from pprint import pprint
from matplotlib import pyplot as plt
from PyMuTT.io_.thermdat import read_thermdat
from PyMuTT.models.empirical.nasa import Nasa

base_path = os.path.dirname(__file__)
#Thermdat file from http://combustion.berkeley.edu/gri_mech/version30/files30/thermo30.dat
species = read_thermdat('{}/thermdat'.format(base_path))

#Printing information related to each specie
pprint(species)

#Plot an example of an imported NASA polynomial
species[1].plot_empirical()
plt.show()