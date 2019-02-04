import matplotlib.pyplot as plt
import os
from pprint import pprint

from pMuTT.io_.excel import read_excel
from pMuTT.io_.thermdat import write_thermdat
from pMuTT.empirical.nasa import Nasa
from pMuTT.empirical.references import Reference, References

'''
User inputs
'''
base_path = os.path.dirname(__file__)

# Reference information
refs_path = '{}/references.xlsx'.format(base_path)

# Input information
species_path = '{}/thermdat_input.xlsx'.format(base_path)
T_low = 200.
T_high = 1100.  # K

# Output information
thermdat_path = '{}/thermdat'.format(base_path)

# Miscellaneous options
show_plot = True
write_date = True

'''
Processing References
'''
# Import from excel
refs_input = read_excel(io=refs_path)
pprint(refs_input)
refs = References(
        references=[Reference(**ref_input) for ref_input in refs_input],
        descriptor='elements')
print('Reference Input:')
pprint(refs_input)
print('Reference Data:')
pprint(refs.references[0])

'''
Processing Input Species
'''
# Import from excel
species_data = read_excel(io=species_path)
# pprint([specie_data['elements'] for specie_data in species_data])
species = []
for specie_data in species_data:
    species.append(Nasa.from_statmech(references=refs, T_low=T_low,
                                      T_high=T_high, **specie_data))
species = [Nasa.from_statmech(references=refs, T_low=T_low, T_high=T_high,
                              **specie_data) for specie_data in species_data]
print('Species Input:')
pprint(species)

'''
Printing Out Results
'''
write_thermdat(nasa_species=species, filename=thermdat_path,
               write_date=write_date)
if show_plot:
    for specie in species:
        specie.plot_statmech_and_empirical(Cp_units='J/mol/K',
                                           H_units='kJ/mol',
                                           S_units='J/mol/K',
                                           G_units='kJ/mol')
    plt.show()
    