import os
from pprint import pprint

from pmutt import pmutt_list_to_dict
from pmutt.empirical.nasa import Nasa
from pmutt.empirical.references import Reference, References
from pmutt.io.excel import read_excel
from pmutt.io.omkm import write_cti
from pmutt.mixture.cov import PiecewiseCovEffect
from pmutt.omkm.phase import IdealGas, InteractingInterface, StoichSolid
from pmutt.omkm.reaction import SurfaceReaction
from pmutt.omkm.units import Units

os.chdir(os.path.dirname(__file__))

'''Input Parameters'''
input_path = './inputs/NH3_Input_Data.xlsx'
output_path = './outputs/cti.out'
T_low = 298. # K
T_high = 800. # K
use_motz_wise = True


'''Designate Units'''
units = Units(length='cm', quantity='mol', act_energy='kcal/mol', mass='g',
              energy='kcal/mol')

'''Initialize references'''
refs_data = read_excel(io=input_path, sheet_name='refs')
refs = [Reference(**ref_data) for ref_data in refs_data]
refs = References(references=refs)

'''Initialize species'''
species_data = read_excel(io=input_path, sheet_name='species')
species = []
species_phases = {}
for ind_species_data in species_data:
    ind_species = Nasa.from_statmech(T_low=T_low, T_high=T_high,
                                     references=refs, **ind_species_data)
    species.append(ind_species)

    # Group the species by phase
    try:
        species_phases[ind_species.phase].append(ind_species)
    except KeyError:
        species_phases[ind_species.phase] = [ind_species]


'''Initialize phases'''
phases_data = read_excel(io=input_path, sheet_name='phases')
phases_dict = {}
for phase_data in phases_data:
    phase_name = phase_data['name']
    phase_type = phase_data.pop('phase_type')
    phase_data['species'] = species_phases[phase_name]
    if phase_type == 'IdealGas':
        phase = IdealGas(**phase_data)
    elif phase_type == 'StoichSolid':
        phase = StoichSolid(**phase_data)
    elif phase_type == 'InteractingInterface':
        phase = InteractingInterface(**phase_data)
    phases_dict[phase.name] = phase

'''Assign phases to objects requiring other phases (e.g. InteractingInterface)'''
for name, phase in phases_dict.items():
    # See if object accepts other phases as an attribute
    try:
        req_phase_names = phase.phases
    except AttributeError:
        continue

    phases_objs = []
    for req_phase_name in req_phase_names:
        phases_objs.append(phases_dict[req_phase_name])
    phase.phases = phases_objs
# Convert phases_dict to list for write_cti
phases = list(phases_dict.values())

'''Initialize reactions'''
# Convert species to dictionary for easier reaction assignment
species_dict = pmutt_list_to_dict(species)
reactions_data = read_excel(io=input_path, sheet_name='reactions')
reactions = []
for reaction_data in reactions_data:
    reaction = SurfaceReaction.from_string(species=species_dict,
                                           **reaction_data)
    reactions.append(reaction)

'''Initialize lateral interactions'''
interactions = []
interactions_data = read_excel(io=input_path, sheet_name='lateral_interactions')
for interaction_data in interactions_data:
    interaction = PiecewiseCovEffect(**interaction_data)
    interactions.append(interaction)

'''Write CTI file'''
write_cti(reactions=reactions, species=species, phases=phases, units=units,
          lateral_interactions=interactions, filename=output_path,
          use_motz_wise=use_motz_wise)