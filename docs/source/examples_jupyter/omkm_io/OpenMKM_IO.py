#!/usr/bin/env python
# coding: utf-8

# # OpenMKM Input and Output
# This notebook describes pmutt's functionality to write OpenMKM CTI and YAML files. We will use the NH3 formation mechanism as a case study.
# 
# ## Topics Covered
# - Read species *ab-initio* data, reactions, lateral interactions, phases, reactor operating conditions, and desired units from a spreadsheet
# - Write the CTI file that can be read by OpenMKM
# - Write a YAML file that can be read by OpenMKM

# ## Input Spreadsheet
# All the data will be imported from the [`./inputs/NH3_Input_data.xlsx`](https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/omkm_io/inputs/NH3_Input_Data.xlsx) file. There are several sheets:
# 
# 1. `units` contains the units that types of quantities should be written
# 2. `refs` contains *ab-initio* and experimental data for a handful of gas species to calculate references (optional)
# 3. `species` contains *ab-initio* data for each specie
# 4. `beps` contains Bronsted-Evans-Polanyi relationships for reactions (optional)
# 5. `reactions` contains elementary steps 
# 6. `lateral_interactions` contains lateral interactions between species (optional)
# 7. `phases` contains phases for the species
# 8. `reactor` contains reactor operating conditions and solver tolerances
# 
# The ``refs``, ``beps`` and ``lateral_interactions`` sheets can be deleted and the code written below should still work.

# First, we change the working directory to the location of the Jupyter notebook.

import os
from pathlib import Path

import numpy as np
import pandas as pd
from IPython.display import display

from pmutt import pmutt_list_to_dict
from pmutt.empirical.nasa import Nasa
from pmutt.empirical.references import Reference, References
from pmutt.empirical.shomate import Shomate
from pmutt.io.excel import read_excel
from pmutt.io.omkm import (organize_phases, write_cti, write_thermo_yaml,
                           write_yaml)
from pmutt.mixture.cov import PiecewiseCovEffect
from pmutt.omkm.reaction import BEP, SurfaceReaction
from pmutt.omkm.units import Units

# Find the location of Jupyter notebook
# Note that normally Python scripts have a __file__ variable but Jupyter notebook doesn't.
# Using pathlib can overcome this limiation
try:
    notebook_path = os.path.dirname(__file__)
except NameError:
    notebook_path = Path().resolve()
    
os.chdir(notebook_path)
input_path = './inputs/NH3_Input_Data.xlsx'


units_data = read_excel(io=input_path, sheet_name='units')[0]
units = Units(**units_data)


# ### Reading References (optional)
# Second, we will open the input spreadsheet and read the `refs` sheet.
try:
    refs_data = read_excel(io=input_path, sheet_name='refs')
except:
    # If references are not used, skip this section
    print('The "refs" sheet could not be found in {}. Skiping references'.format(input_path))
    refs = None
else:
    refs = [Reference(**ref_data) for ref_data in refs_data]
    refs = References(references=refs)


# ### Reading Species
# 
# Third, we will use the ``refs`` defined before and the ``species`` sheet to convert statistical mechanical data to [``NASA``](https://vlachosgroup.github.io/pMuTT/api/empirical/nasa/pmutt.empirical.nasa.Nasa.html#pmutt.empirical.nasa.Nasa) objects.
# Read the species' data
species_data = read_excel(io=input_path, sheet_name='species')

# Create NASA polynomials from the species
species = [Nasa.from_model(references=refs, **ind_species_data)            for ind_species_data in species_data]


# ### Adding species from other empirical sources (optional)
# 
# Note that OpenMKM also supports [``Shomate``](https://vlachosgroup.github.io/pMuTT/api/empirical/shomate/pmutt.empirical.shomate.Shomate.html#pmutt.empirical.shomate.Shomate) and [``NASA9``](https://vlachosgroup.github.io/pMuTT/api/empirical/nasa/pmutt.empirical.nasa.Nasa9.html) objects. Below, we define a single ``Shomate`` species.
Ar = Shomate(name='Ar', elements={'Ar': 1}, phase='gas', T_low=298., T_high=6000.,
             a=np.array([20.78600, 2.825911e-7, -1.464191e-7, 1.092131e-8, -3.661371e-8, -6.19735, 179.999, 0.]))

species.append(Ar)


# ### Reading BEP (optional)
# 
# Next, we read the BEP relationships to include.
try:
    beps_data = read_excel(io=input_path, sheet_name='beps')
except:
    print('The "beps" sheet could not be found in {}. Skiping BEPs'.format(input_path))
    beps = None
    species_with_beps = species.copy()
else:
    beps = [BEP(**bep_data) for bep_data in beps_data]
    species_with_beps = species + beps


# ### Read reactions
# 
# Then, we read the reactions to include.
# Convert species to dictionary for easier reaction assignment
species_with_beps_dict = pmutt_list_to_dict(species_with_beps)

reactions_data = read_excel(io=input_path, sheet_name='reactions')
reactions = [SurfaceReaction.from_string(species=species_with_beps_dict, **reaction_data)              for reaction_data in reactions_data]


# ### Read lateral interactions (optional)
# 
# After, we read lateral interactions to include.
try:
    interactions_data = read_excel(io=input_path, sheet_name='lateral_interactions')
except:
    # If no lateral interactions exist, skip this section
    print('The "lateral_interactions" sheet could not be found in {}. Skiping lateral interactions'.format(input_path))
    interactions = None
else:
    interactions = [PiecewiseCovEffect(**interaction_data) for interaction_data in interactions_data]


# ### Reading Phases
# 
# Finally, we read the phases data from Excel and organize it for use in OpenMKM.
# Read data from Excel sheet about phases
phases_data = read_excel(io=input_path, sheet_name='phases')
phases = organize_phases(phases_data, species=species, reactions=reactions, interactions=interactions)


# ## Write Reactor YAML File
# 
# The YAML file specifying the reactor configuration can be written using the [``write_yaml``](https://vlachosgroup.github.io/pMuTT/api/kinetic_models/omkm/pmutt.io.omkm.write_yaml.html) function. Note that if:
# - ``units`` is not specified, float values are assumed to be in SI units
# - ``units`` is specified, float values are consistent with ``unit``'s attributes
# - you would like a quantity to have particular units, pass the value as a string with the units  (e.g. "10 cm3/s").
Path('./outputs').mkdir(exist_ok=True)
yaml_path = './outputs/reactor.yaml'
reactor_data = read_excel(io=input_path, sheet_name='reactor')[0]
write_yaml(filename=yaml_path, phases=phases, units=units, **reactor_data)


# If you would prefer to return the file as a string instead of writing it, omit the ``filename``.
print(write_yaml(phases=phases, units=units, **reactor_data))

# ## Write Thermo/Kinetic YAML File
# 
# As of OpenMKM version 0.6.0 onwards, the thermodynamic and kinetic parameters can be written as a YAML file. We recommend using this format over the older CTI format. To generate the Thermo/Kinetic YAML file using pMuTT, use the [``write_thermo_yaml``](https://vlachosgroup.github.io/pMuTT/api/kinetic_models/omkm/pmutt.io.omkm.write_thermo_yaml.html) function
T = reactor_data['T']
write_thermo_yaml(T=T,
                  phases=phases,
                  species=species,
                  reactions=reactions,
                  lateral_interactions=interactions,
                  units=units,
                  filename='./outputs/thermo.yaml')


# Like before, omitting the ``filename`` parameter returns a string
print(write_thermo_yaml(phases=phases,
                        species=species,
                        reactions=reactions,
                        lateral_interactions=interactions,
                        units=units))


# ## Write CTI File
# 
# The CTI file species the thermodynamics and kinetics of the system. It can be written using [``write_cti``](https://vlachosgroup.github.io/pMuTT/api/kinetic_models/omkm/pmutt.io.omkm.write_cti.html#pmutt.io.omkm.write_cti). Note that we take the reactor operating conditions previously read for the YAML file to calculate thermodynamic and kinetic parameters.
cti_path = './outputs/thermo.cti'
use_motz_wise = True

write_cti(reactions=reactions, species=species, phases=phases, units=units,
         lateral_interactions=interactions, filename=cti_path,
         use_motz_wise=use_motz_wise, T=T, P=1.)


# Like before, omitting the ``filename`` parameter returns a string.
print(write_cti(reactions=reactions, species=species, phases=phases, units=units,
               lateral_interactions=interactions, use_motz_wise=use_motz_wise))
