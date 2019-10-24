
# coding: utf-8

# # OpenMKM Input and Output
# This notebook describes pmutt's functionality to read and write OpenMKM CTI files. We will use the NH3 formation mechanism as a case study.
# 
# ## Topics Covered
# - Read species *ab-initio* data, reactions, lateral interactions and phases from a spreadsheet
# - Write the CTI file that can be read by OpenMKM

# ## Input Spreadsheet
# All the data will be imported from the [`./inputs/NH3_Input_data.xlsx`](https://github.com/VlachosGroup/pmutt/blob/master/docs/source/examples_jupyter/openmkm_io/inputs/NH3_Input_Data.xlsx) file. There are five sheets:
# 1. `refs` contains *ab-initio* and experimental data for a handful of gas species to calculate references
# 2. `species` contains *ab-initio* data for each specie
# 3. `reactions` contains elementary steps
# 4. `phases` contains phases for the species
# 5. `lateral_interactions` contains lateral interactions between species
# 
# The contents are displayed below:

# **References**
# 
# | name | elements.N | elements.H | elements.Ru | T_ref  | HoRT_ref     | potentialenergy | symmetrynumber | statmech_model | atoms         | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber |
# |------|------------|------------|-------------|--------|--------------|-----------------|----------------|----------------|---------------|----------------|----------------|----------------|----------------|
# | N2   | 2          | 0          | 0           | 298.15 | 0            | -16.63          | 2              | IdealGas       | ./N2/CONTCAR  | 2744           |
# | NH3  | 1          | 3          | 0           | 298.15 | -18.38025311 | -19.54          | 3              | IdealGas       | ./NH3/CONTCAR | 3534           | 3464           | 1765           | 1139           |
# | H2   | 0          | 2          | 0           | 298.15 | 0            | -6.7700         | 2              | IdealGas       | ./H2/CONTCAR  | 4342           |
# | Ru   | 0          | 0          | 1           | 298.15 | 0.0000       |                 |                | Placeholder    |

# **Species**
# 
# | name       | elements.N | elements.H | elements.Ru | phase   | statmech_model | symmetrynumber | atoms         | potentialenergy | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber |
# |------------|------------|------------|-------------|---------|----------------|----------------|---------------|-----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|
# | N2         | 2          |            |             | gas     | IdealGas       | 2              | ./N2/CONTCAR  | -16.63          | 2744.00        |                |                |                |                |                |                |                |                |                |                |                |
# | NH3        | 1          | 3          |             | gas     | IdealGas       | 3              | ./NH3/CONTCAR | -19.54          | 3534.00        | 3464.00        | 1765.00        | 1139.00        |                |                |                |                |                |                |                |                |
# | H2         |            | 2          |             | gas     | IdealGas       | 2              | ./H2/CONTCAR  | -6.77           | 4342.00        |                |                |                |                |                |                |                |                |                |                |                |
# | N2(S)      | 2          |            |             | terrace | Harmonic       |                |               | -17.24          | 2197.19        | 360.42         | 347.34         | 335.67         | 62.08          | 32.18          |                |                |                |                |                |                |
# | N(S)       | 1          |            |             | terrace | Harmonic       |                |               | -9.34           | 549.11         | 538.56         | 504.32         | 475.81         | 459.08         | 410.02         |                |                |                |                |                |                |
# | H(S)       |            | 1          |             | terrace | Harmonic       |                |               | -4.00           | 1003.51        | 625.55         | 616.29         |                |                |                |                |                |                |                |                |                |
# | NH3(S)     | 1          | 3          |             | terrace | Harmonic       |                |               | -20.43          | 3491.09        | 3488.82        | 3364.52        | 1583.52        | 1582.07        | 1124.22        | 570.21         | 567.22         | 333.09         | 122.86         | 83.83          | 70.63          |
# | NH2(S)     | 1          | 2          |             | terrace | Harmonic       |                |               | -16.59          | 3469.30        | 3381.05        | 1503.02        | 698.87         | 625.60         | 615.94         | 475.13         | 298.12         | 153.25         |                |                |                |
# | NH(S)      | 1          | 1          |             | terrace | Harmonic       |                |               | -13.21          | 3403.13        | 718.18         | 710.58         | 528.53         | 415.20         | 410.13         |                |                |                |                |                |                |
# | TS1_NH3(S) | 1          | 3          |             |         | Harmonic       |                |               | -19.24          | 3453.41        | 3355.67        | 1723.85        | 1487.95        | 959.15         | 888.95         | 594.09         | 428.43         | 227.03         | 206.05         | 142.14         |                |
# | TS2_NH2(S) | 1          | 2          |             |         | Harmonic       |                |               | -15.87          | 3426.44        | 1293.72        | 922.83         | 660.97         | 525.60         | 496.84         | 330.67         | 290.28         |                |                |                |                |
# | TS3_NH(S)  | 1          | 1          |             |         | Harmonic       |                |               | -11.93          | 1201.60        | 491.57         | 462.02         | 402.16         | 242.14         |                |                |                |                |                |                |                |
# | TS4_N2(S)  | 2          |            |             |         | Harmonic       |                |               | -14.67          | 485.61         | 392.98         | 386.19         | 280.94         | 168.43         |                |                |                |                |                |                |                |
# | RU(S)      |            |            | 1           | terrace | Placeholder    |                |               |                 |                |                |                |                |                |                |                |                |                |                |                |                |
# | RU(B)      |            |            | 1           | bulk    | Placeholder    |                |               |                 |                |                |                |                |                |                |                |                |                |                |                |                |

# **Phases**
# 
# | name    | phase_type           | density | site_density | reactions | interactions | list.phases | list.phases | note     |
# |---------|----------------------|---------|--------------|-----------|--------------|-------------|-------------|----------|
# | gas     | IdealGas             |         |              |           |              |             |             |          |
# | bulk    | StoichSolid          | 12.4    |              |           |              |             |             | Ru Metal |
# | terrace | InteractingInterface |         | 2.17E-09     | all       | all          | gas         | bulk        | Ru(0001) |

# **Reactions**
# 
# | reaction_str                                        | is_adsorption |
# |-----------------------------------------------------|---------------|
# | H2 + 2RU(S) = 2H(S) + 2RU(B)                        | TRUE          |
# | N2 + RU(S)  = N2(S) + RU(B)                         | TRUE          |
# | NH3 + RU(S) = NH3(S) + RU(B)                        | TRUE          |
# | NH3(S) + RU(S)= TS1_NH3(S) = NH2(S) + H(S) + RU(B)  | FALSE         |
# | NH2(S) + RU(S) = TS2_NH2(S) = NH(S)  + H(S) + RU(B) | FALSE         |
# | NH(S)  + RU(S) = TS3_NH(S) = N(S)   + H(S) + RU(B)  | FALSE         |
# | 2N(S) + RU(B) = TS4_N2(S) = N2(S)  +  RU(S)         | FALSE         |

# **Lateral Interactions**
# 
# | name_i | name_j | list.intervals | list.slopes |
# |--------|--------|----------------|-------------|
# | N(S)   | N(S)   | 0              | -52.6       |
# | N(S)   | H(S)   | 0              | -17.7       |
# | H(S)   | N(S)   | 0              | -17.7       |
# | H(S)   | H(S)   | 0              | -3          |
# | NH2(S) | N(S)   | 0              | -20.7       |

# ## Designate Units
# First, we will designate the units to write the CTI file.

# In[1]:


from pmutt.omkm.units import Units

units = Units(length='cm', quantity='mol', act_energy='kcal/mol', mass='g', energy='kcal/mol')


# ## Reading data
# Before we can initialize our species, we need the references.
# 
# ### Reading References
# We will open the [input spreadsheet](https://github.com/VlachosGroup/pmutt/blob/master/docs/source/examples_jupyter/openmkm_io/inputs/NH3_Input_Data.xlsx) and read the `refs` sheet.

# In[2]:


import os
from pathlib import Path

from pmutt.io.excel import read_excel
from pmutt.empirical.references import Reference, References

# Find the location of Jupyter notebook
# Note that normally Python scripts have a __file__ variable but Jupyter notebook doesn't.
# Using pathlib can overcome this limiation
#notebook_path = Path().resolve()
notebook_path = os.path.dirname(__file__)
os.chdir(notebook_path)
input_path = './inputs/NH3_Input_Data.xlsx'

refs_data = read_excel(io=input_path, sheet_name='refs')
refs = [Reference(**ref_data) for ref_data in refs_data]
refs = References(references=refs)


# ### Reading Species

# In[3]:


from pmutt.empirical.nasa import Nasa

# Lower and higher temperatures
T_low = 298. # K
T_high = 800. # K

species_data = read_excel(io=input_path, sheet_name='species')
species = []
species_phases = {}
for ind_species_data in species_data:
    # Initialize NASA from statistical mechanical data
    ind_species = Nasa.from_model(T_low=T_low, T_high=T_high, references=refs, **ind_species_data)
    species.append(ind_species)

    # Group the species by phase for later use
    try:
        species_phases[ind_species.phase].append(ind_species)
    except KeyError:
        species_phases[ind_species.phase] = [ind_species]


# ### Adding species from other empirical sources

# In[4]:


import numpy as np
from pmutt.empirical.shomate import Shomate

Ar = Shomate(name='Ar', elements={'Ar': 1}, phase='gas', T_low=298., T_high=6000.,
             a=np.array([20.78600, 2.825911e-7, -1.464191e-7, 1.092131e-8, -3.661371e-8, -6.19735, 179.999, 0.]))

species.append(Ar)
species_phases['gas'].append(Ar)


# ### Reading BEP

# In[5]:


from pmutt.omkm.reaction import BEP

beps_data = read_excel(io=input_path, sheet_name='beps')
beps = []
for bep_data in beps_data:
    beps.append(BEP(**bep_data))

# Combine species and BEPs to make reactions
species_with_beps = species + beps


# ### Read reactions

# In[6]:


from pmutt import pmutt_list_to_dict
from pmutt.omkm.reaction import SurfaceReaction

# Convert species to dictionary for easier reaction assignment
species_with_beps_dict = pmutt_list_to_dict(species_with_beps)
reactions_data = read_excel(io=input_path, sheet_name='reactions')
reactions = []
# Store information about phases for later retrieval
reaction_phases = {}
for reaction_data in reactions_data:
    reaction = SurfaceReaction.from_string(species=species_with_beps_dict,
                                           **reaction_data)
    reactions.append(reaction)
    # Assign phase information
    reaction_species = reaction.get_species(include_TS=True)
    for ind_species in reaction_species:
        try:
            phase = species_with_beps_dict[ind_species].phase
        except AttributeError:
            pass
        # Assign if key already exists
        if phase in reaction_phases:
            if reaction not in reaction_phases[phase]:
                reaction_phases[phase].append(reaction)
        else:
            reaction_phases[phase] = [reaction]


# ### Read lateral interactions

# In[7]:


from pmutt.mixture.cov import PiecewiseCovEffect

interactions = []
interactions_data = read_excel(io=input_path, sheet_name='lateral_interactions')
interaction_phases = {}
for interaction_data in interactions_data:
    interaction = PiecewiseCovEffect(**interaction_data)
    interactions.append(interaction)

    # Assign phase information
    phase = species_with_beps_dict[interaction.name_i].phase
    # Assign if key already exists
    if phase in interaction_phases:
        if interaction not in interaction_phases[phase]:
            interaction_phases[phase].append(interaction)
    else:
        interaction_phases[phase] = [interaction]


# ### Reading Phases

# In[8]:


from pmutt.omkm.phase import IdealGas, InteractingInterface, StoichSolid

phases_data = read_excel(io=input_path, sheet_name='phases')
phases = []
for phase_data in phases_data:
    # Pre-processing relevant data
    phase_name = phase_data['name']
    phase_type = phase_data.pop('phase_type')
    phase_data['species'] = species_phases[phase_name]

    # Create the appropriate object
    if phase_type == 'IdealGas':
        phase = IdealGas(**phase_data)
    elif phase_type == 'StoichSolid':
        phase = StoichSolid(**phase_data)
    elif phase_type == 'InteractingInterface':
        phase_data['reactions'] = reaction_phases[phase_name]
        phase_data['interactions'] = interaction_phases[phase_name]
        phase = InteractingInterface(**phase_data)
    phases.append(phase)


# ## Write CTI File

# In[9]:


from pmutt.io.omkm import write_cti

output_path = './outputs/input.cti'
use_motz_wise = True

write_cti(reactions=reactions, species=species, phases=phases, units=units,
          lateral_interactions=interactions, filename=output_path,
          use_motz_wise=use_motz_wise)


# ## Output CTI File
# If you would prefer to return the file as a string instead of writing it, omit the ``filename``.

# In[10]:


print(write_cti(reactions=reactions, species=species, phases=phases, units=units,
                lateral_interactions=interactions, use_motz_wise=use_motz_wise))

