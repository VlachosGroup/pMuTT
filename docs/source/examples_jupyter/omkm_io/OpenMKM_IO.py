
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
# 3. `beps` contains Bronsted-Evans-Polanyi relationships for reactions
# 4. `reactions` contains elementary steps
# 5. `lateral_interactions` contains lateral interactions between species
# 6. `phases` contains phases for the species

# First, we change the working directory to the location of the Jupyter notebook.

# In[1]:


import os
from pathlib import Path

# Find the location of Jupyter notebook
# Note that normally Python scripts have a __file__ variable but Jupyter notebook doesn't.
# Using pathlib can overcome this limiation
try:
    notebook_path = os.path.dirname(__file__)
except NameError:
    notebook_path = Path().resolve()
    
os.chdir(notebook_path)
input_path = './inputs/NH3_Input_Data.xlsx'


# Below is a helper function to print tables easily.

# In[2]:


import pandas as pd
from IPython.display import display

def disp_data(io, sheet_name):
    data = pd.read_excel(io=io, sheet_name=sheet_name, skiprows=[1])
    data = data.fillna(' ')
    display(data)    


# **References**

# In[3]:


disp_data(io=input_path, sheet_name='refs')


# **Species**

# In[4]:


disp_data(io=input_path, sheet_name='species')


# **BEPs**

# In[5]:


disp_data(io=input_path, sheet_name='beps')


# **Reactions**

# In[6]:


disp_data(io=input_path, sheet_name='reactions')


# **Lateral Interactions**

# In[7]:


disp_data(io=input_path, sheet_name='lateral_interactions')


# **Phases**

# In[8]:


disp_data(io=input_path, sheet_name='phases')


# ## Designate Units
# First, we will designate the units to write the CTI file.

# In[9]:


from pmutt.omkm.units import Units

units = Units(length='cm', quantity='mol', act_energy='kcal/mol', mass='g', energy='kcal/mol')


# ## Reading data
# Before we can initialize our species, we need the references.
# 
# ### Reading References (optional)
# We will open the [input spreadsheet](https://github.com/VlachosGroup/pmutt/blob/master/docs/source/examples_jupyter/openmkm_io/inputs/NH3_Input_Data.xlsx) and read the `refs` sheet.

# In[10]:


from pmutt.io.excel import read_excel
from pmutt.empirical.references import Reference, References

try:
    refs_data = read_excel(io=input_path, sheet_name='refs')
except:
    refs = None
else:
    refs = [Reference(**ref_data) for ref_data in refs_data]
    refs = References(references=refs)


# ### Reading Species

# In[11]:


from pmutt.empirical.nasa import Nasa

# Lower and higher temperatures
T_low = 298. # K
T_high = 800. # K

species_data = read_excel(io=input_path, sheet_name='species')
species = []
species_phases = {}
for ind_species_data in species_data:
    # Initialize NASA from statistical mechanical data
    ind_species = Nasa.from_model(T_low=T_low, T_high=T_high, references=refs,
                                  **ind_species_data)
    species.append(ind_species)

    # Group the species by phase for later use
    try:
        species_phases[ind_species.phase].append(ind_species)
    except KeyError:
        species_phases[ind_species.phase] = [ind_species]


# ### Adding species from other empirical sources (optional)

# In[12]:


import numpy as np
from pmutt.empirical.shomate import Shomate

Ar = Shomate(name='Ar', elements={'Ar': 1}, phase='gas', T_low=298., T_high=6000.,
             a=np.array([20.78600, 2.825911e-7, -1.464191e-7, 1.092131e-8, -3.661371e-8, -6.19735, 179.999, 0.]))

species.append(Ar)
species_phases['gas'].append(Ar)


# ### Reading BEP (optional)

# In[13]:


from pmutt.omkm.reaction import BEP

try:
    beps_data = read_excel(io=input_path, sheet_name='beps')
except:
    beps = None
    species_with_beps = species.copy()
else:
    beps = []
    for bep_data in beps_data:
        beps.append(BEP(**bep_data))

    # Combine species and BEPs to make reactions
    species_with_beps = species + beps


# ### Read reactions

# In[14]:


from pmutt import pmutt_list_to_dict
from pmutt.omkm.reaction import SurfaceReaction

# Convert species to dictionary for easier reaction assignment
species_with_beps_dict = pmutt_list_to_dict(species_with_beps)
try:
    reactions_data = read_excel(io=input_path, sheet_name='reactions')
except:
    reactions = None
else:
    reactions = []
    # Store information about phases for later retrieval
    reaction_phases = {}
    for reaction_data in reactions_data:
        reaction = SurfaceReaction.from_string(species=species_with_beps_dict,
                                               **reaction_data)
        reactions.append(reaction)
        # Assign phase information
        reaction_species = reaction.get_species(include_TS=True)
        for ind_species in reaction_species.values():
            try:
                phase = ind_species.phase
            except AttributeError:
                pass
            # Assign if key already exists
            if phase in reaction_phases:
                if reaction not in reaction_phases[phase]:
                    reaction_phases[phase].append(reaction)
            else:
                reaction_phases[phase] = [reaction]


# ### Read lateral interactions (optional)

# In[15]:


from pmutt.mixture.cov import PiecewiseCovEffect

try:
    interactions_data = read_excel(io=input_path, sheet_name='lateral_interactions')
except:
    interactions = None
else:
    interactions = []
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

# In[16]:


from pmutt.omkm.phase import IdealGas, InteractingInterface, StoichSolid

try:
    phases_data = read_excel(io=input_path, sheet_name='phases')
except:
    phases = None
else:
    phases = []
    # Group data related to previously collected data
    additional_fields = {'species': species_phases,
                         'reactions': reaction_phases,
                         'interactions': interaction_phases}
    for phase_data in phases_data:
        # Pre-processing relevant data
        phase_name = phase_data['name']
        phase_type = phase_data.pop('phase_type')

        # Add additional fields to phase data if present
        for field, phase_dict in additional_fields.items():
            try:
                phase_data[field] = phase_dict[phase_name]
            except (NameError, KeyError):
                pass

        # Create the appropriate object
        if phase_type == 'IdealGas':
            # Special rule where reactions are only in the gas phase if
            # all species belong to the gas phase
            del_indices = []
            for i, reaction in enumerate(phase_data['reactions']):
                # Reaction will be deleted if any of the species are a different phase
                valid_rxn = True
                for ind_species in reaction.get_species(include_TS=False).values():
                    try:
                        ind_species_phase = ind_species.phase
                    except AttributeError:
                        valid_rxn = False
                    else:
                        if ind_species_phase != phase_name:
                            valid_rxn = False
                    # Record reaction index if not valid
                    if not valid_rxn:
                        del_indices.append(i)
                        break
            # Delete reactions that do not qualify
            if len(del_indices) == len(phase_data['reactions']):
                phase_data.pop('reactions')
            else:
                for del_i in sorted(del_indices, reverse=True):
                    del phase_data['reactions'][del_i]
            phase = IdealGas(**phase_data)
        elif phase_type == 'StoichSolid':
            phase = StoichSolid(**phase_data)
        elif phase_type == 'InteractingInterface':
            phase = InteractingInterface(**phase_data)
        phases.append(phase)


# ## Write CTI File

# In[17]:


from pmutt.io.omkm import write_cti

cti_path = './outputs/thermo.cti'
use_motz_wise = True

write_cti(reactions=reactions, species=species, phases=phases, units=units,
          lateral_interactions=interactions, filename=cti_path,
          use_motz_wise=use_motz_wise)


# If you would prefer to return the file as a string instead of writing it, omit the ``filename``.

# In[18]:


print(write_cti(reactions=reactions, species=species, phases=phases, units=units,
                lateral_interactions=interactions, use_motz_wise=use_motz_wise))


# ## Write YAML File
# 
# The YAML file specifying the reactor configuration can also be written using the ``write_yaml`` function. Note that if:
# - ``units`` is not specified, float values are assumed to be in SI units
# - ``units`` is specified, float values are consistent with ``unit``'s attributes
# - you would like a quantity to have particular units, pass the value as a string with the units  (e.g. 10 cm3/s).

# In[19]:


from pmutt.io.omkm import write_yaml

yaml_path = './outputs/cstr.yaml'

write_yaml(filename=yaml_path, reactor_type='cstr', mode='isothermal',
           V=1., T=900., P=1., cat_abyv=1500, end_time=50, flow_rate=1.,
           transient=True, stepping='logarithmic', init_step=1e-15, atol=1e-15,
           rtol=1e-10, output_format='csv', phases=phases, units=units)


# If you would prefer to return the file as a string instead of writing it, omit the ``filename``.

# In[20]:


print(write_yaml(reactor_type='cstr', mode='isothermal', V=1., T=900., P=1., cat_abyv=1500,
                 end_time=50, flow_rate=1., transient=True, stepping='logarithmic',
                 init_step=1e-15, atol=1e-15, rtol=1e-10, output_format='csv', phases=phases,
                 units=units))

