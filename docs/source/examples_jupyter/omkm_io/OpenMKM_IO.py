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
# 2. `refs` contains *ab-initio* and experimental data for a handful of gas species to calculate references
# 3. `species` contains *ab-initio* data for each specie
# 4. `beps` contains Bronsted-Evans-Polanyi relationships for reactions
# 5. `reactions` contains elementary steps
# 6. `lateral_interactions` contains lateral interactions between species
# 7. `phases` contains phases for the species
# 8. `reactor` contains reactor operating conditions and solver tolerances
# 

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


# And we define a helper function to print data from the Excel spreadsheet easily.

# In[2]:


import pandas as pd
from IPython.display import display

def disp_data(io, sheet_name):
    try:
        data = pd.read_excel(io=io, sheet_name=sheet_name, skiprows=[1])
    except:
        print('Sheet could not be found.')
    else:
        data = data.fillna(' ')
        display(data)


# Below, we show the contents of the Excel sheets.

# **Units**

# In[3]:


disp_data(io=input_path, sheet_name='units')


# **References**

# In[4]:


disp_data(io=input_path, sheet_name='refs')


# **Species**
# 
# In this mechanism, we have species existing on terraces and steps. We also define the transition states of species. 
# 
# Note that later we will use BEPs to calculate barriers. Hence, these transition state species will actually be ignored. You should use either transition state species or BEP relationships to calculate barriers.

# In[5]:


disp_data(io=input_path, sheet_name='species')


# **BEPs**

# In[6]:


disp_data(io=input_path, sheet_name='beps')


# **Reactions**
# 
# Note that reactions with two '=' signs indicate it has a transition state.

# In[7]:


disp_data(io=input_path, sheet_name='reactions')


# **Lateral Interactions**
# 
# Currently we use piece-wise linear equations for lateral interactions. Here, we only define one line but additional ``list.intervals`` and ``list.slopes`` columns can be added for more complicated behavior.

# In[8]:


disp_data(io=input_path, sheet_name='lateral_interactions')


# **Phases**

# In[9]:


disp_data(io=input_path, sheet_name='phases')


# **Reactor**

# In[10]:


disp_data(io=input_path, sheet_name='reactor')


# ## Reading data
# 
# Throughout this exercise, we will use ``pmutt.io.read_excel`` to extract the data from the Excel spreadsheet.

# ### Designate Units
# First, we will designate the units to write the CTI and YAML file.

# In[11]:


from pmutt.io.excel import read_excel
from pmutt.omkm.units import Units

units_data = read_excel(io=input_path, sheet_name='units')[0]
units = Units(**units_data)


# ### Reading References (optional)
# We will open the [input spreadsheet](https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/omkm_io/inputs/NH3_Input_Data.xlsx) and read the `refs` sheet.

# In[12]:


from pmutt.empirical.references import Reference, References

try:
    refs_data = read_excel(io=input_path, sheet_name='refs')
except:
    # If references are not used, skip this section
    refs = None
else:
    refs = [Reference(**ref_data) for ref_data in refs_data]
    refs = References(references=refs)


# ### Reading Species

# In[13]:


from pmutt.empirical.nasa import Nasa

# Read the species' data
species_data = read_excel(io=input_path, sheet_name='species')

# Create NASA polynomials from the species
species = [Nasa.from_model(references=refs, **ind_species_data)            for ind_species_data in species_data]


# ### Adding species from other empirical sources (optional)

# In[14]:


import numpy as np
from pmutt.empirical.shomate import Shomate

Ar = Shomate(name='Ar', elements={'Ar': 1}, phase='gas', T_low=298., T_high=6000.,
             a=np.array([20.78600, 2.825911e-7, -1.464191e-7, 1.092131e-8, -3.661371e-8, -6.19735, 179.999, 0.]))

species.append(Ar)


# ### Reading BEP (optional)

# In[15]:


from pmutt.omkm.reaction import BEP

try:
    beps_data = read_excel(io=input_path, sheet_name='beps')
except:
    beps = None
    species_with_beps = species.copy()
else:
    beps = [BEP(**bep_data) for bep_data in beps_data]
    species_with_beps = species + beps


# ### Read reactions

# In[16]:


from pmutt import pmutt_list_to_dict
from pmutt.omkm.reaction import SurfaceReaction

# Convert species to dictionary for easier reaction assignment
species_with_beps_dict = pmutt_list_to_dict(species_with_beps)

reactions_data = read_excel(io=input_path, sheet_name='reactions')
reactions = [SurfaceReaction.from_string(species=species_with_beps_dict, **reaction_data)              for reaction_data in reactions_data]


# ### Read lateral interactions (optional)

# In[17]:


from pmutt.mixture.cov import PiecewiseCovEffect

try:
    interactions_data = read_excel(io=input_path, sheet_name='lateral_interactions')
except:
    # If no lateral interactions exist, skip this section
    interactions = None
else:
    interactions = [PiecewiseCovEffect(**interaction_data) for interaction_data in interactions_data]


# ### Reading Phases

# In[18]:


from pmutt.io.omkm import organize_phases

# Read data from Excel sheet about phases
phases_data = read_excel(io=input_path, sheet_name='phases')
phases = organize_phases(phases_data, species=species, reactions=reactions, interactions=interactions)


# If you would prefer to return the file as a string instead of writing it, omit the ``filename``.

# ## Write YAML File
# 
# The YAML file specifying the reactor configuration can also be written using the ``write_yaml`` function. Note that if:
# - ``units`` is not specified, float values are assumed to be in SI units
# - ``units`` is specified, float values are consistent with ``unit``'s attributes
# - you would like a quantity to have particular units, pass the value as a string with the units  (e.g. "10 cm3/s").

# In[19]:


from pmutt.io.omkm import write_yaml

yaml_path = './outputs/reactor.yaml'
reactor_data = read_excel(io=input_path, sheet_name='reactor')[0]
write_yaml(filename=yaml_path, phases=phases, units=units, **reactor_data)


# If you would prefer to return the file as a string instead of writing it, omit the ``filename``.

# In[20]:


print(write_yaml(phases=phases, units=units, **reactor_data))


# ## Write CTI File

# In[21]:


from pmutt.io.omkm import write_cti

cti_path = './outputs/thermo.cti'
use_motz_wise = True
T = reactor_data['T']
P = reactor_data['P']

write_cti(reactions=reactions, species=species, phases=phases, units=units,
          lateral_interactions=interactions, filename=cti_path,
          use_motz_wise=use_motz_wise, T=T, P=P)


# In[22]:


print(write_cti(reactions=reactions, species=species, phases=phases, units=units,
                lateral_interactions=interactions, use_motz_wise=use_motz_wise))

