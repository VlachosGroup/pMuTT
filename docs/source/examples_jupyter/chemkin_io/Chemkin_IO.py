#!/usr/bin/env python
# coding: utf-8

# # Chemkin Input and Output
# This notebook describes pmutt's functionality to read and write Chemkin files. We will use the NH3 formation mechanism as a case study.
# 
# ## Topics Covered
# - Read species *ab-initio* data, reactions, and catalyst sites from a spreadsheet
# - Write the thermdat, gas.inp, surf.inp, T_flow.inp, EAg.inp, EAs.inp, tube_mole.inp files

# ## Input Spreadsheet
# All the data will be imported from the [`./inputs/NH3_Input_Data.xlsx`](https://github.com/VlachosGroup/pmutt/blob/master/docs/source/examples_jupyter/chemkin_io/inputs/NH3_Input_Data.xlsx) file. There are four sheets:
# 1. `cat_sites` contains catalyst site properties for the adsorbed species
# 2. `refs` contains *ab-initio* and experimental data for a handful of gas species to calculate references
# 3. `species` contains *ab-initio* data for each specie
# 4. `reactions` contains elementary steps

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
excel_path = './inputs/NH3_Input_Data.xlsx'


# Below is a helper function to print tables easily.

# In[2]:


import pandas as pd
from IPython.display import display

def disp_data(io, sheet_name):
    data = pd.read_excel(io=io, sheet_name=sheet_name, skiprows=[1])
    data = data.fillna(' ')
    display(data)


# **Catalytic Sites**

# In[3]:


disp_data(io=excel_path, sheet_name='cat_sites')


# **References**

# In[4]:


disp_data(io=excel_path, sheet_name='refs')


# **Species**

# In[5]:


disp_data(io=excel_path, sheet_name='species')


# **Reactions**

# In[6]:


disp_data(io=excel_path, sheet_name='reactions')


# ## Reading data
# Before we can initialize our species, we need the catalytic sites and the references.
# 
# ### Reading Catalytic Sites

# In[7]:


import os
from pprint import pprint
from pathlib import Path
from pmutt.io.excel import read_excel
from pmutt.chemkin import CatSite

cat_site_data = read_excel(io=excel_path, sheet_name='cat_sites')[0]
cat_site = CatSite(**cat_site_data)

# Print the properties of the catalyst site
pprint(cat_site.to_dict())


# ### Reading reference species

# In[8]:


from pmutt.empirical.references import Reference, References

references_data = read_excel(io=excel_path, sheet_name='refs')

# Convert data to Reference objects and put them in a list
refs_list = [Reference(**ref_data) for ref_data in references_data]

# Assign the Reference objects to a References object so offsets can be calculated
refs = References(references=refs_list)

# Print out the offsets calculated
print(refs.offset)


# ### Reading species

# In[9]:


from pmutt.empirical.nasa import Nasa

# Range of data to fit the Nasa polynomials
T_low = 298. # K
T_high = 800. # K

species_data = read_excel(io=excel_path, sheet_name='species')
species = []
for specie_data in species_data:
    specie = Nasa.from_model(T_low=T_low, T_high=T_high, references=refs, **specie_data)
    # If the species is a surface species, assign the catalyst site specified above
    if specie.phase.lower() == 's':
        specie.cat_site = cat_site
        specie.n_sites = 1
    species.append(specie)


# The warning above is typical when empirical objects are fitting to `StatMech` objects with the `placeholder` preset.

# ### Reading reactions

# In[10]:


from pmutt import pmutt_list_to_dict
from pmutt.reaction import ChemkinReaction, Reactions

# Convert list of Nasa polynomials to dictionary of Nasa polynomials
species_dict = pmutt_list_to_dict(species)

reactions_data = read_excel(io=excel_path, sheet_name='reactions')
reactions_list = []
for reaction_data in reactions_data:
    reaction = ChemkinReaction.from_string(species=species_dict, **reaction_data)
    reactions_list.append(reaction)
reactions = Reactions(reactions=reactions_list)


# ## Writing Chemkin files
# Now that we have all the required objects, we can write the output files. All outputs can be found in the [./outputs folder](https://github.com/VlachosGroup/pmutt/blob/master/docs/source/examples_jupyter/chemkin_io/outputs).

# ### Writing thermdat

# In[11]:


from pmutt.io.thermdat import write_thermdat

write_thermdat(filename='./outputs/thermdat', nasa_species=species)


# The thermdat file can be return as a string by omitting ``filename``.

# In[12]:


thermdat_str = write_thermdat(nasa_species=species)
print(thermdat_str)


# ### Writing gas.inp and surf.inp

# In[13]:


from pmutt.io import chemkin as ck_io

ck_io.write_gas(filename='./outputs/gas.inp',
                nasa_species=species,
                reactions=reactions,
                act_method_name='get_G_act',
                act_unit='kcal/mol')
ck_io.write_surf(filename='./outputs/surf.inp',
                 reactions=reactions,
                 act_method_name='get_G_act',
                 ads_act_method='get_H_act',
                 act_unit='kcal/mol')


# <a id='act_method_name_explanation'></a>
# Note that `act_method_name` is 'get_G_act'. We use this formalism here since we do not include entropic effects in the preexponential factor.
# 
# Similarly to ``write_thermdat``, the gas.inp and surf.inp file can written as a string by omitting the filename. Note there are no gas-phase reactions.

# In[14]:


gas_file = ck_io.write_gas(nasa_species=species,
                           reactions=reactions,
                           act_method_name='get_G_act',
                           act_units='kcal/mol')
print(gas_file)


# In[15]:


surf_file = ck_io.write_surf(reactions=reactions,
                             act_method_name='get_G_act',
                             ads_act_method='get_H_act',
                             act_unit='kcal/mol')
print(surf_file)


# ### Writing T_flow.inp

# In[16]:


# Conditions used to write files
T = [300., 400., 500.] # Temperature in K
P = [1., 2., 3.] # Pressure in atm
Q = [10., 20., 30.] # Standard volumetric flow rate in cm3
abyv= [100., 50., 25.] # Catalyst surface area to reactor volume in 1/cm

ck_io.write_T_flow(filename='./outputs/T_flow.inp', T=T, P=P, Q=Q, abyv=abyv)


# As shown before, we can return T_flow as a string by omitting the filename.

# In[17]:


T_flow_str = ck_io.write_T_flow(T=T, P=P, Q=Q, abyv=abyv)
print(T_flow_str)


# ### Writing EAg.inp and EAs.inp

# In[18]:


# Convert T_flow inputs into list of dictionaries that can be used by write_EA.
# In the future, this will be replaced by a function
conditions = []
for T_i, P_i, Q_i, abyv_i in zip(T, P, Q, abyv):
    condition = {
        'T': T_i,
        'P': P_i,
        'Q': Q_i,
        'abyv': abyv}
    conditions.append(condition)

ck_io.write_EA(filename='./outputs/EAs.inp',
               reactions=reactions,
               write_gas_phase=False,
               act_method_name='get_GoRT_act',
               ads_act_method='get_HoRT_act',
               conditions=conditions)
ck_io.write_EA(filename='./outputs/EAg.inp',
               reactions=reactions,
               write_gas_phase=True,
               act_method_name='get_GoRT_act',
               conditions=conditions)


# Reminder that we use `act_method_name` as 'get_GoRT_act' for the [reason described above](#act_method_name_explanation).

# In[19]:


EAg_file = ck_io.write_EA(reactions=reactions,
                          write_gas_phase=False,
                          act_method_name='get_GoRT_act',
                          conditions=conditions)
print(EAg_file)


# In[20]:


EAs_file = ck_io.write_EA(reactions=reactions,
                          write_gas_phase=True,
                          act_method_name='get_GoRT_act',
                          conditions=conditions)
print(EAs_file)


# ### Writing tube_mole.inp

# In[21]:


import numpy as np

# Generating a list of conditions to input
mole_frac_conditions = []
for x_N2 in np.linspace(0., 0.25, 3):
    x_H2 = x_N2*3.
    x_NH3 = 1. - x_N2 - x_H2
    mole_fractions = {'N2': x_N2, 'H2': x_H2, 'NH3': x_NH3, 'RU(S)': 1.}
    mole_frac_conditions.append(mole_fractions)
    
# Write the tube_mole.inp file
ck_io.write_tube_mole(mole_frac_conditions=mole_frac_conditions, 
                      nasa_species=species, 
                      filename='./outputs/tube_mole.inp')


# Below is the output of the function.

# In[22]:


tube_mole_str = ck_io.write_tube_mole(mole_frac_conditions=mole_frac_conditions, 
                                      nasa_species=species)
print(tube_mole_str)

