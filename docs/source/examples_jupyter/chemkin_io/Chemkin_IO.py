
# coding: utf-8

# # Chemkin Input and Output
# This notebook describes pMuTT's functionality to read and write Chemkin files. We will use the NH3 formation mechanism as a case study.
# 
# ## Topics Covered
# - Read species *ab-initio* data, reactions, and catalyst sites from a spreadsheet
# - Write the thermdat, gas.inp, surf.inp, T_flow.inp, EAg.inp, EAs.inp, tube_mole.inp files

# ## Input Spreadsheet
# All the data will be imported from the [`./inputs/NH3_Input_data.xlsx`](https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/chemkin_io/inputs/NH3_Input_Data.xlsx) file. There are four sheets:
# 1. `cat_sites` contains catalyst site properties for the adsorbed species
# 2. `refs` contains *ab-initio* and experimental data for a handful of gas species to calculate references
# 3. `species` contains *ab-initio* data for each specie
# 4. `reactions` contains elementary steps
# 
# The contents are displayed below

# **Catalytic Sites**
# 
# | name   | site_density | density | bulk_specie |
# |--------|--------------|---------|-------------|
# | RU0001 | 2.1671E-09   | 12.2    | RU(B)       |

# **References**
# 
# | name | elements.N | elements.H | T_ref  | HoRT_ref     | potentialenergy | symmetrynumber | statmech_model | atoms         | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber |
# |------|------------|------------|--------|--------------|-----------------|----------------|----------------|---------------|----------------|----------------|----------------|----------------|
# | N2   | 2          | 0          | 298.15 | 0            | -16.63          | 2              | IdealGas       | ./N2/CONTCAR  | 2744           |
# | NH3  | 1          | 3          | 298.15 | -18.38025311 | -19.54          | 3              | IdealGas       | ./NH3/CONTCAR | 3534           | 3464           | 1765           | 1139           |
# | H2   | 0          | 2          | 298.15 | 0            | -6.7700         | 2              | IdealGas       | ./H2/CONTCAR  | 4342           |

# **Species**
# 
# | name | elements.N | elements.H | phase | statmech_model | symmetrynumber | atoms | potentialenergy | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber | vib_wavenumber |
# |------------|------------|------------|-------|----------------|----------------|---------------|-----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|
# | N2 | 2 | 0 | G | IdealGas | 2 | ./N2/CONTCAR | -16.63 | 2744.00 |
# | NH3 | 1 | 3 | G | IdealGas | 3 | ./NH3/CONTCAR | -19.54 | 3534.00 | 3464.00 | 1765.00 | 1139.00 |
# | H2 | 0 | 2 | G | IdealGas | 2 | ./H2/CONTCAR | -6.77 | 4342.00 |
# | N2(S) | 2 | 0 | S | Harmonic |  |  | -17.24 | 2197.19 | 360.42 | 347.34 | 335.67 | 62.08 | 32.18 |
# | N(S) | 1 | 0 | S | Harmonic |  |  | -9.34 | 549.11 | 538.56 | 504.32 | 475.81 | 459.08 | 410.02 |
# | H(S) | 0 | 1 | S | Harmonic |  |  | -4.00 | 1003.51 | 625.55 | 616.29 |
# | NH3(S) | 1 | 3 | S | Harmonic |  |  | -20.43 | 3491.09 | 3488.82 | 3364.52 | 1583.52 | 1582.07 | 1124.22 | 570.21 | 567.22 | 333.09 | 122.86 | 83.83 | 70.63 |
# | NH2(S) | 1 | 2 | S | Harmonic |  |  | -16.59 | 3469.30 | 3381.05 | 1503.02 | 698.87 | 625.60 | 615.94 | 475.13 | 298.12 | 153.25 |
# | NH(S) | 1 | 1 | S | Harmonic |  |  | -13.21 | 3403.13 | 718.18 | 710.58 | 528.53 | 415.20 | 410.13 |
# | TS1_NH3(S) | 1 | 3 | S | Harmonic |  |  | -19.53 | 3453.41 | 3355.67 | 1723.85 | 1487.95 | 959.15 | 888.95 | 594.09 | 428.43 | 227.03 | 206.05 | 142.14 |  |
# | TS2_NH2(S) | 1 | 2 | S | Harmonic |  |  | -16.09 | 3426.44 | 1293.72 | 922.83 | 660.97 | 525.60 | 496.84 | 330.67 | 290.28 |
# | TS3_NH(S) | 1 | 1 | S | Harmonic |  |  | -12.91 | 1201.60 | 491.57 | 462.02 | 402.16 | 242.14 |
# | TS4_N2(S) | 2 | 0 | S | Harmonic |  |  | -16.54 | 485.61 | 392.98 | 386.19 | 280.94 | 168.43 |
# | RU(S) | 0 | 0 | S | Placeholder |
# | RU(B) | 0 | 0 | S | Placeholder |

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

# ## Reading data
# Before we can initialize our species, we need the catalytic sites and the references.
# 
# ### Reading Catalytic Sites

# In[1]:


import os
from pprint import pprint
from pathlib import Path
from pMuTT.io.excel import read_excel
from pMuTT.chemkin import CatSite

# Find the location of Jupyter notebook
# Note that normally Python scripts have a __file__ variable but Jupyter notebook doesn't.
# Using pathlib can overcome this limiation
notebook_path = Path().resolve()
os.chdir(notebook_path)
excel_path = './inputs/NH3_Input_data.xlsx'

cat_site_data = read_excel(io=excel_path, sheet_name='cat_sites')[0]
cat_site = CatSite(**cat_site_data)

# Print the properties of the catalyst site
pprint(cat_site.to_dict())


# ### Reading reference species

# In[2]:


from pMuTT.empirical.references import Reference, References

references_data = read_excel(io=excel_path, sheet_name='refs')

# Convert data to Reference objects and put them in a list
refs_list = [Reference(**ref_data) for ref_data in references_data]

# Assign the Reference objects to a References object so offsets can be calculated
refs = References(references=refs_list)

# Print out the offsets calculated
print(refs.offset)


# ### Reading species

# In[3]:


from pMuTT.empirical.nasa import Nasa

# Range of data to fit the Nasa polynomials
T_low = 298. # K
T_high = 800. # K

species_data = read_excel(io=excel_path, sheet_name='species')
species = []
for specie_data in species_data:
    specie = Nasa.from_statmech(T_low=T_low, T_high=T_high, references=refs, **specie_data)
    # If the species is a surface species, assign the catalyst site specified above
    if specie.phase.lower() == 's':
        specie.cat_site = cat_site
        specie.n_sites = 1
    species.append(specie)


# The warning above is typical when empirical objects are fitting to `StatMech` objects with the `placeholder` preset.

# ### Reading reactions

# In[4]:


from pMuTT import pMuTT_list_to_dict
from pMuTT.reaction import ChemkinReaction, Reactions

# Convert list of Nasa polynomials to dictionary of Nasa polynomials
species_dict = pMuTT_list_to_dict(species)

reactions_data = read_excel(io=excel_path, sheet_name='reactions')
reactions_list = []
for reaction_data in reactions_data:
    reaction = ChemkinReaction.from_string(species=species_dict, **reaction_data)
    reactions_list.append(reaction)
reactions = Reactions(reactions=reactions_list)


# ## Writing Chemkin files
# Now that we have all the required objects, we can write the output files. All outputs can be found in the [./outputs folder](https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/chemkin_io/outputs).

# ### Writing thermdat

# In[5]:


from pMuTT.io.thermdat import write_thermdat

write_thermdat(filename='./outputs/thermdat', nasa_species=species)


# The thermdat file written in shown below.
# ```
# THERMO ALL
#        100       500      1500
# N2              20190425N   2               G298.0     800.0     523.4         1
#  3.75207309E+00-1.30159027E-03 1.83848785E-06-1.17413430E-10-3.65391117E-13    2
# -1.67464905E+02 2.02487460E+00 3.41792195E+00 8.91688129E-04-3.46864753E-06    3
#  5.44691795E-09-2.46804017E-12-1.27217456E+02 3.46925670E+00                   4
# NH3             20190425N   1H   3          G298.0     800.0     533.6         1
#  3.29427463E+00 2.57481517E-03 5.20883059E-07-1.00678313E-09 3.59356180E-13    2
# -8.39739662E+03 3.59518301E+00 4.57124680E+00-6.89724277E-03 2.70972941E-05    3
# -3.44331645E-08 1.62568430E-11-8.53630325E+03-1.78236359E+00                   4
# H2              20190425H   2               G298.0     800.0     574.6         1
#  3.37903251E+00 7.92301907E-04-1.88662952E-06 1.85974220E-09-5.68447254E-13    2
#  1.70231785E+03-3.75384879E+00 3.50717925E+00-8.50397135E-05 3.79570576E-07    3
# -7.57943298E-10 5.72423558E-13 1.68725520E+03-4.30359617E+00                   4
# N2(S)           20190425N   2               S298.0     800.0     482.4         1
#  3.42605836E+00 5.03710708E-03-6.92689694E-06 6.03246113E-09-2.18862675E-12    2
# -7.10922368E+03-1.39236755E+01 4.40553502E-01 2.97759950E-02-8.47698148E-05    3
#  1.16195177E-07-6.12938589E-11-6.81726811E+03-1.67561821E+00                   4
# N(S)            20190425N   1               S298.0     800.0     502.9         1
# -8.29107669E-01 2.64469843E-02-4.38832016E-05 3.47600781E-08-1.07819310E-11    2
# -1.10025913E+04 5.84749603E-01-5.23174593E+00 6.11537184E-02-1.47762309E-04    3
#  1.74618161E-07-8.21997526E-11-1.05502378E+04 1.88654280E+01                   4
# H(S)            20190425H   1               S298.0     800.0     461.9         1
# -2.11302416E+00 1.72572991E-02-2.60046478E-05 1.92388088E-08-5.68363784E-12    2
# -6.07693323E+03 8.35415551E+00-1.41307531E+00 1.00687718E-02 1.09758988E-06    3
# -2.55216914E-08 2.17566538E-11-6.12970520E+03 5.64702945E+00                   4
# NH3(S)          20190425N   1H   3          S298.0     800.0     461.9         1
#  1.58312977E+00 1.57753291E-02-1.71587153E-05 1.15414255E-08-3.21958978E-12    2
# -1.42569335E+04-6.35113154E+00 1.16759533E+00 2.05732170E-02-3.66958772E-05    3
#  4.56112325E-08-2.49628744E-11-1.42313136E+04-4.80541265E+00                   4
# NH2(S)          20190425N   1H   2          S298.0     800.0     523.4         1
# -6.86530719E-01 2.43298826E-02-3.68705157E-05 2.88483783E-08-8.86963660E-12    2
# -1.19475429E+04 1.03190122E+00-2.56305506E+00 3.84726874E-02-7.73503504E-05    3
#  8.09903199E-08-3.43595838E-11-1.17458710E+04 8.90805309E+00                   4
# NH(S)           20190425N   1H   1          S298.0     800.0     543.9         1
# -1.33660645E+00 2.33106600E-02-3.77180253E-05 2.97486446E-08-9.12157176E-12    2
# -1.48903408E+04 3.58875312E+00-3.73687435E+00 4.05226889E-02-8.44226964E-05    3
#  8.65932632E-08-3.53036247E-11-1.46202194E+04 1.37781593E+01                   4
# TS1_NH3(S)      20190425N   1H   3          S298.0     800.0     451.7         1
#  5.21221040E-02 2.11615812E-02-2.49542203E-05 1.66204937E-08-4.59897868E-12    2
# -5.86628985E+03-1.43112446E+00 4.96418515E-02 2.15674671E-02-2.75604959E-05    3
#  2.23225323E-08-8.79255935E-12-5.87008224E+03-1.46494399E+00                   4
# TS2_NH2(S)      20190425N   1H   2          S298.0     800.0     554.1         1
# -1.48538501E+00 2.67287608E-02-3.86606417E-05 2.81558194E-08-8.15140597E-12    2
# -8.43975871E+03 4.23478568E+00-2.58881424E+00 3.43781302E-02-5.87419101E-05    3
#  5.18321984E-08-1.87340797E-11-8.31127804E+03 8.95651586E+00                   4
# TS3_NH(S)       20190425N   1H   1          S298.0     800.0     482.4         1
# -1.24242592E-01 1.65859456E-02-2.43704704E-05 1.77594067E-08-5.19878424E-12    2
# -1.39060398E+04-1.18117850E+00-1.80701218E+00 3.07479132E-02-6.96816870E-05    3
#  8.30020912E-08-4.08170225E-11-1.37437409E+04 5.69757432E+00                   4
# TS4_N2(S)       20190425N   2               S298.0     800.0     482.4         1
#  1.52431744E+00 1.40315985E-02-2.40821843E-05 1.96281604E-08-6.24114074E-12    2
# -1.24945441E+02-8.60931129E+00-1.75107618E+00 4.12189598E-02-1.09828935E-04    3
#  1.41320046E-07-7.17404084E-11 1.95040067E+02 4.82355594E+00                   4
# RU(S)           20190425                    S298.0     800.0     554.1         1
#  0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
#  0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    3
#  0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00                   4
# RU(B)           20190425                    S298.0     800.0     554.1         1
#  0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
#  0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    3
#  0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00                   4
# END
# ```

# ### Writing gas.inp and surf.inp

# In[6]:


from pMuTT.io import chemkin as ck_io

ck_io.write_gas(filename='./outputs/gas.inp', nasa_species=species, reactions=reactions, act_method_name='get_GoRT_act')
ck_io.write_surf(filename='./outputs/surf.inp', reactions=reactions, act_method_name='get_GoRT_act')


# <a id='act_method_name_explanation'></a>
# Note that `act_method_name` is 'get_GoRT_act'. We use this formalism here since we do not include entropic effects in the preexponential factor.
# 
# The gas.inp file written is shown below. Note there are no gas-phase reactions.
# 
# ```
# !File generated by pMuTT on 2019-04-25 16:05:46.954918
# !Elements present in gas and surface species
# ELEMENTS
# H
# N
# END
# 
# !Gas-phase species
# SPECIES
# N2
# NH3
# H2
# END
# 
# !Gas-phase reactions. The rate constant expression is:
# !k = kb/h * (T)^beta * exp(-Ea/RT)
# !Each line has 4 columns:
# !- Reaction reactants and products separated by <=>
# !- Preexponential factor, kb/h
# !- Beta (power to raise T in rate constant expression)
# !- Ea (Activation Energy or Gibbs energy of activation in kcal/mol
# REACTIONS
# END
# ```
# 
# The surf.inp file written in shown below.
# ```
# !File generated by pMuTT on 2019-04-25 16:05:46.977586
# !Surface species
# !Each catalyst site has the following format:
# !SITE/[Site name]/      SDEN/[Site density in mol/cm2]/
# ![Adsorbate Name]/[# of Sites occupied]/ (for every adsorbate)
# !BULK [Bulk name]/[Bulk density in g/cm3]
# SITE/RU0001/       SDEN/2.16710E-09/
# 
#   RU(S)/1/
#   H(S)/1/
#   N2(S)/1/
#   NH3(S)/1/
#   NH2(S)/1/
#   NH(S)/1/
#   N(S)/1/
# 
# BULK RU(B)/12.2/
# END
# 
# !Surface-phase reactions.
# !The reaction line has the following format:
# !REACTIONS  MW[ON/OFF]   [Ea units]
# !where MW stands for Motz-Wise corrections and if the Ea
# !units are left blank, then the activation energy should be in cal/mol
# !The rate constant expression is:
# !k = kb/h/site_den^(n-1) * (T)^beta * exp(-Ea/RT)
# !where site_den is the site density and is the number of surface species (including empty sites)
# !Each line has 4 columns:
# !- Reaction reactants and products separated by =
# !- Preexponential factor, kb/h/site_den^(n-1), or 
# !  sticking coefficient if adsorption reaction
# !- Beta (power to raise T in rate constant expression)
# !- Ea (Activation Energy or Gibbs energy of activation in specified units
# !Adsorption reactions can be represented using the STICK keyword
# REACTIONS  MWON   
# H2+2RU(S)=2H(S)+2RU(B)           1.000E+00   1.000E+00   0.000E+00
# STICK
# N2+RU(S)=N2(S)+RU(B)             1.000E+00   1.000E+00   0.000E+00
# STICK
# NH3+RU(S)=NH3(S)+RU(B)           1.000E+00   1.000E+00   0.000E+00
# STICK
# NH3(S)+RU(S)=NH2(S)+H(S)+RU(B)   9.615E+18   1.000E+00   2.972E+01
# NH2(S)+RU(S)=NH(S)+H(S)+RU(B)    9.615E+18   1.000E+00   1.199E+01
# NH(S)+RU(S)=N(S)+H(S)+RU(B)      9.615E+18   1.000E+00   3.203E+00
# 2N(S)+RU(B)=N2(S)+RU(S)          9.615E+18   1.000E+00   7.317E+01
# END
# ```

# ### Writing T_flow.inp

# In[7]:


# Conditions used to write files
T = [300., 400., 500.] # Temperature in K
P = [1., 2., 3.] # Pressure in atm
Q = [10., 20., 30.] # Standard volumetric flow rate in cm3
abyv= [100., 50., 25.] # Catalyst surface area to reactor volume in 1/cm

ck_io.write_T_flow(filename='./outputs/T_flow.inp', T=T, P=P, Q=Q, abyv=abyv)


# The T_flow.inp file written is shown below.
# 
# ```
# !File generated by pMuTT on 2019-04-25 16:05:47.025434
# !Conditions for each reaction run
# !Only used when MultiInput in tube.inp is set to "T"
# !T[K]      P[atm]     Q[cm3/s]   abyv[cm-1]  Run #
# 3.000E+02  1.000E+00  1.000E+01  1.000E+02  !1  
# 4.000E+02  2.000E+00  2.000E+01  5.000E+01  !2  
# 5.000E+02  3.000E+00  3.000E+01  2.500E+01  !3  
# EOF
# ```

# ### Writing EAg.inp and EAs.inp

# In[8]:


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
               conditions=conditions)
ck_io.write_EA(filename='./outputs/EAg.inp',
               reactions=reactions,
               write_gas_phase=True,
               act_method_name='get_GoRT_act',
               conditions=conditions)


# Reminder that we use `act_method_name` as 'get_GoRT_act' for the [reason described above](#act_method_name_explanation).
# 
# The EAg.inp file is shown below.
# ```
# !File generated by pMuTT on 2019-04-25 16:05:47.088278
# !The first line is the number of reactions. Subsequent lines follow the format
# !of rxn (from surf.out) followed by the EA/RT value at each run condition.
# !There may be one slight deviation from surf.out: any repeated species should
# !be included in the reaction string with a stoichiometric coefficient equal to
# !the number of times the species appears in the reaction. If not using
# !MultiInput, then only the first value is used.
#   0  !Number of reactions
# !           1          2          3
# EOF
# ```
# 
# The EAs.inp file is shown below.
# ```
# !File generated by pMuTT on 2019-04-25 16:05:47.043381
# !The first line is the number of reactions. Subsequent lines follow the format
# !of rxn (from surf.out) followed by the EA/RT value at each run condition.
# !There may be one slight deviation from surf.out: any repeated species should
# !be included in the reaction string with a stoichiometric coefficient equal to
# !the number of times the species appears in the reaction. If not using
# !MultiInput, then only the first value is used.
#   7  !Number of reactions
# !                                           1          2          3
# H2+2RU(S)<=>2H(S)+2RU(B)             0.00E+00   0.00E+00   0.00E+00
# N2+RU(S)<=>N2(S)+RU(B)               0.00E+00   0.00E+00   2.49E+00
# NH3+RU(S)<=>NH3(S)+RU(B)             0.00E+00   1.61E+00   4.47E+00
# NH3(S)+RU(S)<=>NH2(S)+H(S)+RU(B)     2.95E+01   2.28E+01   1.88E+01
# NH2(S)+RU(S)<=>NH(S)+H(S)+RU(B)      1.19E+01   9.12E+00   7.46E+00
# NH(S)+RU(S)<=>N(S)+H(S)+RU(B)        3.18E+00   2.23E+00   1.65E+00
# 2N(S)+RU(B)<=>N2(S)+RU(S)            7.27E+01   5.49E+01   4.45E+01
# EOF
# ```

# ### Writing tube_mole.inp

# In[9]:


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


# The tube_mole.inp file is shown below.
# 
# ```
# !File generated by pMuTT on 2019-04-25 16:05:47.114596
# !Specify the 'species/phase/' pair /(in quotes!)/ and the associated
# !composition values. If the composition does not sum to 1 for each phase or
# !site type, it will be renormalized to 1. At the end of a calculation, a
# !file containing the complete composition and mass flux (the last entry) will
# !be generated. This file's format is completely compatible with the current
# !input file and can be used to restart that calculation.
# 0       itube_restart -- will be >0 if a restart file is used or 0 for the first run
# 4      Number of nonzero species
# !                       1       2       3
# 'N2/GAS/'           0.000   0.125   0.250
# 'NH3/GAS/'          1.000   0.500   0.000
# 'H2/GAS/'           0.000   0.375   0.750
# 'RU(S)/RU0001/'     1.000   1.000   1.000
# EOF
# ```
