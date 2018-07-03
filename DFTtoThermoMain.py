# -*- coding: utf-8 -*-
"""
Created on Fri Apr 07 09:39:17 2017

@author: wittregr
"""

import DFT_to_Thermochemistry as _thermo
import os
import platform

if platform.system() == 'Windows':
    Base_path = 'C:\Users\gerha\Documents\Python Scripts'
else:
    Base_path = '/home/1729'
Input = 'Input'
Output = 'Output'

'''
 Read reference species and calculate all thermodynamic quantities
'''

filepath = os.path.join(Base_path, Input, 'Reference_set_info.txt')
#filepath = os.path.join(Base_path, Input, 'Zacros_Species_Energy.txt')
[lines, dict] = _thermo.DFTFileRead(filepath)
T_ref = []
Path = os.path.join(Base_path, Input)
for s in lines[3:]:
    T_ref.append(_thermo.Reference(s.split('\t'), dict, Path, 298.15))

'''
 Read target species and calculate all thermodynamic quantities
'''

filepath = os.path.join(Base_path, 'Input', 'Tobe_Referenced2.txt')
[lines, dict] = _thermo.DFTFileRead(filepath)
T_target = []
for s in lines[3:]:
    T_target.append(_thermo.Target(s.split('\t'), dict, Base_path))

'''
Create basis set from reference molecules
'''
Basis = _thermo.Reference.BasisSet(T_ref)

'''
Read surface slab energy from input file
'''
filepath = os.path.join(Base_path, Input, 'Surface_energy_data.txt')
[lines, dict] = _thermo.DFTFileRead(filepath)
T_surface = []
for s in lines[3:]:
    T_surface.append(_thermo.Surface(s.split('\t'), dict))

'''
Apply reference molecule correct to thermodynamic data
'''
T_target = _thermo.Target.ReferenceDFT(T_target, T_surface, Basis)

'''
Output thermat file from calculated thermodynamic data
'''
T_target = _thermo.Target.CreateThermdat(T_target, Base_path, Output)

'''
Output Zacros input file from thermodynamic data
'''
pass

'''


Test output (To be removed in final version)


'''
print '\n\nHeat Capacities [cal/mol K] over a range of Temperatures [K]'
print '------------------------------------------------------------\n'
print 'Name             ',
for y in range(0, len(T_target[0].Cp_Range)):
    print ' %4d ' % (T_target[0].Cp_Range[y]),
print '\n',
print '---------------- ',
for y in range(0, len(T_target[0].Cp_Range)):
    print ' %s ' % ('----'),
print '\n',
for x in range(0, len(T_target)):
    print '%16s ' % (T_target[x].name),
    for y in range(0, len(T_target[x].Cp)):
        print ' %5.2f' % (T_target[x].Cp[y]),
    print '\n',

print '\n\nStandard Heats of Formation and Dispersion'
print '------------------------------------------\n'
print '                   hf(%6.2f)     edisp      S(%6.2f)' % (T_target[0].Tstp,
                                                            T_target[0].Tstp)
print 'Name               [kcal/mol]   [kcal/mol]  [cal/mol K]'
print '----------------   ----------   ----------  -----------'
for x in range(0, len(T_target)):
    print '%-16s     %7.3f      %6.3f       %6.3f' % (T_target[x].name,
                                                 T_target[x].hf_Tstp,
                                                 T_target[x].convedisp,
                                                 T_target[x].S_Tstp)
print '\n\n'
