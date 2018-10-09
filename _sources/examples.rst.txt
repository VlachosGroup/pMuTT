.. _examples:

Examples
========

*Ab-initio* Data to Empirical Data
----------------------------------
If you want a NASA polynomial (or other empirical form) and have *ab-initio* data. 
In this example, the data was manually inputted into a spreadsheet.
This example can be found under `PyMuTT.examples.VASP_to_thermdat.example1`_


.. figure:: flowchart_dft_to_empirical.png
   :alt: Flowchart of DFT data to empirical

1. **Read the ab-initio data.**
   Different thermodynamic models will require different inputs. 
   For example, the ``HarmonicThermo`` model only requires potentialenergy 
   and vibrational wavenumbers whereas ``IdealGasThermo`` requires the previously
   mentioned parameters and an atoms object, geometry, symmetry number, and spin.

.. code:: python

   from pprint import pprint
   from PyMuTT.io_.excel import read_excel

   species_data = read_excel(io=species_path)
   print('Input data imported from Excel')
   pprint(species_data)

One element of species_data is shown below. 
Note that ``read_excel`` returns a dictionary instead of any type of object. 
This makes it easy to initialize any type of object.::

   {'atoms': Atoms(symbols='OH2', pbc=True, cell=[20.0, 21.0, 22.0]),
     'elements': {'H': 2, 'O': 1, 'Pt': 0},
     'geometry': 'nonlinear',
     'name': 'H2O',
     'phase': 'G',
     'potentialenergy': -14.2209,
     'spin': 0.0,
     'symmetrynumber': 2.0,
     'thermo_model': <class 'PyMuTT.models.statmech.idealgasthermo.IdealGasThermo'>,
     'vib_energies': [0.47429336414391626,
                      0.460014128927786,
                      0.19619656143825398]}

2. **(Optional) Read the references data.**
   Reference data is needed to adjust DFT enthalpies to real-world enthalpies. 
   It is typical to use gas-phase molecules since high-accuracy experimental data is readily available. 
   The references are calculated for each atom type so the number of references should be greater than 
   or equal to the number of elements. See the :ref:`referencing` for more information.

.. code:: python

   from PyMuTT.models.empirical import BaseThermo
   from PyMuTT.models.empirical.references import References
   
   refs_input = read_excel(io=refs_path)
   print('Reference data imported from Excel')
   pprint(refs_input)
   refs = References([BaseThermo(**ref_input) for ref_input in refs_input])
   print('Reference data')
   pprint(refs)

Similar to above, the result of ``pprint(refs_input)`` for one reference is shown below.::

   {'HoRT_ref': 0.0,
     'T_ref': 298,
     'atoms': Atoms(symbols='H2', pbc=True, cell=[20.0, 21.0, 22.0]),
     'elements': {'H': 2, 'O': 0},
     'geometry': 'linear',
     'name': 'H2',
     'phase': 'G',
     'potentialenergy': -6.7598,
     'spin': 0,
     'symmetrynumber': 2,
     'thermo_model': <class 'PyMuTT.models.statmech.idealgasthermo.IdealGasThermo'>,
     'vib_energies': [0.5338981843116086]},

The result for ``pprint(refs)`` for the same reference is shown below.::

   BaseThermo object for Name: H2
           phase: G
           elements: {'H': 2, 'O': 0}
           T_ref: 298
           references: None
           notes: None
           thermo_model: <PyMuTT.models.statmech.idealgasthermo.IdealGasThermo object at 0x0000021043554518>
           HoRT_dft: -249.34037553721055
           HoRT_ref: 0.0

3. **Fit to the empirical form (e.g. Nasa polynomials)**

.. code:: python

   from PyMuTT.models.empirical.nasa import Nasa
   from PyMuTT import constants as c

   T_low = 200.
   T_high = 1100.
   species = [Nasa.from_statmech(references=refs, T_low=T_low, T_high=T_high, T_ref=c.T0('K'), **specie_data) for specie_data in species_data]

4. **Write to the appropriate file (e.g. thermdat)**

.. code:: python

   from PyMuTT.io_.thermdat import write_thermdat
   
   write_thermdat(nasa_species=species, filename='thermdat', write_date=True)

The thermdat file produced is shown below.::

   THERMO ALL
          100       500      1500
   H2O             20180824H   2O   1          G200.0     1100.0    493.9         1
    3.65264072E+00 1.06108515E-03 3.83455580E-08 3.84923664E-10-2.13953966E-13    2
   -3.02204928E+04 1.60236266E+00 3.99524709E+00 5.18551442E-04-5.53026360E-06    3
    1.85895538E-08-1.55138452E-11-3.02807840E+04-7.89384507E-02                   4
   H2              20180824H   2               G200.0     1100.0    677.6         1
    3.47458128E+00 3.33323856E-04-1.07413869E-06 1.23229932E-09-3.87505553E-13    2
   -1.04411622E+03-4.17284186E+00 3.50368711E+00-5.09520583E-05 2.59076698E-07    3
   -5.76347123E-10 4.75276554E-13-1.04320885E+03-4.26615761E+00                   4
   O2              20180824O   2               G200.0     1100.0    604.1         1
    3.81428705E+00-2.43962827E-03 5.76514137E-06-4.29570548E-09 1.11527700E-12    2
   -5.09263346E+03 3.47144761E+00 3.42023428E+00 1.11878623E-03-5.58175407E-06    3
    1.12396705E-08-6.69234827E-12-5.06166125E+03 5.03398906E+00                   4
   MO(S)           20180824O   1Pt  1          S200.0     1100.0    604.1         1
    7.58148263E-01 6.93537690E-03-9.09332498E-06 5.64388225E-09-1.36222297E-12    2
    5.46713976E+04-4.75898748E+00-1.85092759E+00 2.48466946E-02-5.62021000E-05    3
    6.17887348E-08-2.68761443E-11 5.49823404E+04 6.44810333E+00                   4
   MO(B)           20180824Pt  1               S200.0     1100.0    604.1         1
    7.58148263E-01 6.93537690E-03-9.09332498E-06 5.64388225E-09-1.36222297E-12    2
    7.87256857E+02-4.75898748E+00-1.85092759E+00 2.48466946E-02-5.62021000E-05    3
    6.17887348E-08-2.68761443E-11 1.09819962E+03 6.44810333E+00                   4
   V-MO(S)         20180824Pt  1               S200.0     1100.0    659.2         1
    0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
    8.12316548E+04 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    3
    0.00000000E+00 0.00000000E+00 8.12316548E+04 0.00000000E+00                   4
   MO_bulk(S)      20180824O   1Pt  1          S200.0     1100.0    604.1         1
    7.58148263E-01 6.93537690E-03-9.09332498E-06 5.64388225E-09-1.36222297E-12    2
    5.46713976E+04-4.75898748E+00-1.85092759E+00 2.48466946E-02-5.62021000E-05    3
    6.17887348E-08-2.68761443E-11 5.49823404E+04 6.44810333E+00                   4
   MO_bulk(B)      20180824Pt  1               S200.0     1100.0    604.1         1
    7.58148263E-01 6.93537690E-03-9.09332498E-06 5.64388225E-09-1.36222297E-12    2
    7.87256857E+02-4.75898748E+00-1.85092759E+00 2.48466946E-02-5.62021000E-05    3
    6.17887348E-08-2.68761443E-11 1.09819962E+03 6.44810333E+00                   4
   V-MO_bulk(S)    20180824Pt  1               S200.0     1100.0    659.2         1
    0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
    8.12316548E+04 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    3
    0.00000000E+00 0.00000000E+00 8.12316548E+04 0.00000000E+00                   4
   END

5. **(Optional) Compare the statistical mechanic model to empirical fit.**

.. code:: python

   import matplotlib.pyplot as plt

   for specie in species:
       specie.plot_thermo_model_and_empirical(Cp_units='J/mol/K', H_units='kJ/mol', S_units='J/mol/K', G_units='kJ/mol')
   plt.show()

.. figure:: H2_VASP_to_thermdat_example.png

Experimental Data to Empirical Data
-----------------------------------
If you want a NASA polynomial (or other empirical form) and have the heat capacity at various temperatures, enthalpy of formation and entropy of formation. 
This example can be found under `PyMuTT.examples.Expt_data_to_thermdat`_.

1. Import the input data. In this example, the data is `methanol gas-phase data from NIST`_. This includes:
- heat capacity and corresponding temperatures
- enthalpy of formation 
- standard entropy

.. code:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from PyMuTT import constants as c
   from PyMuTT.models.empirical.nasa import Nasa

   # Gas phase heat capacity data (in J/mol/K) for CH3OH from NIST
   # https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Units=SI&Mask=1#Thermo-Gas
   T = np.array([50, 100, 150, 200, 273.15, 298.15, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2250, 2500, 2750, 3000])
   Cp = np.array([34.00, 36.95, 38.64, 39.71, 42.59, 44.06, 44.17, 51.63, 59.7, 67.19, 73.86, 79.76, 84.95, 89.54, 93.57, 97.12, 100.24, 102.98, 105.4, 110.2, 113.8, 116.5, 118.6, 120, 121])
   CpoR = Cp/c.R('J/mol/K')

   #Enthalpy of Formation for CH3OH (in kJ/mol) from NIST
   T_ref = c.T0('K')
   H_ref = -205.
   HoRT_ref = H_ref/c.R('kJ/mol/K')/T_ref
   #Standard molar entropy (in J/mol/K) from Wikipedia, https://en.wikipedia.org/wiki/Methanol_(data_page)
   S_ref = 239.9
   SoR_ref = S_ref/c.R('J/mol/K')

2. Fit to a NASA polynomial

.. code:: python

   #Input the experimental data and fitting to a NASA polynomial
   CH3OH_nasa = Nasa.from_data(name ='CH3OH', T=T, CpoR=CpoR, T_ref=T_ref, HoRT_ref=HoRT_ref, SoR_ref=SoR_ref)

3. (Optional) Compare the Nasa polynomial to the experimental data

.. code:: python

   #Units to plot the figure
   Cp_units = 'J/mol/K'
   H_units = 'kJ/mol'
   S_units = 'J/mol/K'
   G_units = 'kJ/mol'

   #Compare the Nasa polynomial to the input data
   fig, axes = CH3OH_nasa.plot_empirical(Cp_units=Cp_units, H_units=H_units, S_units=S_units, G_units=G_units)
   axes[0].plot(T, Cp, 'ko')
   axes[1].plot(T_ref, H_ref, 'ko')
   axes[2].plot(T_ref, S_ref, 'ko')
   axes[3].plot(T_ref, H_ref - T_ref * S_ref * c.convert_unit(from_='J', to='kJ'), 'ko')
   plt.show()

.. figure:: expt_data_to_thermdat_example.png

Read Nasa from Thermdat
-----------------------

If you want to read a thermdat file directly and convert it into a list of Nasa objects.
This example can be found under `PyMuTT.examples.read_nasa_from_thermdat`_.

1. Read the data from the thermdat file. In this example, a ``thermdat`` file is located in the same folder as the script.

.. code:: python

   from PyMuTT.io_.thermdat import read_thermdat
   
   species = read_thermdat('thermdat')

2. (Optional) Print out the parameters to ensure it was read properly

.. code:: python

   from pprint import pprint
   pprint(species)
   
A snippet of the output is shown below:

::

   Nasa object for Name: O2
          phase: G
          elements: {'O': 2}
          T_ref: 298.15
          references: None
          notes: TPIS89
          thermo_model: None
          HoRT_dft: None
          HoRT_ref: None
          a_low: [ 3.78245636e+00 -2.99673416e-03  9.84730201e-06 -9.68129509e-09 3.24372837e-12 -1.06394356e+03  3.65767573e+00]
          a_high: [ 3.28253784e+00  1.48308754e-03 -7.57966669e-07  2.09470555e-10 -2.16717794e-14 -1.08845772e+03  5.45323129e+00]
          T_low: 200.0
          T_high: 3500.0
          T_mid: 1000.0,

3. (Optional) Plot each species to ensure reasonable values. The below example shows how to print out a single species.

.. code:: python

   species[1].plot_empirical()
   plt.show()
   

The output for O2 is shown below:

.. figure:: read_nasa_from_thermdat_example_O2.png

.. _`PyMuTT.examples.VASP_to_thermdat.example1`: https://github.com/VlachosGroup/PyMuTT/tree/master/examples/VASP_to_thermdat/example1
.. _`PyMuTT.examples.Expt_data_to_thermdat`: https://github.com/VlachosGroup/PyMuTT/blob/master/examples/Expt_data_to_thermdat/Expt_data_to_thermdat.py
.. _`methanol gas-phase data from NIST`: https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Units=SI&Mask=1#Thermo-Gas
.. _`PyMuTT.examples.read_nasa_from_thermdat`: https://github.com/VlachosGroup/PyMuTT/tree/master/examples/read_nasa_from_thermdat