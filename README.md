# Thermochemistry
This repository contains a Python library for Thermochemistry developed by the Vlachos Research Group at the University of Delaware.

## Implemented features
1. DFT-generated data can be input using Excel and special column header names.
2. Ideal Gas and Harmonic Thermo statistical thermodynamic models have been implemented. Framework to create new models is straightforward (classes only need to: state explicitly attributes required in the ```__init__``` method, contain the methods ```get_CpoR```, ```get_HoRT```, ```get_SoR```, and ```get_GoRT```).
3. Thermdat empirical model has been implemented. Given Cp as a function of temperature, a reference temperature and associated enthalpy and entropy, NASA polynomials can be generated.
4. Framework to convert statistical thermodynamic model to empirical models implemented. See [Thermochemistry.examples.VASP_to_thermdat](https://github.com/VlachosGroup/Thermochemistry/tree/master/examples/VASP_to_thermdat).

## Planned features
1. Automatically read DFT-generated data directly from:
   - VASP
   - Gaussian
2. More complex statistical thermodynamic models like the Rigid Rotor model
3. Read and write to other empirical formats formats, including:
   - Shomate polynomials

## Developers
- Gerhard Wittreich, P.E. (wittregr@udel.edu)
- Jonathan Lym (jlym@udel.edu)

## Dependencies
- Python3
- [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/): Used for I/O operations and to calculate thermodynamic properties
- [Numpy](http://www.numpy.org/): Used for vector and matrix operations
- [Pandas](https://pandas.pydata.org/): Used to import data from Excel files
- [SciPy](https://www.scipy.org/): Used for fitting heat capacities.
- [Matplotlib](https://matplotlib.org/): Used for plotting thermodynamic data

## Getting Started
1. Download the repository to your local machine
2. Add to parent folder to PYTHONPATH

## Outline to convert DFT-generated data to empirical models
*Code can be found in [Thermochemistry.examples.VASP_to_thermdat](https://github.com/VlachosGroup/Thermochemistry/tree/master/examples/VASP_to_thermdat)*
1. **Read the Reference data.** Reference data is needed to adjust DFT enthalpies to real-world enthalpies. It is typical to use gas-phase molecules since experimental data is readily available. The number of references should be greater than or equal to the number of elements. 

For this example, our references are specified in an Excel file. The experimental data is specified using the reference dimensionless enthalpy of formation field, ```HoRT_ref```. The computational data is specified using the fields: ```potentialenergy``` (electronic energy in eV), ```vib_wavenumber``` (vibrational frequencies in 1/cm), ```geometry``` (linear, nonlinear, monatomic), ```atoms``` (the location of the CONTCAR files), and the ```thermo_model``` to use (IdealGasThermo). ```BaseThermo``` is the parent class of empirical classes. For referencing, ```BaseThermo``` is sufficient but any empirical model could have been used.
```python
from Thermochemistry.io_.excel import read_excel
from Thermochemistry.models.empirical import BaseThermo
from Thermochemistry.models.empirical.references import References

refs_path = './references.xlsx'
refs_input = read_excel(io=refs_path)
refs = References([BaseThermo(**ref_input) for ref_input in refs_input])
```

2. **Use the references to calculate the offset between DFT data and experimental data.** The offset is calculated for each element. Therefore, the element composition should be specified in the Excel file. The elemental composition can be specified using the ```formula``` keyword where the formula unit can be typed (e.g. H2O) or by using the ```elements``` keyword where each element has a column indicating the number in a formula unit. [Read below if you're interested on how referencing is done](#referencing).
```python
refs.calc_offset()
```

3. **Read the DFT-generated data for the interested species.** This is similar to Step 1 but the ```T_ref``` and ```HoRT_ref``` fields are not required.
```python
thermdats_data = read_excel(io=thermdats_in_path)
thermdats = [Thermdat(**thermdat_data) for thermdat_data in thermdats_data]
```

4. **Calculate the reference enthalpy using the offset calculated in Step 2 and fit to the desired empirical formula.**
```python
from Thermochemistry import constants as c
from Thermochemistry.models.empirical.thermdat import Thermdat

for thermdat_specie in thermdats:
	thermdat_specie.T_ref = c.T0('K')
	thermdat_specie.HoRT_ref = thermdat_specie.thermo_model.get_HoRT(Ts=c.T0('K')) + refs.get_specie_offset(thermdat_specie.elements)
	thermdat_specie.calc_nasa(T_low=T_low, T_high=T_high, T_ref=c.T0('K'))
```

5. **Save the coefficients of the empirical form.**
```python
from Thermochemistry.io_.thermdat import write_thermdat

write_thermdat(thermdats=thermdats, filename=thermdats_out_path)
```

## Referencing
Enthalpies calculated using VASP (and some other computational methods) have different references than standard references (i.e. the enthalpy of formation of pure substances, like O<sub>2</sub> or Pt, is not necessarily zero). This difference makes it difficult to ensure thermodynamic consistency for our mechanisms since we may be mixing experimental gas thermodynamics with computational surface thermodynamics. In order to make the references consistent, we find a correction factor for each element by solving the equation:

![Eq1](README_Eq1.gif)

where M is the number of reference species, N is the number of elements, H<sup>expt</sup> is the experimental standard enthalpies, H<sup>DFT</sup> is the standard enthalpies calculated using DFT, x is a matrix that describes the composition of the references (each row represents a specie, each column represents an element), and θ is the correction for each element.

The equation can be solved using a Least Squares approach. The correction factor can then be added to subsequent species calculated through DFT to ensure consistent references.