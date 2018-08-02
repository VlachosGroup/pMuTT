# **Py**thon **Mu**ltiscale **T**hermochemistry **T**oolbox (PyMuTT)
This repository contains a Python library for Thermochemistry developed by the Vlachos Research Group at the University of Delaware. This code was originally developed to convert *ab-initio* data from DFT to observable thermodynamic properties such as heat capacity, enthalpy, entropy, and Gibbs energy. These properties can be fit to empirical equations and written to different formats. Check the [Wiki](https://github.com/VlachosGroup/Thermochemistry/wiki) for more detailed explanations.

## Useful Topics
- [Outline to convert DFT data to empirical forms](https://github.com/VlachosGroup/PyMuTT/wiki/Converting-DFT-Generated-Data-to-Empirical-Models)
- [Explanation of enthalpy referencing](https://github.com/VlachosGroup/PyMuTT/wiki/Referencing)
- [Supported IO operations](https://github.com/VlachosGroup/PyMuTT/wiki/Input-and-Output)
- [Examples](https://github.com/VlachosGroup/PyMuTT/tree/master/examples)
- [How to contribute](https://github.com/VlachosGroup/PyMuTT/wiki/Contributing)

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
1. Install the dependencies
2. Download the repository to your local machine
3. Add to parent folder to PYTHONPATH
4. Run the tests by navigating to the [Thermochemistry/tests](https://github.com/VlachosGroup/Thermochemistry/tree/master/tests) folder and inputting the following command:
```
python -m unittest
```

The expected output is shown below. The number of tests will not necessarily be the same.
```
.........................
----------------------------------------------------------------------
Ran 25 tests in 0.020s

OK
```
5. Look at [examples using the code](https://github.com/VlachosGroup/PyMuTT/tree/master/examples)

## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/VlachosGroup/Thermochemistry/blob/master/LICENSE.md) file for details.
