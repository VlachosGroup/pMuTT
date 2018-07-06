# Thermochemistry
This repository contains a Python library for Thermochemistry developed by the Vlachos Research Group at the University of Delaware.

## Planned features
1. Read electronic energies and frequencies from a variety of formats, including:
   - VASP
   - Gaussian
   - CSV files of particular format
2. Calculate thermodynamic properties, such as heat capacity (C<sub>P</sub>), enthalpy (H), entropy (S), Gibbs energy (G), using DFT-generated data.
3. Read and write to a variety of formats, including:
   - NASA polynomials for Chemkin
   - Shomate polynomials

## Developers
- Gerhard Wittreich, P.E. (wittregr@udel.edu)
- Jonathan Lym (jlym@udel.edu)

## Dependencies
- [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/): Used for I/O operations and to calculate thermodynamic properties
- [Numpy](http://www.numpy.org/): Used for vector and matrix operations
   
## Getting Started
1. Download the repository to your local machine
2. Add to PYTHONPATH