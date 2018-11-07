.. _release_notes:

Release Notes
*************

Version 1.1.0
-------------
Nov. 7, 2018
- Fixed bugs in ``Shomate`` object for ``get_HoRT`` and ``get_SoR`` where one temperature would return a 1x1 vector
  instead of a float
- Fixed bug in ``Zacros`` object where it expected vibrational energies instead of wavenumbers.

Version 1.1.0
-------------
Oct. 26, 2018

- Updated Reaction class to parse strings
- New Shomate class
- New eos module (Equations of state)
- New PhaseDiagram class
- New EinsteinVib class
- New functions to read species and reactions from Chemkin surf.inp and gas.inp files
