.. _release_notes:

Release Notes
*************

Git Master Branch
-----------------
`Github Master Branch`_

- Fixed bug related to :class:`~pMuTT.empirical.references.References` and 
  :class:`~pMuTT.empirical.references.Reference` objects not JSON-write 
  compatible.
- Fixed bug related to referencing in :class:`~pMuTT.empirical.shomate.Shomate`
  class

Version 1.2.1
-------------
Dec. 17, 2018

- Added ``vib_outcar`` special rule for :func:`~pMuTT.io_.excel.read_excel` and
  :func:`~pMuTT.io_.vasp.set_vib_wavenumbers_from_outcar` to get vibrational 
  frequencies directly from OUTCAR file.
- Added ``get_X`` methods to :class:`~pMuTT.empirical.nasa.Nasa`, 
  :class:`~pMuTT.empirical.shomate.Shomate`, :class:`~pMuTT.statmech.StatMech` 
  and :class:`~pMuTT.reaction.Reaction` to directly calculate thermodynamic 
  properties (such as H, S, F, G) with the appropriate units
- Changed symbol for Hemlholtz energy from A to F

Version 1.2.0
-------------
Dec. 12, 2018

- Restructured code to exclude ``model`` module

Version 1.1.3
-------------
Dec. 11, 2018

- Added :class:`~pMuTT.reaction.bep.BEP` class
- Restructured :class:`~pMuTT.reaction.Reaction` class so reaction states (i.e. reactants, products, transition states) can be calculated separately
- Updated :class:`~pMuTT.empirical.references.References` class to be able reference any attribute
- Added ``placeholder`` entry to :data:`~pMuTT.statmech.presets` dictionary to represent an empty species
- Added correction factor to calculate partition coefficient, q, in :class:`~pMuTT.statmech.elec.IdealElec` class

Version 1.1.2
-------------
Nov. 27, 2018

- Fixed bugs in :class:`~pMuTT.reaction.Reaction` class for calculating pre-exponential factors
- Added methods in :class:`~pMuTT.reaction.Reaction` class to calculate rate constants and activation energy (currently, this just calculates the difference in enthalpy between the reactant/product and the transition state)
- Quality of life improvements such as allowing :class:`~pMuTT.reaction.Reaction` class inputs to be a single pMuTT object instead of expecting a list

Version 1.1.1
-------------
Nov. 7, 2018

- Fixed bugs in :class:`~pMuTT.empirical.shomate.Shomate` class for ``get_HoRT`` and ``get_SoR`` where one temperature would return a 1x1 vector
  instead of a float
- Fixed bug in :class:`~pMuTT.empirical.zacros.Zacros` class where it expected vibrational energies instead of wavenumbers.

Version 1.1.0
-------------
Oct. 26, 2018

- Updated :class:`~pMuTT.reaction.Reaction` class to parse strings
- New :class:`~pMuTT.empirical.shomate.Shomate` class
- New equation of state classes: :class:`~pMuTT.eos.IdealGasEOS`, :class:`~pMuTT.eos.vanDerWaalsEOS`
- New :class:`~pMuTT.reaction.phasediagram.PhaseDiagram` class
- New :class:`~pMuTT.statmech.vib.EinsteinVib` class
- New :func:`~pMuTT.io_.chemkin.read_reactions` function to read species and reactions from Chemkin surf.inp and gas.inp files

.. _`Github Master Branch`: https://github.com/VlachosGroup/pMuTT/commits/master