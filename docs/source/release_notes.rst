.. _release_notes:

Release Notes
*************

Development Branch
------------------
`Development Branch`_

Version 1.4.17
--------------

May 25, 2025

- Added feature where you can specify 'pymatgen' as the symmetry number for
  a gas molecule and Pymatgen will determine the symmetry number
  
- Added a check that forces temperatures submitted for thermodynamic
  calculations to be astype=float to avoid calculation issues encountered
  when temperatures are entered as astype=int.

Version 1.4.15
--------------

April 21, 2024

- Fixed reaction lists in BEP section of thermo.yaml openMKM file
  so reactions are all listed individually and no longer grouped.

Version 1.4.14
--------------

February 8, 2024

- Fix equilibrium unittest warning messages

Version 1.4.13
--------------

February 4, 2024

- Fix equilibrium class-Upper bound solver violation

Version 1.4.11
--------------

February 3, 2024

- Fix unittest for equilibrium class-Missing thermdat file
- Fix equilibrium class-Lower bound solver violation

Version 1.4.10
--------------

February 2, 2024

- Improved NASA9 polynomial fit function
- Overall clean-up of NASA9 deprecated code

Version 1.4.9
-------------

October 25, 2023

- Fixed deprecated numpy function impacting nasa.py and shomate.py
  Conversion of an array with ndim > 0 to a scalar

Version 1.4.7
-------------

August 23, 2023

- Fixed improper reaction and lateral interaction ID's in write_CTI

Version 1.4.6
-------------

August 15, 2023

- Adding missing dependency for openpyxl

Version 1.4.5
-------------

June 23, 2023

- Updated .iteritems dpecricated function in Pandas
- Cleanup unused imports

Version 1.4.4
-------------

May 1, 2023

- CRITICAL UPDATE: Fixed Nasa polynomial fit issue for enthalpy and entropy when T_ref > T_mid

Version 1.4.3
-------------

Mar. 19, 2023

- Updated PIP installation dependency requirements

Version 1.4.2
-------------

Mar. 18, 2023

- Added new :class:`~pmutt.equilibrium.Equilibrium` class that computes the equilibrium
  concentration on a network of species by minimizing the network
  Gibbs free energy while maintaining the total atom balance of
  the initial starting concentration of species
- Major GitHub documentation update including the new functionality
  from v1.3.2 and v1.4.0, corrections noted in the issues database,
  fixes to broken links, addition of a new "Helper Functions. section,
  and misc upgrades.

Version 1.3.2
-------------

Jan. 26, 2023

- Added ``S_elements`` parameter to entropy and Gibbs free energy
  methods allowing you to compute an entropy and Gibbs free energy
  of formation for a single species by including the entropy of the
  elements in the species.
- Misc bug fixes

Version 1.2.21
--------------

Jul. 3, 2020

- Miscellaneous bug fixes for :func:`~pmutt.io.omkm.write_cti`
- Added preliminary support for OpenMKM YAML files.
- Fixed bug where test OUTCAR file was not present.
- Added ``ads_act_method`` to :func:`~pmutt.io.chemkin.write_surf` and
  :func:`~pmutt.io.chemkin.write_EA` to allow users to specify a different
  method to calculate activation energies for surface versus adsorption
  reactions.

Version 1.2.20
--------------

May 12, 2020

- Added sensitivity analysis options in :func:`~pmutt.io.omkm.write_cti`
- Added ability to arbitrarily specify pre-exponential constant and activation
  energy in :class:`~pmutt.omkm.reaction.SurfaceReaction`
- Reduced complexity of assigning reaction IDs to
  :class:`~pmutt.omkm.reaction.SurfaceReaction`. Before, the ID would be
  renamed if it had a BEP associated with it. Changed the behavior to just
  assign numerical values.

Version 1.2.19
--------------

Apr. 8, 2020

- Fixed bug where :class:`~pmutt.statmech.lsr.LSR` could not be imported in
  :func:`~pmutt.io.excel.read_excel`.
- Updated :mod:`~pmutt.io` sections to incorporate Pathlib library
- Added helper functions in :mod:`~pmutt.io.omkm` to organize phases
- Fixed bug where slopes and y intercepts were switched for
  :class:`~pmutt.mixture.cov.PiecewiseCovEffect` when writing CTI files.
- Fixed bug in :func:`~pmutt.io.chemkin.write_surf` where ``n_sites`` could
  be written as a float
- Updated OpenMKM IO example


Version 1.2.18
--------------

Jan. 31, 2020

- Hotfix to correct broken links in documentation.

Version 1.2.17
--------------

Jan. 31, 2020

- Added more descriptive warning messages when incorrect temperature values are
  passed to :class:`~pmutt.empirical.nasa.Nasa`,
  :class:`~pmutt.empirical.nasa.Nasa9`, and
  :class:`~pmutt.empirical.shomate.Shomate`.
- Fixed bug where the conversion factor for Hartrees was incorrect.
- Added extra parameters for OpenMKM IO.
- Added helper functions for OpenMKM IO to assign phases easily.
- Added a helper method in :class:`~pmutt.omkm.cantera.IdealGas`
  and :class:`~pmutt.omkm.cantera.StoichSolid` to only assign a reaction to the
  phase if all the species belong to that phase.
- Fixed outdated code in Chemkin example and OpenMKM example.
- Reorganized documentation to use stubs. Shorter pages should hopefully
  make the documentation easier to navigate.

Version 1.2.16
--------------
Dec. 9, 2019

- Hotfix to correct a typo for PyYAML version required.


Version 1.2.15
--------------
Dec. 5, 2019

- Added :func:`~pmutt.io.omkm.write_yaml` to write YAML files for OpenMKM.
- Added warning for :func:`~pmutt.io.excel.read_excel` if the header is blank
  but the cells are occupied.
- Fixed bug in :func:`~pmutt.io.excel.read_excel` where ``model`` was not
  correctly initialized with :func:`~pmutt.statmech.StatMech`.
- Added the generic method, :func:`~pmutt.io.excel.set_dict_value` to specify
  dictionaries in :func:`~pmutt.io.excel.read_excel`
- Removed redundant statements involving returning dictionaries in functions to
  process Excel data.
- Fixed warning raised whenever CpoR = 0 when fitting empirical polynomials.

Version 1.2.14
--------------
Oct. 25, 2019

- Added functionality to write files (such as 
  :func:`~pmutt.io.thermdat.write_thermdat`) can return a string containing
  the file if ``filename`` is not specified.
- Bug fix where ``from_model`` for :class:`~pmutt.empirical.nasa.Nasa` and
  :class:`~pmutt.empirical.shomate.Shomate` returned errors. The fix was
  related to incorrect datatyping for ``misc_models``.
- Improved :class:`~pmutt.empirical.shomate.Shomate` to allow users to specify
  the units for the polynomial coefficients.
- Energies from Gaussian input functions (:mod:`~pmutt.io.gaussian`)
  was originally in Hartrees. Changed to allow users to specify what unit they
  desire (default in eV).
- Added functionality to write BEP relationships to OpenMKM CTI files.
- Restructured OpenMKM CTI writer to be more robust when specifying custom IDs
- Added functionality to remove leading and trailing spaces when reading from
  Excel sheets since users found this error hard to pick up.

Contributors
^^^^^^^^^^^^
- Qiang Li (lqcata_)

Version 1.2.13
--------------
Oct. 2, 2019

- Fixed bug where small non-zero rotational inertia modes were chosen
  preferentially over larger contributing modes.
- Fixed bug where presets had to be specified before statistical mechanical
  arguments. Now, the preset will not overwrite any previously set values.
- Updated :func:`~pmutt.io.thermdat.read_thermdat` to allow the user to return
  the :class:`~pmutt.empirical.nasa.Nasa` objects as a list, tuple, or
  dictionary.
- Updated :func:`~pmutt.io.thermdat.write_thermdat` to accept a list or a
  dictionary of :class:`~pmutt.empirical.nasa.Nasa` objects
- Implemented `from_model` method in :class:`~pmutt.empirical.nasa.Nasa` and
  :class:`~pmutt.empirical.shomate.Shomate` classes so empirical objects can be
  created from :class:`~pmutt.statmech.StatMech` objects as well as other
  empirical objects. The ``from_statmech`` method is deprecated.
- Added more descriptive warnings and errors.
- Created :class:`~pmutt.empirical.GasPressureAdj` so entropy and Gibbs energy
  of gas-phase empirical objects (like :class:`~pmutt.empirical.shomate.Shomate`
  and :class:`~pmutt.empirical.nasa.Nasa`) are dependent on pressure. This
  object is assigned automatically to ``misc_models`` if ``phase`` is 'g' or
  'gas' and the ``add_gas_P_adj`` can be set to False if users do not wish to
  assign this object automatically.
- Thermodynamic quantities of individual species can also be calculated on a
  per mass basis (i.e. users can calculate quantities in J/g, cal/kg, etc.).
  The object must contain a dictionary of its composition in ``elements`` for
  this functionality.
- Fixed broken hyperlinks.

Contributors
^^^^^^^^^^^^
- Geun Ho Gu (googhgoo_)

Version 1.2.12
--------------
Aug. 22, 2019

- Refactored :class:`~pmutt.io.thermdat.write_thermdat` so that it is simpler
  to understand
- Implemented :class:`~pmutt.empirical.nasa.Nasa9` and 
  :class:`~pmutt.empirical.nasa.SingleNasa9` polynomials
- Added preliminary CTI file writer for Cantera and OpenMKM
- Added Binder notebooks to Examples page so users can try pMuTT before
  installing
- Fixed bug where :class:`~pmutt.statmech.StatMech` was not passed when
  modes were specified indivudally in spreadsheets.

Contributors
^^^^^^^^^^^^
Xenhua Zhang (xenhua_)

Version 1.2.11
--------------
Jun. 18, 2019

- Added xlrd dependency so spreadsheets can be read using pandas
- Updated documentation page with NAM 2019 instructions.

Version 1.2.10
--------------
Jun. 13, 2019

- Another hotfix to fix a bug where the version was not incremented correctly

Version 1.2.9
-------------
Jun. 13, 2019

- Hotfix where pypi created the folder in the old case (pMuTT) instead of
  lower case (pmutt)

Version 1.2.8
-------------
Jun. 13, 2019

- Importing from pMuTT is now all in lowercase. (i.e. ``import pmutt`` instead
  of ``import pMuTT``)

Version 1.2.7
-------------
Jun. 11, 2019

- Added documentation page for more verbose installation instructions.
- Updated :class:`~pmutt.reaction.network` to use graph theory approach using
  states as nodes
- Bug fix for :class:`~pmutt.statmech.lsr.LSR` to handle inputs that are not
  pmutt model objects
- Added ability to create interactive plots with Pygal
- Updated :class:`~pmutt.statmech.elec.GroundStateElec` to read
  ``potentialenergy`` from inputted ``Atoms`` object.

Version 1.2.6
-------------
Apr. 26, 2019

- Moved ``references`` attribute from empirical classes to
  :class:`~pmutt.statmech.StatMech`
- Changed ``mix_models`` attribute to ``misc_models`` in  indicating any model
  object can be used
- Implemented :class:`~pmutt.statmech.vib.DebyeVib` and
  :class:`~pmutt.statmech.ConstantMode` classes
- Restructured :class:`~pmutt.reaction.bep.BEP` object to act as a transition
  state species in :class:`~pmutt.reaction.Reaction` objects
- Implemented :class:`~pmutt.empirical.lsr.LSR` object
- Added option to calculate pre-exponential factor using ratio of partition
  functions or entropy of activation
- Added option to use electronic energy as descriptor for
  :class:`~pmutt.reaction.bep.BEP` object
- Added some imperial unit functionality to ``pmutt.constants`` module
- Renamed ``from_`` parameter and ``to`` parameter in 
  :func:`pmutt.constants.convert_unit` to ``initial`` and ``final``
- Added ability to import individual translational, rotational, vibrational,
  electronic and nuclear modes to Excel
- Renamed ``pmutt.statmech.trans.IdealTrans`` to
  :class:`~pmutt.statmech.trans.FreeTrans`
- Renamed ``pmutt.statmech.elec.IdealElec`` to
  :class:`~pmutt.statmech.elec.GroundStateElec`
- Renamed ``pmutt.statmech.nucl.IdealNucl`` to
  :class:`~pmutt.statmech.nucl.EmptyNucl`

Version 1.2.5
-------------
Mar. 21, 2019

- Renamed ``pmutt.io_`` module to ``pmutt.io``
- Renamed ``pmutt.io_.jsonio`` module to ``pmutt.io.json``
- Added preliminary IO support for MongoDB in module: ``pmutt.io.db``
- Bug fixes for Chemkin IO behavior

Version 1.2.4
-------------
Mar. 11, 2019

- Hotfix to correct Chemkin IO behavior

Version 1.2.3
-------------
Feb. 25, 2019

- Added ``smiles`` attribute to :class:`~pmutt.statmech.StatMech` and 
  :class:`~pmutt.empirical.EmpiricalBase` classes
- Added functions to write Chemkin surf.inp, gas.inp, and EAs.inp files
- Added :class:`~pmutt.mixture.cov.CovEffect` class to model coverage effects
  and integrated it with :class:`~pmutt.statmech.StatMech` and 
  :class:`~pmutt.empirical.EmpiricalBase` classes
- Added ``include_ZPE`` parameter to ``get_EoRT``, ``get_E``, ``get_delta_EoRT``
  and ``get_delta_E`` for the :class:`~pmutt.statmech.StatMech` class and
  :class:`~pmutt.reaction.Reaction` class to add zero-point energy in
  calculations
- Renamed private methods ``_get_delta_quantity`` and ``_get_state_quantity`` to
  public methods ``get_delta_quantity`` and ``get_state_quantity`` in
  :class:`~pmutt.reaction.Reaction` class
- Added generic method ``get_quantity`` to :class:`~pmutt.statmech.StatMech`
  class so any method can be evaluated. It takes the parameters ``raise_error``
  and ``raise_warning`` so the user has the ability to ignore modes if they do
  not have the desired properties
- Added ``plot_coordinate_diagram`` method to the 
  :class:`~pmutt.reaction.Reactions` class to plot coordinate diagrams.
- Added ``get_EoRT`` and ``get_E`` methods to :class:`~pmutt.statmech.StatMech`
  class to calculate electronic contribution to thermodynamic properties
- Added ``get_EoRT_state`` and ``get_delta_EoRT`` methods to 
  :class:`~pmutt.reaction.Reaction` to calculate electronic contribution to
  reaction properties
- Added an optional parameter, ``activation``, to ``get_delta_X`` methods to 
  specify the difference between the reactants/products and the transition
  state. 
- Added ``pmutt.constants.symmetry_dict`` to allow easy look up of common
  symmetry numbers
- Fixed bug where specie-specific arguments were not passed correctly for
  :class:`~pmutt.reaction.Reaction` class

Version 1.2.2
-------------
Jan. 18, 2019

- Added option to extract imaginary frequencies from VASP's OUTCAR files
- Added support for imaginary frequencies for 
  :class:`~pmutt.statmech.vib.HarmonicVib` and 
  :class:`~pmutt.statmech.vib.QRRHOVib` classes
- Restructured :class:`~pmutt.statmech.vib.HarmonicVib` and 
  :class:`~pmutt.statmech.vib.QRRHOVib` classes to calculate vibrational 
  temperatures, scaled wavenumbers and scaled inertia when methods are called 
  (rather than at initialization) to prevent incorrect calculations due to 
  changes in the vibrational wavenumbers.
- Fixed unit test names
- Added ``get_species`` to :class:`~pmutt.reaction.Reaction` and 
  :class:`~pmutt.reaction.Reactions`
- Fixed bug related to :class:`~pmutt.empirical.references.References` and 
  :class:`~pmutt.empirical.references.Reference` objects not JSON-write 
  compatible.
- Fixed bug related to referencing in :class:`~pmutt.empirical.shomate.Shomate`
  class

Version 1.2.1
-------------
Dec. 17, 2018

- Added ``vib_outcar`` special rule for :func:`~pmutt.io.excel.read_excel` and
  :func:`~pmutt.io.vasp.set_vib_wavenumbers_from_outcar` to get vibrational 
  frequencies directly from VASP's OUTCAR file.
- Added ``get_X`` methods to :class:`~pmutt.empirical.nasa.Nasa`, 
  :class:`~pmutt.empirical.shomate.Shomate`, :class:`~pmutt.statmech.StatMech` 
  and :class:`~pmutt.reaction.Reaction` to directly calculate thermodynamic 
  properties (such as H, S, F, G) with the appropriate units
- Changed symbol for Hemlholtz energy from A to F

Contributors
^^^^^^^^^^^^
- Himaghna Bhattacharjee (himaghna_)

Version 1.2.0
-------------
Dec. 12, 2018

- Restructured code to exclude ``model`` module

Version 1.1.3
-------------
Dec. 11, 2018

- Added :class:`~pmutt.reaction.bep.BEP` class
- Restructured :class:`~pmutt.reaction.Reaction` class so reaction states (i.e.
  reactants, products, transition states) can be calculated separately
- Updated :class:`~pmutt.empirical.references.References` class to be able
  reference any attribute
- Added ``placeholder`` entry to :data:`~pmutt.statmech.presets` dictionary to
  represent an empty species
- Added correction factor to calculate partition coefficient, q, in
  :class:`~pmutt.statmech.elec.IdealElec` class

Version 1.1.2
-------------
Nov. 27, 2018

- Fixed bugs in :class:`~pmutt.reaction.Reaction` class for calculating
  pre-exponential factors
- Added methods in :class:`~pmutt.reaction.Reaction` class to calculate rate
  constants and activation energy (currently, this just calculates the 
  difference in enthalpy between the reactant/product and the transition state)
- Quality of life improvements such as allowing
  :class:`~pmutt.reaction.Reaction` class inputs to be a single pmutt object
  instead of expecting a list

Version 1.1.1
-------------
Nov. 7, 2018

- Fixed bugs in :class:`~pmutt.empirical.shomate.Shomate` class for ``get_HoRT``
  and ``get_SoR`` where one temperature would return a 1x1 vector instead of a
  float
- Fixed bug in :class:`~pmutt.empirical.zacros.Zacros` class where it expected
  vibrational energies instead of wavenumbers.

Version 1.1.0
-------------
Oct. 26, 2018

- Updated :class:`~pmutt.reaction.Reaction` class to parse strings
- New :class:`~pmutt.empirical.shomate.Shomate` class
- New equation of state classes: :class:`~pmutt.eos.IdealGasEOS`,
  :class:`~pmutt.eos.vanDerWaalsEOS`
- New :class:`~pmutt.reaction.phasediagram.PhaseDiagram` class
- New :class:`~pmutt.statmech.vib.EinsteinVib` class
- New :func:`~pmutt.io.chemkin.read_reactions` function to read species and
  reactions from Chemkin surf.inp and gas.inp files

.. _`Development Branch`: https://github.com/VlachosGroup/pmutt/commits/development
.. _himaghna: https://github.com/himaghna
.. _xenhua: https://github.com/xenhua
.. _googhgoo: https://github.com/googhgoo
.. _lqcata: https://github.com/lqcata