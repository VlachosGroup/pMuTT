.. _kinetic_models:

Kinetic Models
**************

.. currentmodule:: pmutt

Chemkin
=======

Classes and functions related to writing Chemkin files. Please note these
modules were written specifically for the Vlachos group's in-house version.
Features may or may not be supported in the commercial version.

--------------------------------------------------------------------------------

Classes
-------

.. autosummary::
   :toctree: chemkin
   :nosignatures:

   chemkin.CatSite
   reaction.ChemkinReaction

Input and Output
----------------

.. autosummary::
   :nosignatures:

   io.chemkin.read_reactions
   io.chemkin.write_EA
   io.chemkin.write_gas
   io.chemkin.write_surf
   io.chemkin.write_T_flow
   io.chemkin.write_tube_mole

Examples
--------

- :ref:`chemkin_io_example`

--------------------------------------------------------------------------------

Cantera
=======

Classes and functions related to `Cantera`_. This functionality is still in
its early stages.

Phases
------

.. autosummary::
   :toctree: cantera_phase
   :nosignatures:

   cantera.phase.Phase
   cantera.phase.IdealGas
   cantera.phase.StoichSolid

Units
-----

.. autosummary::
   :toctree: cantera_units
   :nosignatures:

   cantera.units.Units

Input and Output
----------------

.. autosummary::
   :toctree: cantera_io
   :nosignatures:

   io.cantera.obj_to_cti

--------------------------------------------------------------------------------

OpenMKM
=======

Classes and functions related to OpenMKM. This functionality is still in its
early stages.

Phases
------

.. autosummary::
   :toctree: omkm_phases
   :nosignatures:

   omkm.phase.IdealGas
   omkm.phase.StoichSolid
   omkm.phase.InteractingInterface

Reactions
---------

.. autosummary::
   :toctree: omkm_reactions
   :nosignatures:

   omkm.reaction.SurfaceReaction
   omkm.reaction.BEP

Units
-----

.. autosummary::
   :toctree: omkm_units
   :nosignatures:

   omkm.units.Units

Input and Output
----------------

.. autosummary::
   :toctree: omkm
   :nosignatures:

   io.omkm.write_cti
   io.omkm.write_yaml
   io.omkm.get_species_phases
   io.omkm.get_reactions_phases
   io.omkm.get_interactions_phases
   io.omkm.organize_phases

Examples
--------
- :ref:`openmkm_io_example`

--------------------------------------------------------------------------------

Zacros
======

Classes related to `Zacros`_.

Classes
-------

.. autosummary::
   :toctree: zacros
   :nosignatures:

   empirical.zacros.Zacros

.. _`Cantera`: https://cantera.org/
.. _`Zacros`: http://zacros.org/