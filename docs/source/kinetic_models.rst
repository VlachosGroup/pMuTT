.. _kinetic_models:

Kinetic Models
**************

Chemkin
=======
Classes and functions related to writing Chemkin files. Please note these
modules were written specifically for the Vlachos group's in-house version.
Features may or may not be supported in the commercial version.

Classes
-------

.. autoclass:: pmutt.chemkin.CatSite
   :members:

.. autoclass:: pmutt.reaction.ChemkinReaction
   :members:

Input and Output
----------------

.. automodule:: pmutt.io.chemkin
   :members:

Examples
--------
- :ref:`chemkin_io_example`

Cantera
=======
Classes and functions related to `Cantera`_. This functionality is still in
its early stages.

Classes
-------

.. automodule:: pmutt.cantera.phase
   :members:

.. automodule:: pmutt.cantera.reaction
   :members:

.. automodule:: pmutt.cantera.units
   :members:

Input and Output
----------------

.. automodule:: pmutt.io.cantera
   :members:

OpenMKM
=======

Classes and functions related to OpenMKM. This functionality is still in its
early stages.

Classes
-------

.. automodule:: pmutt.omkm.phase
   :members:

.. automodule:: pmutt.omkm.reaction
   :members:

.. automodule:: pmutt.omkm.units
   :members:

Input and Output
----------------

.. automodule:: pmutt.io.omkm
   :members:

Examples
--------
- :ref:`openmkm_io_example`

Zacros
======

Classes related to `Zacros`_.

Classes
-------

.. autoclass:: pmutt.empirical.zacros.Zacros
   :members:

.. _`Cantera`: https://cantera.org/
.. _`Zacros`: http://zacros.org/