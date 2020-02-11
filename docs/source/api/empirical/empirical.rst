.. _empirical:

Empirical Models
****************

Empirical models predict thermodynamic properties of species using simple
polynomials. These can be evaluated more quickly than statistical mechanical
equations and tend to be used by kinetic modeling software.

When fitting an empirical relationship, be mindful that the resulting polynomial
will have the same reference state as the original data. Databases like the
`NIST Chemistry WebBook`_ and `Burcat`_ tend to use the standard reference
(i.e. pure elements in their natural state at 1 atm have a standard enthalpy of
formation of 0). However, data from DFT software (and DFT-based references) may
use different references. See :ref:`referencing` for more information.

.. currentmodule:: pmutt.empirical

Parent Class
------------

.. autosummary::
   :toctree: empirical_base
   :nosignatures:

   EmpiricalBase

--------------------------------------------------------------------------------

Nasa
----

.. autosummary::
   :toctree: nasa
   :nosignatures:

   nasa.Nasa
   nasa.SingleNasa9
   nasa.Nasa9

--------------------------------------------------------------------------------

Shomate
-------

.. autosummary::
   :toctree: shomate
   :nosignatures:

   shomate.Shomate

--------------------------------------------------------------------------------

Referencing
-----------

Referencing is used to bridge the gap between *ab-initio* code reference states
and the standard state reference. See :ref:`referencing` for more information.

.. autosummary::
   :toctree: references
   :nosignatures:

   references.Reference
   references.References

--------------------------------------------------------------------------------

Misc.
-----

Miscellaneous models associated with the empirical models are located here.

.. autosummary::
   :toctree: misc
   :nosignatures:

   GasPressureAdj


.. _`NIST Chemistry WebBook`: https://webbook.nist.gov/chemistry/
.. _`Burcat`: http://garfield.chem.elte.hu/Burcat/burcat.html