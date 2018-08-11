Referencing
-----------

Enthalpies calculated using VASP (and some other computational methods)
have different references than standard references (i.e. the enthalpy of
formation of pure substances, like O2 or Pt, is not necessarily zero).
This difference makes it difficult to ensure thermodynamic consistency
for our mechanisms since we may be mixing experimental gas
thermodynamics with computational surface thermodynamics. In order to
make the references consistent, we find a correction factor for each
element by solving the equation:

.. figure:: https://github.com/VlachosGroup/PyMuTT/blob/master/docs/README_Eq1.gif
   :alt: Eq1

   Eq1

where M is the number of reference species, N is the number of elements,
Hexpt is the experimental standard enthalpies, HDFT is the standard
enthalpies calculated using DFT, x is a matrix that describes the
composition of the references (each row represents a specie, each column
represents an element), and θ is the correction for each element.

The equation can be solved using a Least Squares approach. The
correction factor can then be added to subsequent species calculated
through DFT to ensure consistent references. Referencing is handled by
the ```PyMuTT.models.empirical.references.References```_ class.

.. _``PyMuTT.models.empirical.references.References``: https://github.com/VlachosGroup/PyMuTT/blob/master/models/empirical/references.py