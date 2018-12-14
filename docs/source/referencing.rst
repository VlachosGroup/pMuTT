.. _referencing:

Referencing
===========

Enthalpies calculated using VASP (and some other computational methods) 
have different references than standard references (i.e. the enthalpy of 
formation of pure substances, like O\ :sub:`2`\  or Pt, is not necessarily 
zero). This difference makes it difficult to ensure thermodynamic consistency 
for our mechanisms since we may be mixing experimental gas thermodynamics with 
computational surface thermodynamics. In order to make the references 
consistent, we find a correction factor for each element by solving the 
equation:

:math:`\underline {\underline {H}}^{expt}_{[M \times N]} = \underline 
{\underline {H}}^{DFT}_{[M \times N]} + \underline {\underline {x}}_
{[M \times N]} \underline{\theta}_{[N]}`

where M is the number of reference species, N is the number of elements, 
H\ :sup:`expt`\  is the experimental standard enthalpies, 
H\ :sup:`DFT`\  is the standard enthalpies calculated using DFT, 
x is a matrix that describes the composition of the references (each row 
represents a specie, each column represents an element), and θ is the 
correction for each element.

The equation can be solved using a Least Squares approach. 
The correction factor can then be added to subsequent species calculated 
through DFT to ensure consistent references. Referencing is handled by the 
:class:`~pMuTT.empirical.references.References`.
