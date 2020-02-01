.. _reactions:

Reactions
*********

.. currentmodule:: pmutt

Reaction Objects
================

These classes are used to model a single chemical reaction.
:class:`~pmutt.reaction.ChemkinReaction` and
:class:`~pmutt.omkm.reaction.SurfaceReaction` have additional formatting options
for the appropriate modeling software.

.. autosummary::
   :toctree: reaction
   :nosignatures:

   reaction.Reaction
   reaction.ChemkinReaction
   omkm.reaction.SurfaceReaction

Reactions Objects
=================

These classes model multiple reactions. Multiple reactions can be used together
for phase diagrams, reaction coordinate diagrams, and other purposes.

.. autosummary::
   :toctree: reactions
   :nosignatures:

   reaction.Reactions
   reaction.ChemkinReaction
   omkm.reaction.SurfaceReaction

BEP
===

Bronsted Evans Polyani relationships can be added to 
:class:`~pmutt.reaction.Reaction` objects so the activation energy can be 
estimated using the change in enthalpy.

.. autosummary::
   :toctree: bep
   :nosignatures:

   reaction.bep.BEP
   omkm.reaction.BEP