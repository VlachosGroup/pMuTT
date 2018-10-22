.. _statmech:

Statistical Thermodynamic Models
********************************

StatMech Base Class
===================
The ``StatMech`` class stores the classes associated with each mode. Each mode has methods that specifies how to
calculate thermodynamic quantities, such as partition functions, heat capacities, internal energy, enthalpy, entropy,
Helmholtz energy and Gibbs energy.

 .. autoclass:: pMuTT.models.statmech.StatMech
   :members:

The ``StatMech`` object can be initialized in two ways:

1. by each passing each mode as objects or
2. passing each mode as a class and the associated parameters. 

The examples below give an equivalent result:

Example of initialization using objects
---------------------------------------

.. code:: python

   import numpy as np
   from ase.build import molecule
   from pMuTT.models.statmech import StatMech, trans, vib, rot, elec

   atoms = molecule('H2O')
   H2O_trans = trans.IdealTrans(n_degrees=3, atoms=atoms)
   H2O_vib = vib.HarmonicVib(vib_wavenumbers=[3825.434, 3710.2642, 1582.432])
   H2O_rot = rot.RigidRotor(symmetrynumber=2, atoms=atoms)
   H2O_elec = elec.IdealElec(potentialenergy=-14.2209, spin=0)
   H2O_statmech = StatMech(trans_model=H2O_trans,
                           vib_model=H2O_vib,
                           rot_model=H2O_rot,
                           elec_model=H2O_elec)
   
Example of initialization using classes and parameters
------------------------------------------------------

.. code:: python

   import numpy as np
   from ase.build import molecule
   from pMuTT.models.statmech import StatMech, trans, vib, rot, elec
   
   
   H2O_statmech = StatMech(trans_model=trans.IdealTrans,
                           n_degrees=3,
                           vib_model=vib.HarmonicVib,
                           vib_wavenumbers=[3825.434, 3710.2642, 1582.432],
                           rot_model=rot.RigidRotor,
                           symmetrynumber=2,
                           atoms=molecule('H2O'),
                           elec_model=elec.IdealElec,
                           potentialenergy=-14.2209,
                           spin=0)

Presets
=======
If you are using a common model (e.g. the ideal gas model), then you can get the default parameters from 
the dictionary, ``pMuTT.models.statmech.presets``. The same H2O StatMech object can be specified without
the need to pass all the types of modes:

.. code:: python

   from ase.build import molecule
   from pMuTT.models.statmech import StatMech, presets
   
   idealgas_defaults = presets['idealgas']
   H2O_new = StatMech(vib_wavenumbers=[3825.434, 3710.2642, 1582.432],
                      potentialenergy=-14.2209,
                      atoms=molecule('H2O'),
                      spin=0,
                      symmetrynumber=2,
                      **idealgas_defaults)
                      

Currently supported presets are described below. The first table shows the attributes already specified, 
and the second table shows the attributes that are still required, and the third table shows the required 
and optional parameters to calculate a thermodynamic property (where the value in parentheses is the default
value).

Ideal Gas (idealgas)
--------------------

+------------------+-----------------------------------------+
| Set Attributes   | Default Value                           |
+==================+=========================================+
| trans_model      | pMuTT.models.statmech.trans.IdealTrans  |
+------------------+-----------------------------------------+
| n_degrees        | 3                                       |
+------------------+-----------------------------------------+
| vib_model        | pMuTT.models.statmech.vib.HarmonicVib   |
+------------------+-----------------------------------------+
| elec_model       | pMuTT.models.statmech.elec.IdealElec    |
+------------------+-----------------------------------------+
| rot_model        | pMuTT.models.statmech.rot.RigidRotor    |
+------------------+-----------------------------------------+

+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Required Attributes | Description                                                                                                                                                |
+=====================+============================================================================================================================================================+
| molecular_weight    | (float) Molecular weight in g/mol                                                                                                                          |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| vib_wavenumbers     | (list of float) Vibrational wavenumbers in 1/cm                                                                                                            |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| potentialenergy     | (float) Electronic potential energy in eV                                                                                                                  |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| spin                | (float) Electron spin. 0 if all electrons are paired (e.g. N2), 0.5 if the specie is a radical (e.g. CH3.), 1 if the specie exists as a triplet (e.g. O2). |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| geometry            | (str) Geometry of molecule                                                                                                                                 |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| rot_temperatures    | (list of float) Rotational temperatures in K                                                                                                               |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| symmetrynumber      | (int) Symmetry number                                                                                                                                      |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| atoms               | (ase.Atoms object) Optional. If this parameter is specified, ``molecular_weight``, ``geometry``, and ``rot_temperatures`` do not have to be specified.     |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+

+-------------------------+---------------------+-----------------------+
| Thermodynamic Quantity  | Expected Parameters | Optional Parameters   |
+=========================+=====================+=======================+
| :math:`q`               | T                   | ignore_q_elec (False) |
|                         |                     | P (1.01325 bar)       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {C_V} {R}` | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {C_P} {R}` | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {U} {RT}`  | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {H} {RT}`  | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {S} {R}`   | T                   | P (1.01325 bar)       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {A} {RT}`  | T                   | P (1.01325 bar)       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {G} {RT}`  | T                   | P (1.01325 bar)       |
+-------------------------+---------------------+-----------------------+

Harmonic Approximation (harmonic)
---------------------------------

Typically used to model adsorbates.

+-------------+-----------------------------------------+
| Parameter   | Default Value                           |
+=============+=========================================+
| vib_model   | pMuTT.models.statmech.vib.HarmonicVib   |
+-------------+-----------------------------------------+
| elec_model  | pMuTT.models.statmech.elec.IdealElec    |
+-------------+-----------------------------------------+

+---------------------+--------------------------------------------------------------------------------------------------------------------------------------------------+
| Required Parameters | Description                                                                                                                                      |
+=====================+==================================================================================================================================================+
| vib_wavenumbers     | (list of float) Vibrational wavenumbers in 1/cm                                                                                                  |
+---------------------+--------------------------------------------------------------------------------------------------------------------------------------------------+
| potentialenergy     | (float) Electronic potential energy in eV                                                                                                        |
+---------------------+--------------------------------------------------------------------------------------------------------------------------------------------------+
| spin                | (float) Electron spin                                                                                                                            |
+---------------------+--------------------------------------------------------------------------------------------------------------------------------------------------+

+-------------------------+---------------------+-----------------------+
| Thermodynamic Quantity  | Expected Parameters | Optional Parameters   |
+=========================+=====================+=======================+
| :math:`q`               | T                   | ignore_q_elec (False) |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {C_V} {R}` | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {C_P} {R}` | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {U} {RT}`  | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {H} {RT}`  | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {S} {R}`   | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {A} {RT}`  | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {G} {RT}`  | T                   |                       |
+-------------------------+---------------------+-----------------------+


Electronic (electronic)
-----------------------

+-------------+-----------------------------------------+
| Parameter   | Default Value                           |
+=============+=========================================+
| elec_model  | pMuTT.models.statmech.elec.IdealElec    |
+-------------+-----------------------------------------+

+---------------------+--------------------------------------------------------------------------------------------------------------------------------------------------+
| Required Parameters | Description                                                                                                                                      |
+=====================+==================================================================================================================================================+
| potentialenergy     | (float) Electronic potential energy in eV                                                                                                        |
+---------------------+--------------------------------------------------------------------------------------------------------------------------------------------------+
| spin                | (float) Electron spin                                                                                                                            |
+---------------------+--------------------------------------------------------------------------------------------------------------------------------------------------+

+-------------------------+---------------------+-----------------------+
| Thermodynamic Quantity  | Expected Parameters | Optional Parameters   |
+=========================+=====================+=======================+
| :math:`q`               | T                   | ignore_q_elec (False) |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {C_V} {R}` |                     |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {C_P} {R}` |                     |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {U} {RT}`  | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {H} {RT}`  | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {S} {R}`   |                     |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {A} {RT}`  | T                   |                       |
+-------------------------+---------------------+-----------------------+
| :math:`\frac {G} {RT}`  | T                   |                       |
+-------------------------+---------------------+-----------------------+


The ``presets`` dictionary is flexible where you can create a new entry if you will use a model often.

.. autodata:: pMuTT.models.statmech.presets


Translational Models
====================

IdealTrans
----------
.. autoclass:: pMuTT.models.statmech.trans.IdealTrans
   :members:

Vibrational Models
==================

HarmonicVib
-----------
.. autoclass:: pMuTT.models.statmech.vib.HarmonicVib
   :members:

QRRHOVib
--------
.. autoclass:: pMuTT.models.statmech.vib.QRRHOVib
   :members:

EinsteinVib
-----------
.. autoclass:: pMuTT.models.statmech.vib.EinsteinVib
   :members:
   
Rotational Models
=================

RigidRotor
----------
.. autoclass:: pMuTT.models.statmech.rot.RigidRotor
   :members:

Electronic Models
=================

IdealElec
---------
.. autoclass:: pMuTT.models.statmech.elec.IdealElec
   :members:

Nuclear Models
==============
Typically these are unimportant for chemical reactions, but the module is present in case
nuclear modes become important in the future.

.. automodule:: pMuTT.models.statmech.nucl
   :members:

Creating New StatMech Models
============================
If you would like to create your own class to be stored by ``StatMech``, it should have the methods: 
``get_q``, ``get_CvoR``, ``get_CpoR``, ``get_UoRT``, ``get_HoRT``, ``get_SoR``, ``get_AoRT``, ``get_GoRT``,
``from_dict``, and ``to_dict``.
