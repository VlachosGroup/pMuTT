.. _io:

Input and Output
****************

Input and output to different forms is an active area of development for pMuTT.

.. currentmodule:: pmutt.io

Excel
=====

.. autosummary::
   :toctree: excel
   :nosignatures:

   excel.read_excel

Special Rules
-------------

Special rules can be defined in the :func:`~pmutt.io.excel.read_excel` 
function to process inputs differently. Currently supported special rules are 
listed below.

.. autosummary::
   :toctree: excel
   :nosignatures:

   excel.set_element
   excel.set_formula
   excel.set_atoms
   excel.set_statmech_model
   excel.set_trans_model
   excel.set_vib_model
   excel.set_rot_model
   excel.set_elec_model
   excel.set_nucl_model
   excel.set_vib_wavenumbers
   excel.set_rot_temperatures
   excel.set_nasa_a_low
   excel.set_nasa_a_high
   excel.set_list_value
   excel.set_dict_value

--------------------------------------------------------------------------------

Thermdat
========

This is the output format used for Chemkin. A list of NASA objects can be 
written to a thermdat file.

.. autosummary::
   :toctree: thermdat
   :nosignatures:

   thermdat.read_thermdat
   thermdat.write_thermdat

--------------------------------------------------------------------------------

JSON
====

`JavaScript Object Notation (JSON)`_ is a format that is easily read and
written by both humans and machines. All pmutt objects with the methods 
``to_dict`` and ``from_dict`` are JSON compatible.

.. autosummary::
   :toctree: json
   :nosignatures:

   json.pmuttEncoder
   json.json_to_pmutt

Examples
--------

Saving pmutt objects can be done by using 
:func:`~pmutt.io.json.pmuttEncoder`.

.. code:: python

   import json
   from pmutt.io.json import pmuttEncoder
   
   with open(json_path, 'w') as f_ptr:
       json.dump(pmutt_obj, f_ptr, cls=pmuttEncoder, indent=True)
   
Loading pmutt objects can be done by using the object hook: 
:func:`~pmutt.io.json.json_to_pmutt`.

.. code:: python

   import json
   from pmutt.io.json import json_to_pmutt

   with open(json_path, 'r') as f_ptr:
       pmutt_obj = json.load(f_ptr, object_hook=json_to_pmutt)

Sample JSON File
----------------

JSON writes in a human-readable syntax. An example showing a H2 
:class:`~pmutt.empirical.shomate.Shomate` object in JSON format is shown below.
::

   {
    "class": "<class 'pmutt.empirical.shomate.Shomate'>",
    "type": "shomate",
    "name": "H2",
    "phase": "G",
    "elements": {
     "H": 2
    },
    "notes": null,
    "smiles": null,
    "model": null,
    "misc_models": [
     {
      "class": "<class 'pmutt.empirical.GasPressureAdj'>"
     }
    ],
    "a": [
     33.066178,
     -11.363417,
     11.432816,
     -2.772874,
     -0.158558,
     -9.980797,
     172.707974,
     0.0
    ],
    "T_low": 298.0,
    "T_high": 1000.0,
    "units": "J/mol/K"
   }

Creating New pmutt Classes
---------------------------

Encoding
^^^^^^^^

To ensure your new class can be encoded using the ``pmuttEncoder``, the 
``to_dict()`` method should be implemented. One of the entries of the 
dictionary should be ``'class': str(self.__class__)`` so that it can be decoded 
later. The other elements should be the attributes that can be used to 
reinitialize the object and must be JSON-supported objects. A simple example 
using :class:`~pmutt.statmech.trans.FreeTrans` is shown below.

.. code:: python

   def to_dict(self):
       return {'class': str(self.__class__),
               'n_degrees': self.n_degrees, 
               'molecular_weight': self.molecular_weight}
   

If the attributes are not supported by JSON (such as other pmutt objects), use 
their ``to_dict()`` methods to convert to JSON-supported objects. An example 
using :class:`~pmutt.statmech.StatMech` is shown below.

.. code:: python

   def to_dict(self):
       return {'class': str(self.__class__),
               'trans_model': self.trans_model.to_dict(),
               'vib_model': self.vib_model.to_dict(),
               'rot_model': self.rot_model.to_dict(),
               'elec_model': self.elec_model.to_dict(),
               'nucl_model': self.nucl_model.to_dict()}

Decoding
^^^^^^^^

To ensure your object can be decoded using the ``json_to_pmutt`` object hook, 
add an entry to the dictionary in the ``pmutt.io.json.type_to_class`` method.
The key should be the type of your object in string format (i.e. the result of 
``str(self.__class__)``). Your class should also have the ``from_dict()`` class 
method to reinitialize your object. A simple example using 
:class:`~pmutt.statmech.trans.FreeTrans` is shown below.

.. code:: python

   from pmutt.io.json import remove_class

   @classmethod
   def from_dict(cls, json_obj):
       json_obj = remove_class(json_obj)
       return cls(**json_obj)

Similarly to encoding, sometimes your object contains pmutt objects. You can 
use the ``json_to_pmutt`` object hook to remake these objects. An example using 
:class:`~pmutt.statmech.StatMech` is shown below.

.. code:: python

   from pmutt.io import json as json_pmutt

   @classmethod
   def from_dict(cls, json_obj):
       json_obj = remove_class(json_obj)
       trans_model = json_pmutt.json_to_pmutt(json_obj['trans_model'])
       vib_model = json_pmutt.json_to_pmutt(json_obj['vib_model'])
       rot_model = json_pmutt.json_to_pmutt(json_obj['rot_model'])
       elec_model = json_pmutt.json_to_pmutt(json_obj['elec_model'])
       nucl_model = json_pmutt.json_to_pmutt(json_obj['nucl_model'])

       return cls(trans_model=trans_model, 
                  vib_model=vib_model, 
                  rot_model=rot_model,
                  elec_model=elec_model,
                  nucl_model=nucl_model)

--------------------------------------------------------------------------------

YAML
====

`YAML Ain't Markup Language (YAML)`_ is a human friendly data serialization
standard for all programming languages. All pmutt objects are natively supported
by YAML.

Examples
--------

Saving pmutt objects can be done by using 
:func:`~pmutt.io.json.pmuttEncoder`.

.. code:: python

   import yaml
   
   with open(yaml_path, 'w') as f_ptr:
       yaml.dump(pmutt_obj, f_ptr)
   
Loading pmutt objects can be done by using the object hook: 
:func:`~pmutt.io.json.json_to_pmutt`.

.. code:: python

   import yaml

   with open(yaml_path, 'r') as f_ptr:
       pmutt_obj = yaml.load(f_ptr)

Sample YAML File
----------------

YAML writes in a human-readable syntax. An example showing a H2 
:class:`~pmutt.empirical.shomate.Shomate` object in YAML format is shown below.
::

   !!python/object:pmutt.empirical.shomate.Shomate
   T_high: 1000.0
   T_low: 298.0
   _units: J/mol/K
   a: !!python/object/apply:numpy.core.multiarray._reconstruct
     args:
     - !!python/name:numpy.ndarray ''
     - !!python/tuple [0]
     - !!binary |
       Yg==
     state: !!python/tuple
     - 1
     - !!python/tuple [8]
     - !!python/object/apply:numpy.dtype
       args: [f8, 0, 1]
       state: !!python/tuple [3, <, null, null, null, -1, -1, 0]
     - false
     - !!binary |
       e9tMhXiIQEDxngPLEbomwP9eCg+a3SZA/G1PkNguBsAZdELooEvEv6MHPgYr9iPAYw0XuaeWZUAA
       AAAAAAAAAA==
   elements: {H: 2}
   misc_models:
   - !!python/object:pmutt.empirical.GasPressureAdj {}
   model: null
   n_sites: null
   name: H2
   notes: null
   phase: G
   smiles: null

--------------------------------------------------------------------------------

VASP
====

.. autosummary::
   :toctree: vasp
   :nosignatures:

   vasp.set_vib_wavenumbers_from_outcar

--------------------------------------------------------------------------------

Gaussian
========

.. autosummary::
   :toctree: gaussian
   :nosignatures:

   gaussian.read_pattern
   gaussian.read_zpe
   gaussian.read_electronic_and_zpe
   gaussian.read_frequencies
   gaussian.read_rotational_temperatures
   gaussian.read_molecular_mass
   gaussian.read_rot_symmetry_num

--------------------------------------------------------------------------------

Chemkin
=======

.. autosummary::
   :toctree: chemkin
   :nosignatures:

   chemkin.read_reactions
   chemkin.write_EA
   chemkin.write_gas
   chemkin.write_surf
   chemkin.write_T_flow
   chemkin.write_tube_mole

--------------------------------------------------------------------------------

OpenMKM
=======

.. autosummary::
   :toctree: omkm
   :nosignatures:

   omkm.write_cti
   omkm.write_yaml
   omkm.get_species_phases
   omkm.get_reactions_phases
   omkm.get_interactions_phases
   omkm.organize_phases

.. _`pmutt.examples.read_nasa_from_thermdat`: https://github.com/VlachosGroup/pmutt/tree/master/examples/read_nasa_from_thermdat 
.. _`YAML Ain't Markup Language (YAML)`: https://yaml.org/
.. _`JavaScript Object Notation (JSON)`: https://www.json.org/