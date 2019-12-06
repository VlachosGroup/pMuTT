.. _io:

Input and Output
****************

Input and output to different forms is an active area of development for pmutt

Excel
=====

.. autofunction:: pmutt.io.excel.read_excel

Special Rules
-------------

Special rules can be defined in the :func:`~pmutt.io.excel.read_excel` 
function to process inputs differently. Currently supported special rules are 
listed below.

.. autofunction:: pmutt.io.excel.set_element

.. autofunction:: pmutt.io.excel.set_formula

.. autofunction:: pmutt.io.excel.set_atoms

.. autofunction:: pmutt.io.excel.set_statmech_model

.. autofunction:: pmutt.io.excel.set_trans_model

.. autofunction:: pmutt.io.excel.set_vib_model

.. autofunction:: pmutt.io.excel.set_rot_model

.. autofunction:: pmutt.io.excel.set_elec_model

.. autofunction:: pmutt.io.excel.set_nucl_model

.. autofunction:: pmutt.io.excel.set_vib_wavenumbers

.. autofunction:: pmutt.io.excel.set_rot_temperatures

.. autofunction:: pmutt.io.excel.set_nasa_a_low

.. autofunction:: pmutt.io.excel.set_nasa_a_high

.. autofunction:: pmutt.io.excel.set_list_value

.. autofunction:: pmutt.io.excel.set_dict_value

Examples
--------

.. _DFT-Input-Example:

DFT Input Example
^^^^^^^^^^^^^^^^^^

This example uses data found in :ref:`excel_to_empirical_example`. 
Below, we show the contents of the references.xlsx spreadsheet. The first row 
corresponds to header labels. The headers may have special processing rules, 
which can be found in the docstring of :func:`~pmutt.io.read_excel`. If no 
special rules are defined, then the output dictionary will use the header as a 
key and field as a value. The second row (only shown in the Excel file) is a 
description of the header. A good description should include units, and 
supported options if the field is discrete. The subsequent rows describe 
the interested species.

+------+-------+------------+------------+----------------+-------+--------------+-----------------+-----------+---------------+----------------+------+----------------+----------------+----------------+
| name | phase | elements.H | elements.O | statmech_model | T_ref | HoRT_ref     | potentialenergy | geometry  | atoms         | symmetrynumber | spin | vib_wavenumber | vib_wavenumber | vib_wavenumber |
+======+=======+============+============+================+=======+==============+=================+===========+===============+================+======+================+================+================+
| H2   | G     | 2          | 0          | IdealGas       | 298   | 0            | -6.7598         | linear    | .\H2\CONTCAR  | 2              | 0    | 4306.1793      |                |                |
+------+-------+------------+------------+----------------+-------+--------------+-----------------+-----------+---------------+----------------+------+----------------+----------------+----------------+
| H2O  | G     | 2          | 1          | IdealGas       | 298   | -97.60604334 | -14.2209        | nonlinear | .\H2O\CONTCAR | 2              | 0    | 3825.434       | 3710.2642      | 1582.432       |
+------+-------+------------+------------+----------------+-------+--------------+-----------------+-----------+---------------+----------------+------+----------------+----------------+----------------+

The :func:`~pmutt.io.excel.read_excel` function returns a list of dictionaries. 
The dictionaries contain field-to-value pairings that can be used to initilize 
objects using the keyword argument syntax (\*\*kwargs). This is shown in code 
below:

.. code:: python

    from pprint import pprint
    from pmutt.io.excel import read_excel
    from pmutt.empirical.references import Reference, References

    refs_path = './references.xlsx'
    refs_input = read_excel(io=refs_path)
    refs = References([Reference(**ref_input) for ref_input in refs_input])

    print('Reference Input:')
    pprint(refs_input)

The output can be shown below::

    [{'atoms': Atoms(symbols='OH2', pbc=True, cell=[20.0, 21.0, 22.0]),
      'elements': {'H': 2, 'O': 1, 'Pt': 0},
      'geometry': 'nonlinear',
      'name': 'H2O',
      'phase': 'G',
      'potentialenergy': -14.2209,
      'spin': 0.0,
      'symmetrynumber': 2.0,
      'statmech_model': <class 'pmutt.statmech.idealgasthermo.IdealGasThermo'>,
      'vib_energies': [0.47429336414391626,
                       0.460014128927786,
                       0.19619656143825398]},
     {'atoms': Atoms(symbols='H2', pbc=True, cell=[20.0, 21.0, 22.0]),
      'elements': {'H': 2, 'O': 0, 'Pt': 0},
      'geometry': 'linear',
      'name': 'H2',
      'phase': 'G',
      'potentialenergy': -6.7598,
      'spin': 0.0,
      'symmetrynumber': 2.0,
      'statmech_model': <class 'pmutt.statmech.idealgasthermo.IdealGasThermo'>,
      'vib_energies': [0.5338981843116086]},
     {'atoms': Atoms(symbols='O2', pbc=True, cell=[20.0, 20.0, 20.0]),
      'elements': {'H': 0, 'O': 2, 'Pt': 0},
      'geometry': 'linear',
      'name': 'O2',
      'phase': 'G',
      'potentialenergy': -9.86,
      'spin': 1.0,
      'symmetrynumber': 2.0,
      'statmech_model': <class 'pmutt.statmech.idealgasthermo.IdealGasThermo'>,
      'vib_energies': [0.2733851552365915]},
     {'elements': {'H': 0, 'O': 1, 'Pt': 1},
      'name': 'MO(S)',
      'phase': 'S',
      'potentialenergy': 0.0,
      'statmech_model': <class 'pmutt.statmech.harmonicthermo.HarmonicThermo'>,
      'vib_energies': [0.07025434894614345,
                       0.06873635809621279,
                       0.034434367577936324]},
     {'elements': {'H': 0, 'O': 0, 'Pt': 1},
      'name': 'MO(B)',
      'phase': 'S',
      'potentialenergy': 0.0,
      'statmech_model': <class 'pmutt.statmech.harmonicthermo.HarmonicThermo'>,
      'vib_energies': [0.07025434894614345,
                       0.06873635809621279,
                       0.034434367577936324]},
     {'elements': {'H': 0, 'O': 0, 'Pt': 1},
      'name': 'V-MO(S)',
      'phase': 'S',
      'potentialenergy': 7.0,
      'statmech_model': <class 'pmutt.statmech.harmonicthermo.HarmonicThermo'>,
      'vib_energies': []},
     {'elements': {'H': 0, 'O': 1, 'Pt': 1},
      'name': 'MO_bulk(S)',
      'phase': 'S',
      'potentialenergy': 0.0,
      'statmech_model': <class 'pmutt.statmech.harmonicthermo.HarmonicThermo'>,
      'vib_energies': [0.07025434894614345,
                       0.06873635809621279,
                       0.034434367577936324]},
     {'elements': {'H': 0, 'O': 0, 'Pt': 1},
      'name': 'MO_bulk(B)',
      'phase': 'S',
      'potentialenergy': 0.0,
      'statmech_model': <class 'pmutt.statmech.harmonicthermo.HarmonicThermo'>,
      'vib_energies': [0.07025434894614345,
                       0.06873635809621279,
                       0.034434367577936324]},
     {'elements': {'H': 0, 'O': 0, 'Pt': 1},
      'name': 'V-MO_bulk(S)',
      'phase': 'S',
      'potentialenergy': 7.0,
      'statmech_model': <class 'pmutt.statmech.harmonicthermo.HarmonicThermo'>,
      'vib_energies': []}]

NASA Polynomial Input Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Due to the special rules defined for NASA parsing, a group of NASA polynomials
can be directly imported using :func:`~pmutt.io.excel.read_excel`.

+------+-------+------------+------------+-------+-------+--------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+--------+
| name | phase | elements.H | elements.O | T_low | T_mid | T_high | nasa.a_low.0 | nasa.a_low.1 | nasa.a_low.2 | nasa.a_low.3 | nasa.a_low.4 | nasa.a_low.5 | nasa.a_low.6 | nasa.a_high.0 | nasa.a_high.1 | nasa.a_high.2 | nasa.a_high.3 | nasa.a_high.4 | nasa.a_high.5 | nasa.a_high.6 | notes  |
+======+=======+============+============+=======+=======+========+==============+==============+==============+==============+==============+==============+==============+===============+===============+===============+===============+===============+===============+===============+========+
| O2   | G     |            | 2          | 200   | 1000  | 3500   | 3.78E+00     | -3.00E-03    | 9.85E-06     | -9.68E-09    | 3.24E-12     | -1.06E+03    | 3.66E+00     | 3.28E+00      | 1.48E-03      | -7.58E-07     | 2.09E-10      | -2.17E-14     | -1.09E+03     | 5.45E+00      | TPIS89 |
+------+-------+------------+------------+-------+-------+--------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+--------+
| H2   | G     | 2          |            | 200   | 1000  | 3500   | 2.34E+00     | 7.98E-03     | -1.95E-05    | 2.02E-08     | -7.38E-12    | -9.18E+02    | 6.83E-01     | 3.34E+00      | -4.94E-05     | 4.99E-07      | -1.80E-10     | 2.00E-14      | -950.158922   | -3.20502331   | TPIS78 |
+------+-------+------------+------------+-------+-------+--------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+--------+

Using a similar set of statements as the :ref:`DFT-Input-Example`.

.. code:: python

    from pmutt.io.excel import read_excel
    from pmutt.empirical.nasa import Nasa
    
    species_data = read_excel('input_data.xlsx')
    species = [Nasa(**specie_data) for specie_data in species_data]
    
    pprint(species_data)

The above code gives the following output::

    [{'T_high': 3500,
      'T_low': 200,
      'T_mid': 1000,
      'a_high': array([ 3.28253784e+00,  1.48308754e-03, -7.57966669e-07,  2.09470555e-10,
           -2.16717794e-14, -1.08845772e+03,  5.45323129e+00]),
      'a_low': array([ 3.78245636e+00, -2.99673416e-03,  9.84730201e-06, -9.68129509e-09,
            3.24372837e-12, -1.06394356e+03,  3.65767573e+00]),
      'elements': {'O': 2.0},
      'name': 'O2',
      'notes': 'TPIS89',
      'phase': 'G',
      'vib_energies': []},
     {'T_high': 3500,
      'T_low': 200,
      'T_mid': 1000,
      'a_high': array([ 3.33727920e+00, -4.94024731e-05,  4.99456778e-07, -1.79566394e-10,
            2.00255376e-14, -9.50158922e+02, -3.20502331e+00]),
      'a_low': array([ 2.34433112e+00,  7.98052075e-03, -1.94781510e-05,  2.01572094e-08,
           -7.37611761e-12, -9.17935173e+02,  6.83010238e-01]),
      'elements': {'H': 2.0},
      'name': 'H2',
      'notes': 'TPIS78',
      'phase': 'G',
      'vib_energies': []}]

Thermdat
========
This is the output format used for Chemkin. A list of NASA objects can be 
written to a thermdat file.

.. automodule:: pmutt.io.thermdat
   :members:

Examples
--------
Reading Thermdat
^^^^^^^^^^^^^^^^
A thermdat file can be read directly by using 
:func:`~pmutt.io.thermdat.read_thermdat`. See the :ref:`overview_example` of
how its done.

.. code:: python

    import os
    from pprint import pprint
    from matplotlib import pyplot as plt
    from pmutt.io.thermdat import read_thermdat
    from pmutt.empirical.nasa import Nasa

    base_path = os.path.dirname(__file__)
    #Thermdat file from http://combustion.berkeley.edu/gri_mech/version30/files30/thermo30.dat
    species = read_thermdat('{}/thermdat'.format(base_path))

    #Printing information related to each specie
    for specie in species:
        print('Name: {}'.format(specie.name))
        for key, val in specie.__dict__.items():
            if key != 'name':
                print('\t{}\t{}'.format(key, val))

    #Plot an example of an imported NASA polynomial
    species[1].plot_empirical()
    plt.show()

A snippet of the species information printed is shown below::

    Name: AR
        phase   G
        elements        {'AR': 1}
        T_ref   298.15
        references      None
        notes   120186
        statmech_model    None
        HoRT_dft        None
        HoRT_ref        None
        a_low   [   2.5      0.       0.       0.       0.    -745.375    4.366]
        a_high  [   2.5      0.       0.       0.       0.    -745.375    4.366]
        T_low   300.0
        T_high  5000.0
        T_mid   1000.0

And a sample plot is shown below

.. image:: read_nasa_from_thermdat_example_O2.png

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

JSON
====
`JavaScript Object Notation (JSON)`_ is a format that is easily read and
written by both humans and machines. All pmutt objects with the methods 
``to_dict`` and ``from_dict`` are JSON compatible.

.. automodule:: pmutt.io.json
   :members:

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

VASP
====

.. autofunction:: pmutt.io.vasp.set_vib_wavenumbers_from_outcar

Gaussian
========

.. automodule:: pmutt.io.gaussian
   :members:

Chemkin
=======

.. automodule:: pmutt.io.chemkin
   :members:

OpenMKM
=======

.. automodule:: pmutt.io.omkm
   :members:

.. _`pmutt.examples.read_nasa_from_thermdat`: https://github.com/VlachosGroup/pmutt/tree/master/examples/read_nasa_from_thermdat 
.. _`JavaScript Object Notation (JSON)`: https://www.json.org/