Input and Output
----------------

Reading from and writing to different formats is an active area for
development. These operations can be found in ```PyMuTT.io_```_. We will
discuss reading DFT data via Excel.

Excel
~~~~~

A spreadsheet with DFT data can be read using
```PyMuTT.io_.excel.read_excel```_ function. It uses the Pandas function
```pandas.read_excel```_ to facilitate flexible reading options (such as
rows to skip or what to recognize as a NaN).

Below, we show the contents of the references.xlsx spreadsheet found in
```PyMuTT.examples.VASP_to_thermdat.example1```_ (the description of the
fields have been omitted here due to difficult of formatting). The first
row corresponds to the field name and the second row by default is
ignored so comments can be inputted here. Some fields do not require any
additional processing, like ``name``, ``potentialenergy``, and
``geometry``. However, in some cases the field does need further
processing, like ``atoms`` (the value is only the location to find a
file to read and not the actual atoms object) or ``elements~H`` (in this
case, the field name needs to be processed to show this header
represents hydrogen atoms in a formula unit). Special parsing
instructions can be found and added to the
```PyMuTT.io_.excel.read_excel`` <https://github.com/VlachosGroup/PyMuTT/blob/master/io_/excel.py#L45>`__
function.

+---+---+---+---+----+---+---+-----+---+----+-----+---+-----+-----+-----+
| n | p | e | e | th | T | H | pot | g | at | sym | s | vib | vib | vib |
| a | h | l | l | er | _ | o | ent | e | om | met | p | _wa | _wa | _wa |
| m | a | e | e | mo | r | R | ial | o | s  | ryn | i | ven | ven | ven |
| e | s | m | m | _m | e | T | ene | m |    | umb | n | umb | umb | umb |
|   | e | e | e | od | f | _ | rgy | e |    | er  |   | er~ | er~ | er~ |
|   |   | n | n | el |   | r |     | t |    |     |   | 1   | 2   | 3   |
|   |   | t | t |    |   | e |     | r |    |     |   |     |     |     |
|   |   | s | s |    |   | f |     | y |    |     |   |     |     |     |
|   |   | ~ | ~ |    |   |   |     |   |    |     |   |     |     |     |
|   |   | H | O |    |   |   |     |   |    |     |   |     |     |     |
+===+===+===+===+====+===+===+=====+===+====+=====+===+=====+=====+=====+
| H | G | 2 | 0 | Id | 2 | 0 | -6. | l | .: | 2   | 0 | 430 |     |     |
| 2 |   |   |   | ea | 9 |   | 759 | i | ra |     |   | 6.1 |     |     |
|   |   |   |   | lG | 8 |   | 8   | n | w- |     |   | 793 |     |     |
|   |   |   |   | as |   |   |     | e | la |     |   |     |     |     |
|   |   |   |   |    |   |   |     | a | te |     |   |     |     |     |
|   |   |   |   |    |   |   |     | r | x: |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | `\ |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | H2 |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | `: |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | ra |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | w- |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | la |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | te |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | x: |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | `\ |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | CO |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | NT |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | CA |     |   |     |     |     |
|   |   |   |   |    |   |   |     |   | R` |     |   |     |     |     |
+---+---+---+---+----+---+---+-----+---+----+-----+---+-----+-----+-----+

| H2O \| G \| 2 \| 1 \| IdealGas \| 298 \| -97.606 \| -14.2209 \|
  nonlinear \|
  .:raw-latex:`\H2O\CONTCAR | 2              | 0    | 3825.434         | 3710.2642        | 1582.432         |`

The ``PyMuTT.io_.excel.read_excel`` function returns a list of
dictionaries. The dictionaries contain field-to-value pairings that can
be used to initilize objects using the keyword argument syntax
(**kwargs). This is shown in code below:

\```python from pprint import pprint from PyMuTT.io_.excel import
read_excel from PyMuTT.models.empirical import BaseThermo from
PyMuTT.models.empirical.references import References

refs_path = ‘./references.xlsx’ refs_input = read_excel(io=refs_path)
refs = References([BaseThermo(**ref_input) for ref_input in refs_in

.. _``PyMuTT.io_``: https://github.com/VlachosGroup/PyMuTT/tree/master/io_
.. _``PyMuTT.io_.excel.read_excel``: https://github.com/VlachosGroup/PyMuTT/blob/master/io_/excel.py
.. _``pandas.read_excel``: https://pandas.pydata.org/pandas-docs/version/0.23/generated/pandas.read_excel.html
.. _``PyMuTT.examples.VASP_to_thermdat.example1``: https://github.com/VlachosGroup/PyMuTT/tree/master/examples/VASP_to_thermdat/example1