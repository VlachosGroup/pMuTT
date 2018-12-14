.. _examples:

Examples
========
This page and its subpages show some examples using the pMuTT code. All the 
codes listed here are available in the `examples folder`_ as Jupyter notebooks 
(*.ipynb) and ordinary Python scripts (*.py). 

Experimental to Empirical
-------------------------
- `Webpage <examples_jupyter/expt_data_to_empirical/expt_data_to_empirical_object.ipynb>`_
- `Jupyter Notebook <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/expt_data_to_empirical/expt_data_to_empirical_object.ipynb>`_
- `Python Script <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/expt_data_to_empirical/expt_data_to_empirical_object.py>`_

Topics Covered
^^^^^^^^^^^^^^

- Using pMuTT's constants for unit conversions
- Create a :class:`~pMuTT.empirical.shomate.Shomate` object from experimental data
- Calculate thermodynamic properties using the :class:`~pMuTT.empirical.shomate.Shomate` object
- Plot the shape of the :class:`~pMuTT.empirical.shomate.Shomate` curve
- Save the :class:`~pMuTT.empirical.shomate.Shomate` object as a JSON file


Excel to Empirical Data
-----------------------
- `Webpage <examples_jupyter/excel_to_empirical/excel_to_empirical.ipynb>`_
- `Jupyter Notebook <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/excel_to_empirical/excel_to_empirical.ipynb>`_
- `Python Script <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/excel_to_empirical/excel_to_empirical.py>`_

Topics Covered
^^^^^^^^^^^^^^

- Reading *ab-initio* data from an Excel file
- Initialize :class:`~pMuTT.empirical.references.Reference` objects and a :class:`~pMuTT.empirical.references.References` object
- Generate a :class:`~pMuTT.empirical.nasa.Nasa` object using :class:`~pMuTT.statmech.StatMech` models
- Write :class:`~pMuTT.empirical.nasa.Nasa` object to a ``thermdat`` file

Reaction
--------
- `Webpage <examples_jupyter/reactions/reactions.ipynb>`_
- `Jupyter Notebook <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/reactions/reactions.ipynb>`_
- `Python Script <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/reactions/reactions.py>`_

Topics Covered
^^^^^^^^^^^^^^

- Read a thermdat file and convert it to a dictionary of Nasa objects
- Initialize a :class:`~pMuTT.reaction.Reaction` object manually and from strings
- Add a BEP relationship to a :class:`~pMuTT.reaction.Reaction` object
- Calculate thermodynamic and kinetic properties using the 
  :class:`~pMuTT.reaction.Reaction` object
- Save the :class:`~pMuTT.reaction.Reaction` object as a ``JSON`` file

.. _`examples folder`: https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter