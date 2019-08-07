.. _examples:

Examples
========
This page and its subpages show some examples using the pMuTT code. All the 
codes listed here are available in the `examples folder`_ as Binder notebooks,
Jupyter notebooks (\*.ipynb) and ordinary Python scripts (\*.py).

Clicking on the binder badge allows you to run the Jupyter notebook without
needing to install pMuTT. If you would prefer to open the Jupyter notebook
locally, `follow these instructions`_.

Overview
--------
- .. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Foverview%2Foverview.ipynb
- `Jupyter Notebook (Overview) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/overview/overview.ipynb>`_
- `Python Script (Overview) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/overview/overview.py>`_


Topics Covered
^^^^^^^^^^^^^^

- Using constants and converting units using the :mod:`~pmutt.constants` module
- Initializing :class:`~pmutt.statmech.StatMech` objects by specifying all modes
  and by using :ref:`presets`
- Initializing empirical objects such as :class:`~pmutt.empirical.nasa.Nasa`
  objects using a :class:`~pmutt.statmech.StatMech` object or from a previously
  generated Nasa polynomial
- Initializing :class:`~pmutt.empirical.references.Reference` and
  :class:`~pmutt.empirical.references.References` objects to adjust DFT's
  reference to more traditional references
- Input (via Excel) and output :class:`~pmutt.empirical.nasa.Nasa` polynomials
  to thermdat format
- Initializing :class:`~pmutt.reaction.Reaction` objects from strings

------------

Experimental to Empirical
-------------------------
- .. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Fexpt_data_to_empirical%2Fexpt_data_to_empirical_object.ipynb
- `Jupyter Notebook (Experimental to Empirical) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/expt_data_to_empirical/expt_data_to_empirical_object.ipynb>`_
- `Python Script (Experimental to Empirical) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/expt_data_to_empirical/expt_data_to_empirical_object.py>`_

Topics Covered
^^^^^^^^^^^^^^

- Using pmutt's constants for unit conversions
- Create a :class:`~pmutt.empirical.shomate.Shomate` object from experimental data
- Calculate thermodynamic properties using the :class:`~pmutt.empirical.shomate.Shomate` object
- Plot the shape of the :class:`~pmutt.empirical.shomate.Shomate` curve
- Save the :class:`~pmutt.empirical.shomate.Shomate` object as a JSON file

------------

Excel to Empirical Data
-----------------------
- .. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Fexcel_to_empirical%2Fexcel_to_empirical.ipynb
- `Jupyter Notebook (Excel to Empirical Data) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/excel_to_empirical/excel_to_empirical.ipynb>`_
- `Python Script (Excel to Empirical Data) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/excel_to_empirical/excel_to_empirical.py>`_

Topics Covered
^^^^^^^^^^^^^^

- Reading *ab-initio* data from an Excel file
- Initialize :class:`~pmutt.empirical.references.Reference` objects and a :class:`~pmutt.empirical.references.References` object
- Generate a :class:`~pmutt.empirical.nasa.Nasa` object using :class:`~pmutt.statmech.StatMech` models
- Write :class:`~pmutt.empirical.nasa.Nasa` object to a ``thermdat`` file

------------

Reaction
--------
- .. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Freactions%2Freactions.ipynb
- `Jupyter Notebook (Reaction) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/reactions/reactions.ipynb>`_
- `Python Script (Reaction) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/reactions/reactions.py>`_

Topics Covered
^^^^^^^^^^^^^^

- Read a thermdat file and convert it to a dictionary of Nasa objects
- Initialize a :class:`~pmutt.reaction.Reaction` object manually and from strings
- Add a BEP relationship to a :class:`~pmutt.reaction.Reaction` object
- Calculate thermodynamic and kinetic properties using the 
  :class:`~pmutt.reaction.Reaction` object
- Save the :class:`~pmutt.reaction.Reaction` object as a ``JSON`` file

------------

Chemkin_IO
----------
- .. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Fchemkin_io%2FChemkin_IO.ipynb
- `Jupyter Notebook (Chemkin_IO) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/chemkin_io/Chemkin_IO.ipynb>`_
- `Python Script (Chemkin_IO) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/chemkin_io/Chemkin_IO.py>`_

Topics Covered
^^^^^^^^^^^^^^
- Read species *ab-initio* data, reactions, and catalyst sites from a
  spreadsheet
- Write the thermdat, gas.inp, surf.inp, T_flow.inp, EAg.inp, EAs.inp, 
  tube_mole.inp files

------------

Phase Diagram
-------------
- .. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Fphase_diagram%2FPhaseDiagram.ipynb
- `Jupyter Notebook (Phase Diagram) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/phase_diagram/PhaseDiagram.ipynb>`_
- `Python Script (Phase Diagram) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/phase_diagram/PhaseDiagram.py>`_

Topics Covered
^^^^^^^^^^^^^^

- Create :class:`~pmutt.empirical.nasa.Nasa` and 
  :class:`~pmutt.statmech.StatMech` objects 
- Initialize :class:`~pmutt.reaction.Reaction` objects to describe the 
  formation reaction of FeOx species
- Generate a 1D phase diagram by varying T
- Generate a 2D phase diagram by varying T and P
- Save the :class:`~pmutt.reaction.phasediagram.PhaseDiagram` object as a 
  ``JSON`` file

Workshops
---------

NAM26 Workshop
--------------
- .. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Fnam2019%2FNAM_2019_Workshop.ipynb
- `Jupyter Notebook (NAM26 Workhop) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/nam2019/NAM_2019_Workshop.ipynb>`_
- `Python Script (NAM26 Workhop) <https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/nam2019/NAM_2019_Workshop.py>`_

Topics Covered
^^^^^^^^^^^^^^

- Using constants and converting units using the :mod:`~pmutt.constants` module
- Initializing :class:`~pmutt.statmech.StatMech` objects by specifying all modes
  and by using :ref:`presets`
- Initializing empirical objects such as :class:`~pmutt.empirical.nasa.Nasa`
  objects using a :class:`~pmutt.statmech.StatMech` object or directly using 
  the polynomial
- Input (via Excel) and output :class:`~pmutt.empirical.nasa.Nasa` polynomials
  to thermdat format
- Initializing :class:`~pmutt.reaction.Reaction` objects from strings

.. _`examples folder`: https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter
.. _`follow these instructions`: https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/