.. _examples:

Examples
========

This page and its subpages show some examples using the pMuTT code. There are
several formats available. 

Running the example with Binder only requires the user to click the logo.
For the other formats, a page of text may appear after clicking the link.
Right-click and choose "Save Page As" or find the file in the ZIP folder.

Downloading the ZIP folder is recommended because the code may not run without
other files (e.g. spreadsheets).


The pros and cons of each format are listed below.

+------------------+----------------------------+-----------------------------+----------------------------------------+------------------------------------------------------------------------+--------------------------------+----------------------+
| Format           | File located in ZIP folder | pMuTT installation required | Jupyter Notebook installation required | Code may require additional files to be downloaded (e.g. spreadsheets) | Non-code elements easy to read | Code easily editable |
+==================+============================+=============================+========================================+========================================================================+================================+======================+
| Binder           | N                          | N                           | N                                      | N                                                                      | Y                              | Y                    |
+------------------+----------------------------+-----------------------------+----------------------------------------+------------------------------------------------------------------------+--------------------------------+----------------------+
| HTML             | Y                          | N                           | N                                      | N                                                                      | Y                              | N                    |
+------------------+----------------------------+-----------------------------+----------------------------------------+------------------------------------------------------------------------+--------------------------------+----------------------+
| Jupyter Notebook | Y                          | Y                           | Y                                      | Y                                                                      | Y                              | Y                    |
+------------------+----------------------------+-----------------------------+----------------------------------------+------------------------------------------------------------------------+--------------------------------+----------------------+
| Python Script    | Y                          | Y                           | N                                      | Y                                                                      | N                              | Y                    |
+------------------+----------------------------+-----------------------------+----------------------------------------+------------------------------------------------------------------------+--------------------------------+----------------------+

.. _overview_example:

Overview
--------

.. image:: ./images/binder_logo.png
  :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Foverview%2Foverview.ipynb
  :width: 15%
.. image:: ./images/html_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/overview/overview.html
  :width: 15%
.. image:: ./images/jupyter_notebook.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/overview/overview.ipynb
  :width: 15%
.. image:: ./images/python_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/overview/overview.py
  :width: 15%
.. image:: ./images/zip_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/overview/overview.zip
  :width: 15%

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

--------------------------------------------------------------------------------

.. _expt_to_empirical_example:

Experimental to Empirical
-------------------------

.. image:: ./images/binder_logo.png
  :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Fexpt_data_to_empirical%2Fexpt_data_to_empirical_object.ipynb
  :width: 15%
.. image:: ./images/html_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/expt_data_to_empirical/expt_data_to_empirical_object.html
  :width: 15%
.. image:: ./images/jupyter_notebook.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/expt_data_to_empirical/expt_data_to_empirical_object.ipynb
  :width: 15%
.. image:: ./images/python_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/expt_data_to_empirical/expt_data_to_empirical_object.py
  :width: 15%
.. image:: ./images/zip_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/expt_data_to_empirical/expt_data_to_empirical.zip
  :width: 15%

Topics Covered
^^^^^^^^^^^^^^

- Using pmutt's constants for unit conversions
- Create a :class:`~pmutt.empirical.shomate.Shomate` object from experimental data
- Calculate thermodynamic properties using the :class:`~pmutt.empirical.shomate.Shomate` object
- Plot the shape of the :class:`~pmutt.empirical.shomate.Shomate` curve
- Save the :class:`~pmutt.empirical.shomate.Shomate` object as a JSON file

--------------------------------------------------------------------------------

.. _excel_to_empirical_example:

Excel to Empirical Data
-----------------------

.. image:: ./images/binder_logo.png
  :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Fexcel_to_empirical%2Fexcel_to_empirical.ipynb
  :width: 15%
.. image:: ./images/html_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/excel_to_empirical/excel_to_empirical.html
  :width: 15%
.. image:: ./images/jupyter_notebook.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/excel_to_empirical/excel_to_empirical.ipynb
  :width: 15%
.. image:: ./images/python_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/excel_to_empirical/excel_to_empirical.py
  :width: 15%
.. image:: ./images/zip_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/excel_to_empirical/excel_to_empirical.zip
  :width: 15%

Topics Covered
^^^^^^^^^^^^^^

- Reading *ab-initio* data from an Excel file
- Initialize :class:`~pmutt.empirical.references.Reference` objects and a :class:`~pmutt.empirical.references.References` object
- Generate a :class:`~pmutt.empirical.nasa.Nasa` object using :class:`~pmutt.statmech.StatMech` models
- Write :class:`~pmutt.empirical.nasa.Nasa` object to a ``thermdat`` file

--------------------------------------------------------------------------------

.. _reaction_example:

Reaction
--------

.. image:: ./images/binder_logo.png
  :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Freactions%2Freactions.ipynb
  :width: 15%
.. image:: ./images/html_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/reactions/reactions.html
  :width: 15%
.. image:: ./images/jupyter_notebook.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/reactions/reactions.ipynb
  :width: 15%
.. image:: ./images/python_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/reactions/reactions.py
  :width: 15%
.. image:: ./images/zip_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/reactions/reactions.zip
  :width: 15%

Topics Covered
^^^^^^^^^^^^^^

- Read a thermdat file and convert it to a dictionary of Nasa objects
- Initialize a :class:`~pmutt.reaction.Reaction` object manually and from strings
- Add a BEP relationship to a :class:`~pmutt.reaction.Reaction` object
- Calculate thermodynamic and kinetic properties using the 
  :class:`~pmutt.reaction.Reaction` object
- Save the :class:`~pmutt.reaction.Reaction` object as a ``JSON`` file

--------------------------------------------------------------------------------

.. _chemkin_io_example:

Chemkin_IO
----------

.. image:: ./images/binder_logo.png
  :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Fchemkin_io%2FChemkin_IO.ipynb
  :width: 15%
.. image:: ./images/html_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/chemkin_io/Chemkin_IO.html
  :width: 15%
.. image:: ./images/jupyter_notebook.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/chemkin_io/Chemkin_IO.ipynb
  :width: 15%
.. image:: ./images/python_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/chemkin_io/Chemkin_IO.py
  :width: 15%
.. image:: ./images/zip_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/chemkin_io/chemkin_io.zip
  :width: 15%

Topics Covered
^^^^^^^^^^^^^^
- Read species *ab-initio* data, reactions, and catalyst sites from a
  spreadsheet
- Write the thermdat, gas.inp, surf.inp, T_flow.inp, EAg.inp, EAs.inp, 
  tube_mole.inp files

--------------------------------------------------------------------------------

.. _openmkm_io_example:

OpenMKM_IO
----------

.. image:: ./images/binder_logo.png
  :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Fomkm_io%2FOpenMKM_IO.ipynb
  :width: 15%
.. image:: ./images/html_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/omkm_io/OpenMKM_IO.html
  :width: 15%
.. image:: ./images/jupyter_notebook.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/omkm_io/OpenMKM_IO.ipynb
  :width: 15%
.. image:: ./images/python_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/omkm_io/OpenMKM_IO.py
  :width: 15%
.. image:: ./images/zip_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/omkm_io/omkm_io.zip
  :width: 15%

Topics Covered
^^^^^^^^^^^^^^
- Read species *ab-initio* data, reactions, lateral interactions and phases
  from a spreadsheet
- Write the CTI input file

--------------------------------------------------------------------------------

.. _phase_diagram_example:

Phase Diagram
-------------

.. image:: ./images/binder_logo.png
  :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Fphase_diagram%2FPhaseDiagram.ipynb
  :width: 15%
.. image:: ./images/html_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/phase_diagram/PhaseDiagram.html
  :width: 15%
.. image:: ./images/jupyter_notebook.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/phase_diagram/PhaseDiagram.ipynb
  :width: 15%
.. image:: ./images/python_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/phase_diagram/PhaseDiagram.py
  :width: 15%
.. image:: ./images/zip_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/phase_diagram/phase_diagram.zip
  :width: 15%

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

.. _nam26_workshop_example:

NAM26 Workshop
--------------

.. image:: ./images/binder_logo.png
  :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Fnam2019%2FNAM_2019_Workshop.ipynb
  :width: 15%
.. image:: ./images/html_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/nam2019/NAM_2019_Workshop.html
  :width: 15%
.. image:: ./images/jupyter_notebook.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/nam2019/NAM_2019_Workshop.ipynb
  :width: 15%
.. image:: ./images/python_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/nam2019/NAM_2019_Workshop.py
  :width: 15%
.. image:: ./images/zip_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/nam2019/nam2019.zip
  :width: 15%

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

.. aiche_2019_workshop_example:

AIChE 2019 Workshop
-------------------

.. image:: ./images/binder_logo.png
  :target: https://mybinder.org/v2/gh/VlachosGroup/pmutt/master?filepath=docs%2Fsource%2Fexamples_jupyter%2Faiche2019%2Fpmutt_aiche2019.ipynb
  :width: 15%
.. image:: ./images/html_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/aiche2019/pmutt_aiche2019.html
  :width: 15%
.. image:: ./images/jupyter_notebook.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/aiche2019/pmutt_aiche2019.ipynb
  :width: 15%
.. image:: ./images/python_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/aiche2019/pmutt_aiche2019.py
  :width: 15%
.. image:: ./images/zip_logo.png
  :target: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/examples_jupyter/aiche2019/aiche2019.zip
  :width: 15%

.. _`examples folder`: https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter
.. _`Jupyter Notebook`: https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/