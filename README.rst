Python Multiscale Thermochemistry Toolbox (pMuTT)
=================================================
  
The **P**\ ython **Mu**\ ltiscale **T**\ hermochemistry **T**\ oolbox
(pMuTT) is a Python library for Thermochemistry developed by the
Vlachos Research Group at the University of Delaware. This code was
originally developed to convert *ab-initio* data from DFT to observable
thermodynamic properties such as heat capacity, enthalpy, entropy, and
Gibbs energy. These properties can be fit to empirical equations and
written to different formats. 

.. image:: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/logos/pmutt_web.png
   :target: https://vlachosgroup.github.io/pMuTT/
   :align: center

Documentation
-------------

See our `documentation page`_ for examples, equations used, and docstrings.

Developers
----------

-  Jonathan Lym, Ph.D. (jonathanalym@gmail.com)
-  Gerhard Wittreich, Ph.D, P.E. (wittregr@udel.edu)

Dependencies
------------

-  Python3
-  `Atomic Simulation Environment`_: Used for I/O operations and to
   calculate some thermodynamic properties
-  `Numpy`_: Used for vector and matrix operations
-  `Pandas`_: Used to import data from Excel files
-  `xlrd`_: Used by Pandas to import Excel files
-  `SciPy`_: Used for fitting heat capacities and generating smooth curves for
   reaction coordinate diagram
-  `Matplotlib`_: Used for plotting thermodynamic data
-  `pyGal`_: Similar to Matplotlib. Used for plotting interactive graphs
-  `PyMongo`_: Used to read/write to databases
-  `dnspython`_: Used to connect to databases
-  `NetworkX`_: Used to plot reaction networks
-  `More Itertools`_: Used for writing ranges for OpenMKM output.
-  `PyYAML`_: Used to write YAML input files for OpenMKM.
-  `pymatgen`_: Used to determine symmetry number for gas molecules.

Getting Started
---------------

1. Install using pip (`see documentation for more thorough instructions`_)::

    pip install pmutt

2. Look at `examples using the code`_

3. Run the `unit tests`_.

License
-------

This project is licensed under the MIT License - see the `LICENSE.md`_
file for details.

Publications
------------

- J. Lym, G.R. Wittreich and D.G. Vlachos, A Python Multiscale Thermochemistry
  Toolbox (pMuTT) for thermochemical and kinetic parameter estimation, Computer
  Physics Communications (2019) 106864,
  https://doi.org/10.1016/j.cpc.2019.106864.

Contributing
------------

If you have a suggestion or find a bug, please post to our `Issues page`_ with 
the |enhancement_label| or |bug_label| tag respectively.

Finally, if you would like to add to the body of code, please:

- fork the development branch
- make the desired changes
- write the appropriate unit tests
- submit a `pull request`_.

Questions
---------

If you are having issues, please post to our `Issues page`_ with the 
|help_wanted_label| or |question_label| tag. We will do our best to assist.

Funding
-------

This material is based upon work supported by the Department of Energy's Office 
of Energy Efficient and Renewable Energy's Advanced Manufacturing Office under 
Award Number DE-EE0007888-9.5.

Special Thanks
--------------

-  Dr. Jeffrey Frey (pip and conda compatibility)
-  Jaynell Keely (Logo design)


.. |bug_label| image:: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/images/labels/bug_small.png
   :height: 20
   :target: https://github.com/VlachosGroup/pMuTT/issues?utf8=%E2%9C%93&q=label%3Abug

.. |enhancement_label| image:: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/images/labels/enhancement_small.png
   :height: 20
   :target: https://github.com/VlachosGroup/pMuTT/issues?utf8=%E2%9C%93&q=label%3Aenhancement

.. |help_wanted_label| image:: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/images/labels/help_wanted_small.png
   :height: 20
   :target: https://github.com/VlachosGroup/pMuTT/issues?utf8=%E2%9C%93&q=label%3A%22help%20wanted%22

.. |question_label| image:: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/images/labels/question_small.png
   :height: 20
   :target: https://github.com/VlachosGroup/pMuTT/issues?utf8=%E2%9C%93&q=label%3Aquestion

.. _`documentation page`: https://vlachosgroup.github.io/pMuTT/
.. _Atomic Simulation Environment: https://wiki.fysik.dtu.dk/ase/
.. _Numpy: http://www.numpy.org/
.. _Pandas: https://pandas.pydata.org/
.. _xlrd: https://xlrd.readthedocs.io/en/latest/
.. _SciPy: https://www.scipy.org/
.. _Matplotlib: https://matplotlib.org/
.. _pyGal: http://www.pygal.org/en/stable/
.. _PyMongo: http://api.mongodb.com/python/current/
.. _dnspython: http://www.dnspython.org/
.. _networkx: https://networkx.github.io/
.. _`More Itertools`: https://more-itertools.readthedocs.io/en/stable/index.html
.. _`PyYAML`: https://pyyaml.org
.. _`pymatgen`: https://pymatgen.org/
.. _tests directory: https://github.com/VlachosGroup/pMuTT/tree/master/pmutt/tests
.. _LICENSE.md: https://github.com/VlachosGroup/pMuTT/blob/master/LICENSE.md
.. _`see documentation for more thorough instructions`: https://vlachosgroup.github.io/pMuTT/install.html
.. _`examples using the code`: https://vlachosgroup.github.io/pMuTT/examples_jupyter/examples.html
.. _`Issues page`: https://github.com/VlachosGroup/pMuTT/issues
.. _`pull request`: https://github.com/VlachosGroup/pMuTT/pulls
.. _`unit tests`: https://vlachosgroup.github.io/pMuTT/install.html#running-unit-tests
