Python Multiscale Thermochemistry Toolbox (pmutt)
==================================================
  
The **P**\ ython **Mu**\ ltiscale **T**\ hermochemistry **T**\ oolbox
(pmutt) is a Python library for Thermochemistry developed by the
Vlachos Research Group at the University of Delaware. This code was
originally developed to convert *ab-initio* data from DFT to observable
thermodynamic properties such as heat capacity, enthalpy, entropy, and
Gibbs energy. These properties can be fit to empirical equations and
written to different formats. 

Documentation
-------------
See our `documentation page`_ for examples, equations used, and docstrings.

Developers
----------

-  Jonathan Lym (jlym@udel.edu)
-  Gerhard Wittreich, P.E. (wittregr@udel.edu)

Dependencies
------------

-  Python3
-  `Atomic Simulation Environment`_: Used for I/O operations and to
   calculate thermodynamic properties
-  `Numpy`_: Used for vector and matrix operations
-  `Pandas`_: Used to import data from Excel files
-  `SciPy`_: Used for fitting heat capacities and generating smooth curves for
   reaction coordinate diagram
-  `Matplotlib`_: Used for plotting thermodynamic data
-  `pyGal`_: Similar to Matplotlib. Used for plotting interactive graphs.
-  `PyMongo`_: Used to read/write to databases
-  `dnspython`_: Used to connect to databases
-  `NetworkX`_: Used to plot reaction networks.

Getting Started
---------------
1. Install using pip::

    pip install --user pmutt
   
2. Run the tests by navigating to the `tests directory`_ in a
   command-line interface and inputting the following command::

    python -m unittest

The expected output is shown below. The number of tests will not
necessarily be the same. ::

    .........................
    ----------------------------------------------------------------------
    Ran 25 tests in 0.020s

    OK

3. Look at `examples using the code`_

License
-------

This project is licensed under the MIT License - see the `LICENSE.md`_
file for details.

Contributing
------------

If you have a suggestion or find a bug, please post to our `Issues page`_ with 
the ``enhancement`` or ``bug`` tag respectively.

Finally, if you would like to add to the body of code, please:

- fork the development branch
- make the desired changes
- write the appropriate unit tests
- submit a `pull request`_.

Questions
---------

If you are having issues, please post to our `Issues page`_ with the 
``help wanted`` or ``question`` tag. We will do our best to assist.

Funding
-------

This material is based upon work supported by the Department of Energy's Office 
of Energy Efficient and Renewable Energy's Advanced Manufacturing Office under 
Award Number DE-EE0007888-9.5.

Special Thanks
--------------

-  Dr. Jeffrey Frey (pip and conda compatibility)
-  Jaynell Keely (Logo design)

.. _`documentation page`: https://vlachosgroup.github.io/pmutt/
.. _Atomic Simulation Environment: https://wiki.fysik.dtu.dk/ase/
.. _Numpy: http://www.numpy.org/
.. _Pandas: https://pandas.pydata.org/
.. _SciPy: https://www.scipy.org/
.. _Matplotlib: https://matplotlib.org/
.. _`pyGal`: http://pygal.org/en/stable/
.. _PyMongo: http://api.mongodb.com/python/current/
.. _dnspython: http://www.dnspython.org/
.. _networkx: https://networkx.github.io/
.. _tests directory: https://github.com/VlachosGroup/pmutt/tree/master/pmutt/tests
.. _LICENSE.md: https://github.com/VlachosGroup/pmutt/blob/master/LICENSE.md
.. _`examples using the code`: https://vlachosgroup.github.io/pmutt/examples.html
.. _`Issues page`: https://github.com/VlachosGroup/pmutt/issues
.. _`pull request`: https://github.com/VlachosGroup/pmutt/pulls