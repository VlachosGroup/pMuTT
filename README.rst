Python Multiscale Thermochemistry Toolbox (PyMuTT)
==================================================

The **Py**\ thon **Mu**\ ltiscale **T**\ hermochemistry **T**\ oolbox
(PyMuTT) is a Python library for Thermochemistry developed by the
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

-  Gerhard Wittreich, P.E. (wittregr@udel.edu)
-  Jonathan Lym (jlym@udel.edu)

Dependencies
------------

-  Python3
-  `Atomic Simulation Environment`_: Used for I/O operations and to
   calculate thermodynamic properties
-  `Numpy`_: Used for vector and matrix operations
-  `Pandas`_: Used to import data from Excel files
-  `SciPy`_: Used for fitting heat capacities.
-  `Matplotlib`_: Used for plotting thermodynamic data

Getting Started
---------------

1. Install the dependencies
2. `Clone the repository`_ to your local machine
3. Add the parent folder to PYTHONPATH. `Windows instructions`_
4. Run the tests by navigating to the `tests directory`_ in a
   command-line interface and inputting the following command:

::

   python -m unittest

The expected output is shown below. The number of tests will not
necessarily be the same.

::

   .........................
   ----------------------------------------------------------------------
   Ran 25 tests in 0.020s

   OK

5. Look at `examples using the code`_

License
-------

This project is licensed under the MIT License - see the `LICENSE.md`_
file for details.

Contributing
------------

If you have a suggestion, please post to our `Issues page`_ with the ``enhancement`` tag. Similarly, if you 
encounter a bug, please post to our `Issues page`_ with the ``bug`` tag. Finally, if you would like to add 
to the body of code, please check our documentation to make sure the new code is consistent with the relevant 
page and submit a `pull request`_.

Questions
---------

If you are having issues, please post to our `Issues page`_ with the ``help wanted`` or ``question`` tag. We 
will do our best to assist.

.. _`Clone the repository`: https://help.github.com/articles/cloning-a-repository/
.. _`Windows instructions`: https://docs.python.org/3/using/windows.html#excursus-setting-environment-variables
.. _`documentation page`: https://vlachosgroup.github.io/PyMuTT/
.. _Atomic Simulation Environment: https://wiki.fysik.dtu.dk/ase/
.. _Numpy: http://www.numpy.org/
.. _Pandas: https://pandas.pydata.org/
.. _SciPy: https://www.scipy.org/
.. _Matplotlib: https://matplotlib.org/
.. _tests directory: https://github.com/VlachosGroup/PyMuTT/tree/master/tests
.. _LICENSE.md: https://github.com/VlachosGroup/PyMuTT/blob/master/LICENSE.md
.. _`examples using the code`: https://github.com/VlachosGroup/PyMuTT/tree/master/examples
.. _`Issues page`: https://github.com/VlachosGroup/PyMuTT/issues
.. _`pull request`: https://github.com/VlachosGroup/PyMuTT/pulls