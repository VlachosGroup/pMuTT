.. _install:

Installation
************

Installing Python
-----------------
Anaconda is the recommended method to install Python for scientific
applications. It is supported on Linux, Windows and Mac OS X.
`Download Anaconda here`_. Note that pmutt runs on Python 3.X.

Installing pmutt using pip
--------------------------
Using pip is the most straightforward way to install pmutt.

1. Open a command prompt with access to Python (if Python is installed via
   Anaconda on Windows, open the Anaconda Prompt from the start menu).

2. Install pMuTT by typing the following in the command prompt:
::

    pip install --user pmutt

The output towards the end should state "Successfully built pMuTT" if the
installation was successful.

Installing pMuTT from source
----------------------------
If you would prefer to install from source or you are interested in development,
follow the instructions below.

1. Clone the repository from GitHub (using the GitHub Desktop makes it easier).
   `See GitHub instructions on cloning repositories here`_.

2. Add the path of the ``pmutt`` folder to the ``PYTHONPATH`` environment
   variable.

Upgrading pMuTT using pip
-------------------------
To upgrade to a newer release, use the --upgrade flag:
::

    pip install --user --upgrade pmutt

Running unit tests
------------------
pMuTT has a suite of unit tests that should be run before committing any code.
To run the tests, navigate to the tests folder (pmutt/tests) via a command line
with access to Python.

Run the following command:
::

     python -m unittest

The expected output is shown below. The number of tests will not
necessarily be the same. ::

    .........................
    ----------------------------------------------------------------------
    Ran 25 tests in 0.020s

    OK

.. _`Download Anaconda here`: https://www.anaconda.com/distribution/#download-section
.. _`See GitHub instructions on cloning repositories here`: https://help.github.com/en/articles/cloning-a-repository
