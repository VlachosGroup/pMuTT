.. _install:

Installation
************

Installing Python
-----------------
Anaconda is the recommended method to install Python for scientific
applications. It is supported on Linux, Windows and Mac OS X.
`Download Anaconda here`_. Note that pMuTT runs on Python 3.X.

Installing pMuTT using pip
--------------------------
Using pip is the most straightforward way to install pMuTT.

1. Open a command prompt with access to Python (if Python is installed via
   Anaconda on Windows, open the Anaconda Prompt from the start menu).

2. Install pMuTT by typing the following in the command prompt:
::

    pip install pmutt

The output towards the end should state "Successfully built pMuTT" if the
installation was successful. 

Receiving installation errors? `Post the error to our Issues page`_.

Common installation errors
^^^^^^^^^^^^^^^^^^^^^^^^^^

PyYAML Uninstallation Error
+++++++++++++++++++++++++++

**Error**::

   ERROR: Cannot uninstall 'PyYAML'. It is a distutils installed project and
   thus we cannot accurately determine which files belong to it which would
   lead to only a partial uninstall.

**Solution**

Append ``--ignore-installed PyYAML`` to pip command.::

   pip install pmutt --ignore-installed PyYAML

See issue regarding the `PyYAML error here`_.

Pip Permission Error
++++++++++++++++++++

**Error**::

   Could not install packages due to an EnvironmentError: [Errno 13] Permission
   denied: '/usr/local/bin/pmutt'

**Solution**

Append ``--user`` to pip command.::

   pip install pmutt --user

See explanation why this `permission error occurs here`_.

Installing pMuTT from source
----------------------------
If you would prefer to install from source or you are interested in development,
follow the instructions below.
::

    pip install git+https://github.com/VlachosGroup/pMuTT.git

Installing the developer branch
-------------------------------
pMuTT has a release roughly once a month. Changes that will be in the next
release will be located in the Developer branch but may have more bugs than
the master branch. You can install using the following:
::


    pip install --upgrade git+https://github.com/VlachosGroup/pMuTT.git@development

Upgrading pMuTT using pip
-------------------------
To upgrade to a newer release, use the --upgrade flag:
::

    pip install --upgrade pmutt

Running unit tests
------------------
pMuTT has a suite of unit tests that should be run before committing any code.
To run the tests, run the following commands in a Python terminal.
::

     import pmutt
     pmutt.run_tests()

The expected output is shown below. The number of tests will not
necessarily be the same. ::

    .........................
    ----------------------------------------------------------------------
    Ran 25 tests in 0.020s

    OK

.. _`Download Anaconda here`: https://www.anaconda.com/distribution/#download-section
.. _`See GitHub instructions on cloning repositories here`: https://help.github.com/en/articles/cloning-a-repository
.. _`PyYAML error here`: https://github.com/VlachosGroup/pMuTT/issues/132
.. _`permission error occurs here`: https://pip.pypa.io/en/stable/user_guide/#user-installs
.. _`Post the error to our Issues page`: https://github.com/VlachosGroup/pMuTT/issues