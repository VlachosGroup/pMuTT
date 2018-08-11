Contributing
============

Adding a Statistical Thermodynamic Model
----------------------------------------

The overall goal is to create a framework that allows models to be added
quickly and easily. For this page, we will be using the Harmonic Thermo
model (located in ``PyMuTT.models.statmech.harmonicthermo``) as an
example. The only requirement for a class is to have the methods:
get_CpoR, get_HoRT, get_SoR, and get_GoRT. ``Ts`` must be an argument
and represents the temperature(s) in K. The initialization of the
HarmonicThermo class is as follows:

.. code:: python

   from ase import thermochemistry

   def __init__(self, vib_energies, potentialenergy=0.0):
       self.model = thermochemistry.HarmonicThermo(
           vib_energies = vib_energies,
           potentialenergy = potentialenergy)

Objects made from this class hold an
```ase.thermochemistry.harmonicthermo```_ object and uses it to
calculate the enthalpy and entropy. In this case, ``calc_HoRT`` calls
the ``get_internal_energy`` method in the ASE object.

.. code:: python

   from PyMuTT import constants as c

   def get_HoRT(self, Ts, verbose = False):
       """
       Returns the dimensionless enthalpy at a given temperature

       Parameters
           Ts - float or (N,) ndarray
               Temperature(s) in K
           verbose - bool
               Whether a table breaking down each contribution should be printed
       Returns
           float or (N,) ndarray
               Dimensionless heat capacity (H/RT) at the specified temperature
       """
       #Testing if Ts is a single value or an array
       try:
           iter(Ts)
       except TypeError:
           #Cannot be iterated over and therefore is a single value
           HoRT = self.model.get_internal_energy(temperature=Ts, verbose=verbose)/(c.kb('eV/K') * Ts)
       else:
           #Can be iterated over and therefore is an array
           HoRT = np.zeros_like(Ts)
           for i, T in enumerate(Ts):
               HoRT[i] = self.model.get_internal_energy(temperature=T, verbose=verbose)/(c.kb('eV/K') * T)
       return HoRT

However, ASE’s HarmonicThermo object does not have a method to calculate
heat capacity so we have have to implement it ourselves.

.. code:: python

   from PyMuTT.models.statmech.heat_capacity import get_CvoR_vib

   def get_CpoR(self, Ts):
       """
       Calculates the dimensionless heat capacity (Cp/R) at a given temperature.
       If you would like to use different behavior from the default, the
       thermo_model used must have the method 'get_CpoR'.

       Parameters
           Ts - float or (N,) ndarray
               Temperature(s) in K
       Returns
           float or (N,) ndarray
               Dimensionless heat capacity (Cp/R)
       """ 
       return get_CvoR_vib(vib_energies=self.model.vib_energies, Ts=Ts)

--------------

Adding an Empirical Model
-------------------------

The empirical models are a child class of the ``BaseThermo`` class
(located in ``Themochemistry.models.empirical.BaseThermo``). First, we
will discuss the ``BaseThermo`` class. Second, we will show the ``Nasa``
class (located in ``PyMuTT.models.empirical.nasa.Nasa``) that builds on
it.

BaseThermo
~~~~~~~~~~

``BaseThermo`` has several attributes that will be applicable to many
thermodynamic problems. The initialization is show

.. _``ase.thermochemistry.harmonicthermo``: https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html#ase.thermochemistry.HarmonicThermo