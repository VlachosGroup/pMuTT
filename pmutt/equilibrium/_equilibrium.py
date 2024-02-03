from scipy.optimize import minimize
import numpy as np
import sys
from itertools import repeat
from pmutt.io.thermdat import read_thermdat
from pmutt import pmutt_list_to_dict
from pmutt import constants as c
from collections import namedtuple


class Equilibrium():
    """Reaction thermodynamic equilibrium.

    Attributes
    ----------
        model : list or dictionary of :class:`~pmutt.empirical` or :class:`~pmutt.statmech` objects
            Thermodynamic object to compute Gibbs free energy for all
            species in 'network'
        network : dictionary object
            List of species to consider in equilibrium calculation including
            initial moles of each species. All species names must match and
            be contained in 'model' data.
    """

    def __init__(self,
                 model,
                 network):
        self.model = model
        self.network = network
        if type(self.model) is list:
            self.model = pmutt_list_to_dict(self.model)
        elif type(self.model) is dict:
            pass
        else:
            sys.exit('model must be list or dict')

        # Build molecule-element configuration matrix from species
        self.elements = []
        self.species = list(self.network.keys())
        feed = np.array(list(network.values()))
        self.mol_elem = np.zeros([len(self.species), 2])
        for i, x in enumerate(self.species):
            ele = self.model[x].elements
            # Read elements in each species
            for y in ele:
                # Check if the current element is in the list
                try:
                    self.elements.index(y)
                except ValueError:
                    # If not add the element to the list
                    self.elements.append(y)
                    # Check if the molecule-element matrix has sufficient
                    # colums to accomodate the aditional element and add a
                    # column if necessary
                    if len(self.elements) > np.size(self.mol_elem, 1):
                        self.mol_elem = np.append(self.mol_elem,
                                                  np.zeros([len(self.species),
                                                            1]), 1)
                # Enter the quantity of the current element of the current
                # molecule into the molecule-element matrix
                self.mol_elem[i, self.elements.index(y)] =\
                    self.model[x].elements[y]

        # Elimnate zero columns
        self.elements = list(np.array(self.elements)
                             [sum(self.mol_elem, 0) > 0])
        self.mol_elem = self.mol_elem[:, sum(self.mol_elem, 0) > 0]
        # Determine the moles of each element in the feed
        self.ele_feed = feed.dot(self.mol_elem)
        self.species_mw = self.mol_elem.dot([c.atomic_weight[x]
                                             for x in self.elements])

    # Objective (Cost) Function: Summ of Gibb's Free Energies

    def _objective(self, x, *args):
        s = 0.0
        nT = sum(x)
        g = np.array(args[0])
        p = args[1]
        for i in range(len(x)):
            s += x[i]*(g[i] + np.log(x[i]*p/nT))
        # Return sum of the Gibb's free energies
        return s

    # Objective (Cost) Function Jacobian: Summ of Gibb's Free Energies

    def _objective_jac(self, x, *args):
        s = np.zeros_like(x)
        nT = sum(x)
        g = np.array(args[0])
        p = args[1]
        for i in range(len(x)):
            s[i] = g[i] + np.log(x[i]*p/nT)
        # Return sum of the Gibb's free energies
        return s

    # Elemental Balance Equality Constraint. The returned value
    # must be = 0

    def _constraints1_eq(self, x):
        s = x.dot(self.mol_elem) - self.ele_feed  # Elemental balances
        # Return slack variable for each element
        return s

    # Elemental Balance Equality Constraint Jacobian.

    def _constraints1_eq_jac(self, x):
        # Return jacobian
        return self.mol_elem.T

    def get_net_comp(self, T, P):
        """Returns the equilibrium composition of the specified molecule
        network.

        Parameters
        ----------
            T : float
                Temperature in K
            P : float
                Pressure in atm
        Returns
        -------
            res : equilibrium._equilibrium.res

            Important attributes are
            .species : list of strings
                list of species in network
            .moles : `numpy.ndarray`_
                Equilibrium moles of each species in the network
            .mole_frac : `numpy.ndarray`_
                Equilibrium mole fraction of each species in the network
            .P : float
                Pressure (atm) used in equilibrium calculation
            .T : float
                Temperature (K) used in equilibrium calculation

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        self.P = P
        self.T = T
        # Model initialization parameters
        # Starting mole guesses = 1
        self.guess = list(repeat(1.0, len(self.species)))
        # Mole value bounds. Lower bound near zero
        b = [1e-20, sum(self.ele_feed)]
        # Upper bound is the total moles of elements
        self.bounds = list(repeat(b, len(self.species)))

        # Designate each constraint as an equality or inequality constraint
        self.con = {'type': 'eq', 'fun': self._constraints1_eq,
                    'jac': self._constraints1_eq_jac}

        self.maxiter = 5000  # Maximum iterations for minimize solver
        # Initialize model results lists

        # t0 = time.time()
        # Calculate species Gibb's at current temperature
        self.gibbs = []
        for x in self.species:
            self.gibbs.append(self.model[x].get_GoRT(T=T))

        # Run solver once and collect data
        sol = minimize(self._objective, self.guess,
                       args=(self.gibbs, self.P*1.01325),
                       jac=self._objective_jac,
                       method='SLSQP',
                       options={'ftol': 1e-14, 'maxiter': self.maxiter},
                       bounds=self.bounds,
                       constraints=self.con)

        res = namedtuple("res", ["species", "moles", "mole_frac", "P", "T"])

        return res(self.species, sol.x, sol.x/np.sum(sol.x), self.P, self.T)

    @classmethod
    def from_thermdat(cls,
                      thermdat,
                      network):
        """Reaction thermodynamic equilibrium.

        Attributes
        ----------
            thermdat : string, filepath to thermdat file
                File in thermdat format containing NASA polynomials for
                species in network to compute Gibbs free energy
            network : dictionary object
                List of species to consider in equilibrium calculation
                including initial moles of each species. All species names
                must match and be contained in 'model' data.
        """
        model = read_thermdat(thermdat, "dict")
        return cls(model=model, network=network)
