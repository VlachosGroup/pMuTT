
import numpy as np
from Thermochemistry import constants as c
from Thermochemistry import get_molecular_weight as mw
from Thermochemistry.models.empirical import BaseThermo

class Zacros(BaseThermo):
	"""
	Stores the information for an individual nasa specie
	Inherits from Thermochemistry.models.empirical.BaseThermo

	"""
	def __init__(self, A_st=None, inertia=None, **kwargs):
		super().__init__(**kwargs)
		self.A_st = A_st
		self.symmetrynumber = symmetrynumber
		self.inertia = inertia
		self.theta = vib_energies * c.kb('eV/K')
		self.q_vib = _np.product(np.divide(1, (1 - np.exp(-self.theta/c.T0('K')))))
		self.MW = mw(self.elements)*c.convert_unit(from_='g', to='kg')
		if self.phase == 'G':
            if self.inertia is not None:
                self.I3 = self.inertia
            else:
                self.I3 = atoms.get_moments_of_inertia()*c.convert_unit(from_='A2', to='m2')*c.convert_unit(from_='amu', to='kg')
            self.T_I = c.h('J S')**2/(8*np.pi**2*c.kb('J/K'))
		if self.phase == 'G':
			Irot = np.max(self.I3)
			if self.islinear == 0:
				self.q_rot = np.sqrt(np.pi*Irot)/self.symmetrynumber*(c.T0('K')/self.T_I)**(3./2.)
			else:
				self.q_rot = (T*Irot/self.symmetrynumber)/self.T_I
		else:
			self.q_rot = 0.
		self.q_rot = 
		if self.A_st is not None:
			self.q_trans2D = self.A_st * (2*np.pi*self.MW*c.kb('J/K')*c.T0('K'))/c.h('J s')**2
