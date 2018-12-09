# -*- coding: utf-8 -*-
"""
pMuTT.test_pMuTT_model_reaction_bep
Tests for pMuTT module
"""
import unittest
import numpy as np
from ase.build import molecule
from pMuTT import constants as c
from pMuTT.models.reaction import Reaction
from pMuTT.models.reaction.bep import BEP
from pMuTT.models.statmech import StatMech, presets

class TestBEP(unittest.TestCase):
    def setUp(self):
        self.T = c.T0('K')
        # Factor to convert potential energy to dimensionless number
        dim_factor = c.R('eV/K')*self.T
        self.m = 10. # BEP Slope 
        self.c = 20. # BEP Intercept
        species = {
            'H2': StatMech(name='H2', potentialenergy=2.*dim_factor, 
                           **presets['electronic']),
            'O2': StatMech(name='O2', potentialenergy=4.*dim_factor, 
                           **presets['electronic']),
            'H2O': StatMech(name='H2O', potentialenergy=3.*dim_factor, 
                            **presets['electronic']),
        }
        self.rxn_delta = Reaction.from_string(
                reaction_str='H2 + 0.5O2 = H2O', 
                species=species,
                descriptor='delta', 
                slope=self.m, 
                intercept=self.c, 
                bep=BEP)
        self.bep_delta = self.rxn_delta.bep

        rxn_rev_delta = Reaction.from_string(
                reaction_str='H2 + 0.5O2 = H2O', 
                species=species,
                descriptor='rev_delta', 
                slope=self.m, 
                intercept=self.c, 
                bep=BEP)
        self.bep_rev_delta = rxn_rev_delta.bep

        rxn_reactants = Reaction.from_string(
                reaction_str='H2 + 0.5O2 = H2O', 
                species=species,
                descriptor='reactants', 
                slope=self.m, 
                intercept=self.c, 
                bep=BEP)
        self.bep_reactants = rxn_reactants.bep

        rxn_products = Reaction.from_string(
                reaction_str='H2 + 0.5O2 = H2O', 
                species=species,
                descriptor='products', 
                slope=self.m, 
                intercept=self.c, 
                bep=BEP)
        self.bep_products = rxn_products.bep

    def test_get_EoRT(self):
        delta_H = self.rxn_delta.get_delta_HoRT(T=self.T)
        exp_EoRT_delta = self.m*delta_H + self.c
        exp_EoRT_delta_rev = (self.m-1.)*delta_H + self.c

        self.assertAlmostEqual(self.bep_delta.get_EoRT_act(T=self.T, rev=False), 
                               exp_EoRT_delta)
        self.assertAlmostEqual(self.bep_delta.get_EoRT_act(T=self.T, rev=True), 
                               exp_EoRT_delta_rev)

        rev_delta_H = self.rxn_delta.get_delta_HoRT(T=self.T, rev=True)
        exp_EoRT_delta = (self.m-1.)*rev_delta_H + self.c
        exp_EoRT_delta_rev = self.m*rev_delta_H + self.c

        self.assertAlmostEqual(self.bep_rev_delta.get_EoRT_act(T=self.T, 
                                                               rev=False), 
                               exp_EoRT_delta)
        self.assertAlmostEqual(self.bep_rev_delta.get_EoRT_act(T=self.T, 
                                                               rev=True), 
                               exp_EoRT_delta_rev)

        reactants_H = self.rxn_delta.get_HoRT_state(state='reactants', T=self.T)
        exp_EoRT_delta = self.m*reactants_H + self.c
        exp_EoRT_delta_rev = (self.m-1.)*reactants_H + self.c

        self.assertAlmostEqual(self.bep_reactants.get_EoRT_act(T=self.T, 
                                                               rev=False), 
                               exp_EoRT_delta)
        self.assertAlmostEqual(self.bep_reactants.get_EoRT_act(T=self.T, 
                                                               rev=True), 
                               exp_EoRT_delta_rev)
 
 
        products_H = self.rxn_delta.get_HoRT_state(state='products', T=self.T)        
        exp_EoRT_delta = self.m*products_H + self.c
        exp_EoRT_delta_rev = (self.m-1.)*products_H + self.c

        self.assertAlmostEqual(self.bep_products.get_EoRT_act(T=self.T, 
                                                               rev=False), 
                               exp_EoRT_delta)
        self.assertAlmostEqual(self.bep_products.get_EoRT_act(T=self.T, 
                                                               rev=True), 
                               exp_EoRT_delta_rev)


if __name__ == '__main__':
    unittest.main()