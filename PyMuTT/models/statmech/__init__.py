# -*- coding: utf-8 -*-
"""
PyMuTT.models.statmech
Vlachos group code for thermodynamic models.
Created on Fri Jul 7 12:40:00 2018
"""
import inspect
import numpy as np
from PyMuTT import _pass_expected_arguments

class EmptyMode:
    """Placeholder mode that returns 0 for all functions."""
    def get_q(self):
        return 1.
    def get_CvoR(self):
        return 0.
    def get_CpoR(self):
        return 0.
    def get_UoRT(self):
        return 0.
    def get_HoRT(self):
        return 0.
    def get_SoR(self):
        return 0.
    def get_AoRT(self):
        return 0.
    def get_GoRT(self):
        return 0.

class StatMech:
    """Base class for statistical mechanic models.

    Attributes
    ----------
        trans_model : `PyMuTT.models.statmech.trans` object
            Deals with translational modes
        vib_model : `PyMuTT.models.statmech.vib` object
            Deals with vibrational modes
        rot_model : `PyMuTT.models.statmech.rot` object
            Deals with rotational modes
        elec_model : `PyMuTT.models.statmech.elec` object
            Deals with electronic modes
        nucl_model : `PyMuTT.models.statmech.nucl` object
            Deals with nuclear modes
    """

    def __init__(self, trans_model=EmptyMode(), vib_model=EmptyMode(), 
                 rot_model=EmptyMode(), elec_model=EmptyMode(), 
                 nucl_model=EmptyMode(), **kwargs):

        # Translational modes
        if inspect.isclass(trans_model):
            self.trans_model = _pass_expected_arguments(trans_model, **kwargs)
        else:
            self.trans_model = trans_model
        
        # Vibrational modes
        if inspect.isclass(vib_model):
            self.vib_model = _pass_expected_arguments(vib_model, **kwargs)
        else:
            self.vib_model = vib_model

        # Rotational modes
        if inspect.isclass(rot_model):
            self.rot_model = _pass_expected_arguments(rot_model, **kwargs)
        else:
            self.rot_model = rot_model

        # Electronic modes
        if inspect.isclass(elec_model):
            self.elec_model = _pass_expected_arguments(elec_model, **kwargs)
        else:
            self.elec_model = elec_model

        # Nuclear modes
        if inspect.isclass(nucl_model):
            self.nucl_model = _pass_expected_arguments(nucl_model, **kwargs)
        else:
            self.nucl_model = nucl_model

    def get_q(self, verbose=False, **kwargs):
        """Partition function

        Parameters
        ----------
            kwargs : key-word arguments
                Parametres passed to each mode
            verbose : bool, optional
                If False, returns the product of partition functions. If True,
                returns contributions of each mode
        Returns
        -------
            q : float or tuple
                Partition function. If tuple returned, contribution to each
                mode are as follows:
                [trans, vib, rot, elec, nucl]
        """
        q = (
            _pass_expected_arguments(self.trans_model.get_q, **kwargs),
            _pass_expected_arguments(self.vib_model.get_q, **kwargs),
            _pass_expected_arguments(self.rot_model.get_q, **kwargs),
            _pass_expected_arguments(self.elec_model.get_q, **kwargs),
            _pass_expected_arguments(self.nucl_model.get_q, **kwargs))
        
        if verbose:
            return q
        else:
            return np.prod(q)

    def get_CvoR(self, verbose=False, **kwargs):
        """Dimensionless heat capacity (constant V)

        Parameters
        ----------
            kwargs : key-word arguments
                Parameters passed to each mode
            verbose : bool, optional
                If False, returns the total heat capacity. If True, returns
                contribution of each mode.
        Returns
        -------
            CvoR : float or tuple
                Dimensionless heat capacity. If tuple returned, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]
        """
        CvoR = (
            _pass_expected_arguments(self.trans_model.get_CvoR, **kwargs),
            _pass_expected_arguments(self.vib_model.get_CvoR, **kwargs),
            _pass_expected_arguments(self.rot_model.get_CvoR, **kwargs),
            _pass_expected_arguments(self.elec_model.get_CvoR, **kwargs),
            _pass_expected_arguments(self.nucl_model.get_CvoR, **kwargs))
        
        if verbose:
            return CvoR
        else:
            return np.sum(CvoR)


    def get_CpoR(self, verbose=False, **kwargs):
        """Dimensionless heat capacity (constant P)

        Parameters
        ----------
            kwargs : key-word arguments
                Parameters passed to each mode
            verbose : bool, optional
                If False, returns the total heat capacity. If True, returns
                contribution of each mode.
        Returns
        -------
            CpoR : float or tuple
                Dimensionless heat capacity. If tuple returned, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]
        """
        CpoR = (
            _pass_expected_arguments(self.trans_model.get_CpoR, **kwargs),
            _pass_expected_arguments(self.vib_model.get_CpoR, **kwargs),
            _pass_expected_arguments(self.rot_model.get_CpoR, **kwargs),
            _pass_expected_arguments(self.elec_model.get_CpoR, **kwargs),
            _pass_expected_arguments(self.nucl_model.get_CpoR, **kwargs))
        
        if verbose:
            return CpoR
        else:
            return np.sum(CpoR)

    def get_UoRT(self, verbose=False, **kwargs):
        """Dimensionless internal energy

        Parameters
        ----------
            kwargs : key-word arguments
                Parameters passed to each mode
            verbose : bool, optional
                If False, returns the total internal energy. If True, returns
                contribution of each mode.
        Returns
        -------
            UoRT : float or tuple
                Dimensionless internal energy. If tuple returned, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]
        """

        UoRT = (
            _pass_expected_arguments(self.trans_model.get_UoRT, **kwargs),
            _pass_expected_arguments(self.vib_model.get_UoRT, **kwargs),
            _pass_expected_arguments(self.rot_model.get_UoRT, **kwargs),
            _pass_expected_arguments(self.elec_model.get_UoRT, **kwargs),
            _pass_expected_arguments(self.nucl_model.get_UoRT, **kwargs))
        
        if verbose:
            return UoRT
        else:
            return np.sum(UoRT)

    def get_HoRT(self, verbose=False, **kwargs):
        """Dimensionless enthalpy

        Parameters
        ----------
            kwargs : key-word arguments
                Parameters passed to each mode
            verbose : bool, optional
                If False, returns the total enthalpy. If True, returns
                contribution of each mode.
        Returns
        -------
            HoRT : float or tuple
                Dimensionless enthalpy. If tuple returned, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]
        """

        HoRT = (
            _pass_expected_arguments(self.trans_model.get_HoRT, **kwargs),
            _pass_expected_arguments(self.vib_model.get_HoRT, **kwargs),
            _pass_expected_arguments(self.rot_model.get_HoRT, **kwargs),
            _pass_expected_arguments(self.elec_model.get_HoRT, **kwargs),
            _pass_expected_arguments(self.nucl_model.get_HoRT, **kwargs))
        
        if verbose:
            return HoRT
        else:
            return np.sum(HoRT)

    def get_SoR(self, verbose=False, **kwargs):
        """Dimensionless entropy

        Parameters
        ----------
            kwargs : key-word arguments
                Parameters passed to each mode
            verbose : bool, optional
                If False, returns the total entropy. If True, returns
                contribution of each mode.
        Returns
        -------
            SoR : float or tuple
                Dimensionless entropy. If tuple returned, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]
        """

        SoR = (
            _pass_expected_arguments(self.trans_model.get_SoR, **kwargs),
            _pass_expected_arguments(self.vib_model.get_SoR, **kwargs),
            _pass_expected_arguments(self.rot_model.get_SoR, **kwargs),
            _pass_expected_arguments(self.elec_model.get_SoR, **kwargs),
            _pass_expected_arguments(self.nucl_model.get_SoR, **kwargs))
        
        if verbose:
            return SoR
        else:
            return np.sum(SoR)

    def get_AoRT(self, verbose=False, **kwargs):
        """Dimensionless Helmholtz energy

        Parameters
        ----------
            kwargs : key-word arguments
                Parameters passed to each mode
            verbose : bool, optional
                If False, returns the total Helmholtz energy. If True, returns
                contribution of each mode.
        Returns
        -------
            AoRT : float or tuple
                Dimensionless Helmoltz energy. If tuple returned, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]
        """

        AoRT = (
            _pass_expected_arguments(self.trans_model.get_AoRT, **kwargs),
            _pass_expected_arguments(self.vib_model.get_AoRT, **kwargs),
            _pass_expected_arguments(self.rot_model.get_AoRT, **kwargs),
            _pass_expected_arguments(self.elec_model.get_AoRT, **kwargs),
            _pass_expected_arguments(self.nucl_model.get_AoRT, **kwargs))
        
        if verbose:
            return AoRT
        else:
            return np.sum(AoRT)

    def get_GoRT(self, verbose=False, **kwargs):
        """Dimensionless Gibbs energy

        Parameters
        ----------
            kwargs : key-word arguments
                Parameters passed to each mode
            verbose : bool, optional
                If False, returns the total Gibbs energy. If True, returns
                contribution of each mode.
        Returns
        -------
            GoRT : float or tuple
                Dimensionless Gibbs energy. If tuple returned, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]
        """

        GoRT = (
            _pass_expected_arguments(self.trans_model.get_GoRT, **kwargs),
            _pass_expected_arguments(self.vib_model.get_GoRT, **kwargs),
            _pass_expected_arguments(self.rot_model.get_GoRT, **kwargs),
            _pass_expected_arguments(self.elec_model.get_GoRT, **kwargs),
            _pass_expected_arguments(self.nucl_model.get_GoRT, **kwargs))
        
        if verbose:
            return GoRT
        else:
            return np.sum(GoRT)

