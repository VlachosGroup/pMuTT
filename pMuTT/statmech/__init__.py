# -*- coding: utf-8 -*-
"""
pMuTT.statmech
Vlachos group code for thermodynamic models.
Created on Fri Jul 7 12:40:00 2018
"""
import inspect
import warnings
import numpy as np
from pMuTT import _pass_expected_arguments
from pMuTT import constants as c
from pMuTT.statmech import trans, vib, elec, rot
from pMuTT.io_ import jsonio as json_pMuTT


class EmptyMode:
    """Placeholder mode that returns 1 for partition function and
    0 for all functions other thermodynamic properties."""
    def __init__(self):
        pass

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

    def get_FoRT(self):
        return 0.

    def get_GoRT(self):
        return 0.

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        return {'class': str(self.__class__)}

    @classmethod
    def from_dict(cls, json_obj):
        """Recreate an object from the JSON representation.

        Parameters
        ----------
            json_obj : dict
                JSON representation
        Returns
        -------
            EmptyMode : EmptyMode object
        """
        return cls()


class StatMech:
    """Base class for statistical mechanic models.

    Attributes
    ----------
        name : str, optional
            Name of the specie. Default is None
        trans_model : :ref:`pMuTT.statmech.trans <trans>` object, optional
            Deals with translational modes. Default is
            :class:`~pMuTT.statmech.EmptyMode`
        vib_model : :ref:`pMuTT.statmech.vib <vib>` object, optional
            Deals with vibrational modes. Default is
            :class:`~pMuTT.statmech.EmptyMode`
        rot_model : :ref:`pMuTT.statmech.rot <rot>` object, optional
            Deals with rotational modes. Default is
            :class:`~pMuTT.statmech.EmptyMode`
        elec_model : :ref:`pMuTT.statmech.elec <elec>` object, optional
            Deals with electronic modes. Default is
            :class:`~pMuTT.statmech.EmptyMode`
        nucl_model : :ref:`pMuTT.statmech.nucl <nucl>` object
            Deals with nuclear modes. Default is
            :class:`~pMuTT.statmech.EmptyMode`
        notes : str, optional
            Any additional details you would like to include such as
            computational set up. Default is None
    """

    def __init__(self, name=None, trans_model=EmptyMode(),
                 vib_model=EmptyMode(), rot_model=EmptyMode(),
                 elec_model=EmptyMode(), nucl_model=EmptyMode(), notes=None,
                 **kwargs):
        self.name = name

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

        self.notes = notes

    def get_quantity(self, method_name, raise_error=True, raise_warning=True, 
                     operation='sum', verbose=False, **kwargs):
        """Generic method to get any quantity from modes.

        Parameters
        ----------
            method_name : str
                Name of method to use to calculate quantity. Calculates any
                quantity as long as the relevant objects have the same method
                name
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            operation : str, optional
                Operation to apply when combining the modes. Supported options
                include:
                
                - sum (Default)
                - prod
            verbose : bool, optional
                If False, returns the total Gibbs energy. If True, returns
                contribution of each mode.
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            quantity : float or (N,) `numpy.ndarray`_
                Dimensionless Gibbs energy. If verbose is True, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        # Get the default value
        operation = operation.lower()
        if operation == 'sum':
            default_value = 0.
        elif operation == 'prod':
            default_value = 1.
        else:
            raise ValueError('Operation: {} not supported'.format(operation))

        # Calculate the quantity for each mode
        quantity = np.array([
            _get_mode_quantity(mode=self.trans_model, method_name=method_name,
                               raise_error=raise_error,
                               raise_warning=raise_warning, 
                               default_value=default_value,
                               **kwargs),
            _get_mode_quantity(mode=self.vib_model, method_name=method_name,
                               raise_error=raise_error,
                               raise_warning=raise_warning, 
                               default_value=default_value,
                               **kwargs),
            _get_mode_quantity(mode=self.rot_model, method_name=method_name,
                               raise_error=raise_error,
                               raise_warning=raise_warning, 
                               default_value=default_value,
                               **kwargs),
            _get_mode_quantity(mode=self.elec_model, method_name=method_name,
                               raise_error=raise_error,
                               raise_warning=raise_warning, 
                               default_value=default_value,
                               **kwargs),
            _get_mode_quantity(mode=self.nucl_model, method_name=method_name,
                               raise_error=raise_error,
                               raise_warning=raise_warning, 
                               default_value=default_value,
                               **kwargs)])
        # Return value
        if verbose:
            return quantity
        elif operation == 'sum':
            return np.sum(quantity)
        elif operation == 'prod':
            return np.prod(quantity)
        else:
            raise ValueError('Operation: {} not supported'.format(operation))

    def get_q(self, verbose=False, raise_error=True, raise_warning=True,
              **kwargs):
        """Partition function

        Parameters
        ----------
            verbose : bool, optional
                If False, returns the product of partition functions. If True,
                returns contributions of each mode
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            q : float or (N,) `numpy.ndarray`_
                Partition function. If verbose is True, contribution to each
                mode are as follows:
                [trans, vib, rot, elec, nucl]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_q', raise_error=raise_error,
                                 raise_warning=raise_warning, operation='prod',
                                 verbose=verbose, **kwargs)

    def get_CvoR(self, verbose=False, raise_error=True, raise_warning=True,
                 **kwargs):
        """Dimensionless heat capacity (constant V)

        Parameters
        ----------
            verbose : bool, optional
                If False, returns the total heat capacity. If True, returns
                contribution of each mode.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            CvoR : float or (N,) `numpy.ndarray`_
                Dimensionless heat capacity. If verbose is True, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_CvoR', operation='sum',
                                 raise_error=raise_error, 
                                 raise_warning=raise_warning,
                                 verbose=verbose, **kwargs)

    def get_Cv(self, units, verbose=False, raise_error=True, raise_warning=True,
               **kwargs):
        """Calculate the heat capacity (constant V)

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            verbose : bool, optional
                If False, returns the total heat capacity. If True, returns
                contribution of each mode.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            Cv : float or (N,) `numpy.ndarray`_
                Heat capacity

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_CvoR(verbose=verbose, raise_error=raise_error, 
                             raise_warning=raise_warning, **kwargs)*c.R(units)

    def get_CpoR(self, verbose=False, raise_error=True, raise_warning=True,
                 **kwargs):
        """Dimensionless heat capacity (constant P)

        Parameters
        ----------
            verbose : bool, optional
                If False, returns the total heat capacity. If True, returns
                contribution of each mode.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            CpoR : float or (N,) `numpy.ndarray`_
                Dimensionless heat capacity. If verbose is True, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_CpoR', operation='sum',
                                 raise_error=raise_error, 
                                 raise_warning=raise_warning,
                                 verbose=verbose, **kwargs)

    def get_Cp(self, units, verbose=False, raise_error=True, raise_warning=True,
               **kwargs):
        """Calculate the heat capacity (constant P)

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            verbose : bool, optional
                If False, returns the total heat capacity. If True, returns
                contribution of each mode.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            Cp : float or (N,) `numpy.ndarray`_
                Heat capacity

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_CpoR(verbose=verbose, raise_error=raise_error,
                             raise_warning=raise_warning, **kwargs)*c.R(units)

    def get_EoRT(self, T=c.T0('K'), include_ZPE=False, 
                 raise_error=True, raise_warning=True, **kwargs):
        """Dimensionless electronic energy

        Parameters
        ----------
            T : float, optional
                Temperature in K. Default is 298.15 K
            include_ZPE : bool, optional
                If True, includes the zero point energy. Default is False
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to electronic mode
        Returns
        -------
            EoRT : float
                Dimensionless electronic energy.
        """
        kwargs['T'] = T
        EoRT = _get_mode_quantity(mode=self.elec_model, 
                                  method_name='get_UoRT',
                                  raise_error=raise_error, 
                                  raise_warning=raise_warning,
                                  default_value=0., **kwargs)
        if include_ZPE:
            EoRT += _get_mode_quantity(mode=self.vib_model, 
                                       method_name='get_ZPE',
                                       raise_error=raise_error, 
                                       raise_warning=raise_warning,
                                       default_value=0., **kwargs)/c.R('eV/K')/T
        return EoRT            

    def get_E(self, units, T=c.T0('K'), raise_error=True, raise_warning=True,
              include_ZPE=False, **kwargs):
        """Calculate the electronic energy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float, optional
                Temperature in K. Default is 298.15 K
            include_ZPE : bool, optional
                If True, includes the zero point energy. Default is False
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            E : float
                Electronic energy
        """
        return self.get_EoRT(T=T, raise_error=raise_error,
                             raise_warning=raise_warning, 
                             include_ZPE=include_ZPE, **kwargs) \
               *T*c.R('{}/K'.format(units))

    def get_UoRT(self, verbose=False, raise_error=True, raise_warning=True,
                 **kwargs):
        """Dimensionless internal energy

        Parameters
        ----------
            verbose : bool, optional
                If False, returns the total internal energy. If True, returns
                contribution of each mode.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            UoRT : float or (N,) `numpy.ndarray`_
                Dimensionless internal energy. If verbose is True, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_UoRT', operation='sum',
                                 raise_error=raise_error, 
                                 raise_warning=raise_warning,
                                 verbose=verbose, **kwargs)

    def get_U(self, units, T=c.T0('K'), verbose=False, raise_error=True,
              raise_warning=True, **kwargs):
        """Calculate the internal energy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float, optional
                Temperature in K. Default is 298.15 K
            verbose : bool, optional
                If False, returns the internal energy. If True, returns
                contribution of each mode.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            U : float or (N,) `numpy.ndarray`_
                Internal energy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_UoRT(verbose=verbose, T=T, raise_error=raise_error,
                             raise_warning=raise_warning, **kwargs) \
               *T*c.R('{}/K'.format(units))

    def get_HoRT(self, verbose=False, raise_error=True, raise_warning=True,
                 **kwargs):
        """Dimensionless enthalpy

        Parameters
        ----------
            verbose : bool, optional
                If False, returns the total enthalpy. If True, returns
                contribution of each mode.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            HoRT : float or (N,) `numpy.ndarray`_
                Dimensionless enthalpy. If verbose is True, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_HoRT', operation='sum',
                                 raise_error=raise_error, 
                                 raise_warning=raise_warning,
                                 verbose=verbose, **kwargs)

    def get_H(self, units, T=c.T0('K'), raise_error=True, raise_warning=True,
              verbose=False, **kwargs):
        """Calculate the enthalpy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float, optional
                Temperature in K. Default is 298.15 K
            verbose : bool, optional
                If False, returns the enthalpy. If True, returns
                contribution of each mode.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            H : float or (N,) `numpy.ndarray`_
                Enthalpy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_HoRT(verbose=verbose, raise_error=raise_error,
                             raise_warning=raise_warning, T=T, **kwargs) \
               *T*c.R('{}/K'.format(units))

    def get_SoR(self, verbose=False, raise_error=True, raise_warning=True,
                **kwargs):
        """Dimensionless entropy

        Parameters
        ----------
            verbose : bool, optional
                If False, returns the total entropy. If True, returns
                contribution of each mode.
            kwargs : key-word arguments
                Parameters passed to each mode
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
        Returns
        -------
            SoR : float or (N,) `numpy.ndarray`_
                Dimensionless entropy. If verbose is True, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_SoR', operation='sum',
                                 raise_error=raise_error, 
                                 raise_warning=raise_warning,
                                 verbose=verbose, **kwargs)


    def get_S(self, units, verbose=False, raise_error=True, 
              raise_warning=True, **kwargs):
        """Calculate the entropy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units.
            verbose : bool, optional
                If False, returns the entropy. If True, returns
                contribution of each mode.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            S : float or (N,) `numpy.ndarray`_
                Entropy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_SoR(verbose=verbose, raise_error=raise_error,
                            raise_warning=raise_warning, **kwargs)*c.R(units)

    def get_FoRT(self, verbose=False, raise_error=True, raise_warning=True,
                 **kwargs):
        """Dimensionless Helmholtz energy

        Parameters
        ----------
            verbose : bool, optional
                If False, returns the total Helmholtz energy. If True, returns
                contribution of each mode.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            FoRT : float or (N,) `numpy.ndarray`_
                Dimensionless Helmoltz energy. If verbose is True, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_FoRT', operation='sum',
                                 raise_error=raise_error, 
                                 raise_warning=raise_warning,
                                 verbose=verbose, **kwargs)


    def get_F(self, units, T=c.T0('K'), verbose=False, raise_error=True,
              raise_warning=True, **kwargs):
        """Calculate the Helmholtz energy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float, optional
                Temperature in K. Default is 298.15 K
            verbose : bool, optional
                If False, returns the Helmholtz energy. If True, returns
                contribution of each mode.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            F : float or (N,) `numpy.ndarray`_
                Helmholtz energy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_FoRT(verbose=verbose, T=T, raise_error=raise_error,
                             raise_warning=raise_warning, **kwargs) \
               *T*c.R('{}/K'.format(units))

    def get_GoRT(self, verbose=False, raise_error=True, raise_warning=True,
                 **kwargs):
        """Dimensionless Gibbs energy

        Parameters
        ----------
            verbose : bool, optional
                If False, returns the total Gibbs energy. If True, returns
                contribution of each mode.
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            GoRT : float or (N,) `numpy.ndarray`_
                Dimensionless Gibbs energy. If verbose is True, contribution
                to each mode are as follows:
                [trans, vib, rot, elec, nucl]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_GoRT', operation='sum',
                                 raise_error=raise_error, 
                                 raise_warning=raise_warning,
                                 verbose=verbose, **kwargs)

    def get_G(self, units, T=c.T0('K'), verbose=False, raise_error=True,
              raise_warning=True, **kwargs):
        """Calculate the Gibbs energy

        Parameters
        ----------
            verbose : bool, optional
                If False, returns the Gibbs energy. If True, returns
                contribution of each mode.
            units : str
                Units as string. See :func:`~pMuTT.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float, optional
                Temperature in K. Default is 298.15 K
            raise_error : bool, optional
                If True, raises an error if any of the modes do not have the 
                quantity of interest. Default is True
            raise_warning : bool, optional
                Only relevant if raise_error is False. Raises a warning if any
                of the modes do not have the quantity of interest. Default is
                True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            G : float or (N,) `numpy.ndarray`_
                Gibbs energy

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
        """
        return self.get_GoRT(verbose=verbose, raise_error=raise_error,
                             raise_warning=raise_warning, T=T, **kwargs) \
               *T*c.R('{}/K'.format(units))

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        return {'class': str(self.__class__),
                'name': self.name,
                'trans_model': self.trans_model.to_dict(),
                'vib_model': self.vib_model.to_dict(),
                'rot_model': self.rot_model.to_dict(),
                'elec_model': self.elec_model.to_dict(),
                'nucl_model': self.nucl_model.to_dict(),
                'notes': self.notes}

    @classmethod
    def from_dict(cls, json_obj):
        """Recreate an object from the JSON representation.

        Parameters
        ----------
            json_obj : dict
                JSON representation
        Returns
        -------
            StatMech : StatMech object
        """
        json_obj = json_pMuTT.remove_class(json_obj)
        name = json_obj['name']
        trans_model = json_pMuTT.json_to_pMuTT(json_obj['trans_model'])
        vib_model = json_pMuTT.json_to_pMuTT(json_obj['vib_model'])
        rot_model = json_pMuTT.json_to_pMuTT(json_obj['rot_model'])
        elec_model = json_pMuTT.json_to_pMuTT(json_obj['elec_model'])
        nucl_model = json_pMuTT.json_to_pMuTT(json_obj['nucl_model'])
        notes = json_obj['notes']

        return cls(name=name,
                   trans_model=trans_model,
                   vib_model=vib_model,
                   rot_model=rot_model,
                   elec_model=elec_model,
                   nucl_model=nucl_model,
                   notes=notes)


presets = {
    'idealgas': {
        'statmech_model': StatMech,
        'trans_model': trans.IdealTrans,
        'n_degrees': 3,
        'vib_model': vib.HarmonicVib,
        'elec_model': elec.IdealElec,
        'rot_model': rot.RigidRotor,
        'required': ('molecular_weight', 'vib_wavenumbers', 'potentialenergy',
                     'spin', 'geometry', 'rot_temperatures', 'symmetrynumber'),
        'optional': ('atoms')
        },
    'harmonic': {
        'statmech_model': StatMech,
        'vib_model': vib.HarmonicVib,
        'elec_model': elec.IdealElec,
        'required': ('vib_wavenumbers', 'potentialenergy', 'spin'),
    },
    'electronic': {
        'statmech_model': StatMech,
        'elec_model': elec.IdealElec,
        'required': ('potentialenergy', 'spin'),
    },
    'placeholder': {
        'statmech_model': StatMech,
        'trans_model': EmptyMode,
        'vib_model': EmptyMode,
        'elec_model': EmptyMode,
        'rot_model': EmptyMode,
        'nucl_model': EmptyMode,
        'required': (),
    }
}
"""Commonly used models. The 'required' and 'optional' fields indicate
parameters that still need to be passed."""

def _get_mode_quantity(mode, method_name, raise_error=True, raise_warning=True,
                       default_value=0., **kwargs):
    """Calculate the quantity from that mode.

    Parameters
    ----------
        mode : ``pMuTT.statmech`` object
            Trans, Vib, Rot, Elec, or Nucl StatMech model
        method_name : str
            Name of method to use to calculate quantity.
        raise_error : bool, optional
            If True, raises an error if any of the modes do not have the 
            quantity of interest. Default is True
        raise_warning : bool, optional
            Only relevant if raise_error is False. Raises a warning if any
            of the modes do not have the quantity of interest. Default is
            True
        default_value : float, optional
            Default value if the object does not contain the method
        verbose : bool, optional
            If False, returns the total Gibbs energy. If True, returns
            contribution of each mode.
        kwargs : key-word arguments
            Parameters passed to each mode
    Returns
    -------
        quantity : float
            Quantity of the mode.
    Raises
    ------
        AttributeError
            If raise_error is True and the mode does not have the method_name
    """
    try:
        method = getattr(mode, method_name)
    except AttributeError as e:
        if raise_error:
            raise e
        elif raise_warning:
            warnings.warn(e, RuntimeWarning)
        quantity = default_value
    else:
        quantity = _pass_expected_arguments(method, **kwargs)
    return quantity
        

