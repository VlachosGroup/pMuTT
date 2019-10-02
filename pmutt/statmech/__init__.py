# -*- coding: utf-8 -*-
"""
pmutt.statmech
Vlachos group code for thermodynamic models.
Created on Fri Jul 7 12:40:00 2018
"""
import inspect
from copy import copy
import numpy as np
from pmutt import (_pass_expected_arguments, _is_iterable, _get_mode_quantity,
                   _get_specie_kwargs, _apply_numpy_operation, _ModelBase,
                   _get_R_adj, parse_formula)
from pmutt import constants as c
from pmutt.statmech import trans, vib, elec, rot
from pmutt.mixture import _get_mix_quantity
from pmutt.io import json as json_pmutt


class EmptyMode(_ModelBase):
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

class ConstantMode(_ModelBase):
    """Mode where thermodynamic properties can be arbitrarily set. Note that
    thermodynamic properties must be explicitly set and are not calculated from
    other properties.
 
    Attributes
    ----------
        q : float, optional
            Partition function. Default is 1
        Cv : float, optional
            Heat capacity at constant volume in eV/K. Default is 0
        Cp : float, optional
            Heat capacity at constant pressure in eV/K. Default is 0
        U : float, optional
            Internal energy in eV. Default is 0
        H : float, optional
            Enthalpy in eV. Default is 0
        S : float, optional
            Entropy in eV/K. Default is 0
        F : float, optional
            Helmholtz energy in eV. Default is 0
        G : float, optional
            Gibbs energy in eV. Default is 0
        notes : str or dict, optional
            Any additional details you would like to include such as
            source of data. Default is None
            
    """
    def __init__(self, q=1., Cv=0., Cp=0., U=0., H=0., S=0., F=0., G=0.,
                 notes=None):
        self.q = q
        self.Cv = Cv
        self.Cp = Cp
        self.U = U
        self.H = H
        self.S = S
        self.F = F
        self.G = G
        self.notes = notes

    def get_q(self):
        """Calculate the partition function
        
        Returns
        -------
            q : float
                Partition function
        """
        return self.q

    def get_CvoR(self):
        """Calculate the dimensionless heat capacity (constant volume)
        
        Returns
        -------
            CvoR : float
                Dimensionless heat capacity (constant volume)
        """
        return self.Cv/c.R('eV/K')

    def get_CpoR(self):
        """Calculate the dimensionless heat capacity (constant pressure)
        
        Returns
        -------
            CpoR : float
                Dimensionless heat capacity (constant pressure)
        """
        return self.Cp/c.R('eV/K')

    def get_UoRT(self, T=c.T0('K')):
        """Calculate the dimensionless internal energy

        Parameters
        ----------
            T : float, optional
                Temperature in K. Default is 298.15 K
        Returns
        -------
            UoRT : float
                Dimensionless internal energy
        """
        return self.U/c.R('eV/K')/T

    def get_HoRT(self, T=c.T0('K')):
        """Calculate the dimensionless enthalpy

        Parameters
        ----------
            T : float, optional
                Temperature in K. Default is 298.15 K
        Returns
        -------
            HoRT : float
                Dimensionless enthalpy
        """
        return self.H/c.R('eV/K')/T

    def get_SoR(self):
        """Calculate the dimensionless entropy
        
        Returns
        -------
            SoR : float
                Dimensionless entropy
        """
        return self.S/c.R('eV/K')

    def get_FoRT(self, T=c.T0('K')):
        """Calculate the dimensionless Helmholtz energy
        
        Parameters
        ----------
            T : float, optional
                Temperature in K. Default is 298.15 K
        Returns
        -------
            FoRT : float
                Dimensionless Helmholtz energy
        """
        return self.F/c.R('eV/K')/T

    def get_GoRT(self, T=c.T0('K')):
        """Calculate the dimensionless Gibbs energy
        
        Parameters
        ----------
            T : float, optional
                Temperature in K. Default is 298.15 K
        Returns
        -------
            GoRT : float
                Dimensionless Gibbs energy
        """
        return self.G/c.R('eV/K')/T

class StatMech(_ModelBase):
    """Base class for statistical mechanic models.

    Attributes
    ----------
        name : str, optional
            Name of the specie. Default is None
        trans_model : :ref:`pmutt.statmech.trans <trans>` object, optional
            Deals with translational modes. Default is
            :class:`~pmutt.statmech.EmptyMode`
        vib_model : :ref:`pmutt.statmech.vib <vib>` object, optional
            Deals with vibrational modes. Default is
            :class:`~pmutt.statmech.EmptyMode`
        rot_model : :ref:`pmutt.statmech.rot <rot>` object, optional
            Deals with rotational modes. Default is
            :class:`~pmutt.statmech.EmptyMode`
        elec_model : :ref:`pmutt.statmech.elec <elec>` object, optional
            Deals with electronic modes. Default is
            :class:`~pmutt.statmech.EmptyMode`
        nucl_model : :ref:`pmutt.statmech.nucl <nucl>` object
            Deals with nuclear modes. Default is
            :class:`~pmutt.statmech.EmptyMode`
        elements : dict, optional
            Composition of the species. Default is None.
            Keys of dictionary are elements, values are stoichiometric values
            in a formula unit.
            e.g. CH3OH can be represented as:
            {'C': 1, 'H': 4, 'O': 1,}.
        references : :class:`~pmutt.empirical.references.References`, optional
            Contains references to adjust ``HoRT``. If not specified
            then HoRT_dft will be used without adjustment. Default is None
        misc_models : list of pmutt model objects, optional
            Extra models that do not fit in the above attributes. Commonly used
            models would be ``pmutt.statmech`` and ``pmutt.mixture`` models
        smiles : str, optional
            Smiles representation of species
        notes : str, optional
            Any additional details you would like to include such as
            computational set up. Default is None
    """

    def __init__(self, name=None, trans_model=EmptyMode(),
                 vib_model=EmptyMode(), rot_model=EmptyMode(),
                 elec_model=EmptyMode(), nucl_model=EmptyMode(),
                 misc_models=None, elements=None, references=None,
                 smiles=None, notes=None, **kwargs):
        self.name = name
        self.smiles = smiles
        self.notes = notes
        self.references = references

        # Assign elements from Atoms
        if elements is None:
            try:
                elements = parse_formula( \
                        kwargs['atoms'].get_chemical_formula('hill'))
            except KeyError:
                pass
        self.elements=elements


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

        # Assign misc models
        # TODO misc models can not be initialized by passing the class
        # because all the models will have the same attributes. Figure out
        # a way to pass them. Perhaps have a dictionary that contains the
        # attributes separated by species
        if not _is_iterable(misc_models) and misc_models is not None:
            misc_models = [misc_models]
        self.misc_models = misc_models

    def get_quantity(self, method_name, raise_error=True, raise_warning=True,
                     operation='sum', verbose=False, use_references=True,
                     **kwargs):
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            verbose : bool, optional
                If False, returns the total Gibbs energy. If True, returns
                contribution of each mode.
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            quantity : float or (N+6,) `numpy.ndarray`_
                Desired quantity. N represents the number of misc models. If
                verbose is True, contribution to each mode are as follows:
                [trans, vib, rot, elec, nucl, references, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        # Get the default value
        operation = operation.lower()
        if operation == 'sum':
            default_value = 0.
        elif operation == 'prod':
            default_value = 1.
        else:
            err_msg = 'Operation: {} not supported'.format(operation)
            raise ValueError(err_msg)

        # Calculate the quantity for each mode
        specie_kwargs = _get_specie_kwargs(specie_name=self.name, **kwargs)
        quantity = np.array([
            _get_mode_quantity(mode=self.trans_model, method_name=method_name,
                               raise_error=raise_error,
                               raise_warning=raise_warning,
                               default_value=default_value,
                               **specie_kwargs),
            _get_mode_quantity(mode=self.vib_model, method_name=method_name,
                               raise_error=raise_error,
                               raise_warning=raise_warning,
                               default_value=default_value,
                               **specie_kwargs),
            _get_mode_quantity(mode=self.rot_model, method_name=method_name,
                               raise_error=raise_error,
                               raise_warning=raise_warning,
                               default_value=default_value,
                               **specie_kwargs),
            _get_mode_quantity(mode=self.elec_model, method_name=method_name,
                               raise_error=raise_error,
                               raise_warning=raise_warning,
                               default_value=default_value,
                               **specie_kwargs),
            _get_mode_quantity(mode=self.nucl_model, method_name=method_name,
                               raise_error=raise_error,
                               raise_warning=raise_warning,
                               default_value=default_value,
                               **specie_kwargs)])
        if use_references and self.references is not None:
            ref_kwargs = copy(specie_kwargs)
            ref_kwargs['descriptors'] = getattr(self,
                                                self.references.descriptor)
            refs_quantity = np.atleast_1d(
                    _get_mode_quantity(mode=self.references,
                                       method_name=method_name,
                                       raise_error=raise_error,
                                       raise_warning=raise_warning,
                                       default_value=default_value,
                                       **ref_kwargs))
        else:
            refs_quantity = np.atleast_1d(default_value)
        # Calculate contribution from misc models if any
        misc_quantity = _get_mix_quantity(misc_models=self.misc_models,
                                         method_name=method_name,
                                         raise_error=raise_error,
                                         raise_warning=raise_warning,
                                         default_value=default_value,
                                         verbose=verbose,
                                         **kwargs)
        # Add misc quantities onto quantity
        quantity = np.concatenate([quantity, refs_quantity, misc_quantity])
        quantity = _apply_numpy_operation(quantity, verbose=verbose,
                                          operation=operation)
        return quantity

    def get_q(self, verbose=False, raise_error=True, raise_warning=True,
              use_references=True, **kwargs):
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            q : float or (N+5,) `numpy.ndarray`_
                Partition function. N represents the number of misc models.
                If verbose is True, contribution to each mode are as follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_q', raise_error=raise_error,
                                 raise_warning=raise_warning, operation='prod',
                                 verbose=verbose, use_references=use_references,
                                 **kwargs)

    def get_CvoR(self, verbose=False, raise_error=True, raise_warning=True,
                 use_references=True, **kwargs):
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            CvoR : float or (N+5,) `numpy.ndarray`_
                Dimensionless heat capacity. N represents the number of misc
                models. If verbose is True, contribution to each mode are as
                follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_CvoR', operation='sum',
                                 raise_error=raise_error,
                                 raise_warning=raise_warning,
                                 verbose=verbose, use_references=use_references,
                                 **kwargs)

    def get_Cv(self, units, verbose=False, raise_error=True,
               raise_warning=True, use_references=True, **kwargs):
        """Calculate the heat capacity (constant V)

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            Cv : float or (N+5,) `numpy.ndarray`_
                Heat capacity. N represents the number of misc models. If
                verbose is True, contribution to each mode are as follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        R_adj = _get_R_adj(units=units, elements=self.elements)
        return self.get_CvoR(verbose=verbose, raise_error=raise_error,
                             raise_warning=raise_warning,
                             use_references=use_references, **kwargs)*R_adj

    def get_CpoR(self, verbose=False, raise_error=True, raise_warning=True,
                 use_references=True, **kwargs):
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            CpoR : float or (N+5,) `numpy.ndarray`_
                Dimensionless heat capacity. N represents the number of misc
                models. If verbose is True, contribution to each mode are as
                follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_CpoR', operation='sum',
                                 raise_error=raise_error,
                                 raise_warning=raise_warning,
                                 verbose=verbose, use_references=use_references,
                                 **kwargs)

    def get_Cp(self, units, verbose=False, raise_error=True,
               raise_warning=True, use_references=True, **kwargs):
        """Calculate the heat capacity (constant P)

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            Cp : float or (N+5,) `numpy.ndarray`_
                Heat capacity. N represents the number of misc models. If
                verbose is True, contribution to each mode are as follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        R_adj = _get_R_adj(units=units, elements=self.elements)
        return self.get_CpoR(verbose=verbose, raise_error=raise_error,
                             raise_warning=raise_warning,
                             use_references=use_references, **kwargs)*R_adj

    def get_EoRT(self, T=c.T0('K'), include_ZPE=False,
                 raise_error=True, raise_warning=True, **kwargs):
        """Dimensionless electronic energy

        Parameters
        ----------
            T : float, optional
                Temperature in K. If the electronic mode is
                :class:`~pmutt.statmech.elec.GroundStateElec`, then the output is
                insensitive to this input. Default is 298.15 K
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
                Dimensionless electronic energy
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
                                       default_value=0.,
                                       **kwargs)/c.R('eV/K')/T
        return EoRT

    def get_E(self, units, T=c.T0('K'), raise_error=True, raise_warning=True,
              include_ZPE=False, **kwargs):
        """Calculate the electronic energy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            T : float, optional
                Temperature in K. If the electronic mode is
                :class:`~pmutt.statmech.elec.GroundStateElec`, then the output is
                insensitive to this input. Default is 298.15 K
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
        units = '{}/K'.format(units)
        R_adj = _get_R_adj(units=units, elements=self.elements)
        return self.get_EoRT(T=T, raise_error=raise_error,
                             raise_warning=raise_warning,
                             include_ZPE=include_ZPE, **kwargs)*T*R_adj
                             
    def get_UoRT(self, verbose=False, raise_error=True, raise_warning=True,
                 use_references=True, **kwargs):
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            UoRT : float or (N+5,) `numpy.ndarray`_
                Dimensionless internal energy. N represents the number of
                misc models. If verbose is True, contribution to each mode
                are as follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_UoRT', operation='sum',
                                 raise_error=raise_error,
                                 raise_warning=raise_warning,
                                 verbose=verbose,
                                 use_references=use_references, **kwargs)

    def get_U(self, units, T=c.T0('K'), verbose=False, raise_error=True,
              raise_warning=True, use_references=True, **kwargs):
        """Calculate the internal energy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            U : float or (N+5,) `numpy.ndarray`_
                Internal energy. N represents the number of misc models. If
                verbose is True, contribution to each mode are as follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        units = '{}/K'.format(units)
        R_adj = _get_R_adj(units=units, elements=self.elements)
        return self.get_UoRT(verbose=verbose, T=T, raise_error=raise_error,
                             raise_warning=raise_warning,
                             use_references=use_references, **kwargs)*T*R_adj

    def get_HoRT(self, verbose=False, raise_error=True, raise_warning=True,
                 use_references=True, **kwargs):
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            HoRT : float or (N+5,) `numpy.ndarray`_
                Dimensionless enthalpy. N represents the number of misc
                models. If verbose is True, contribution to each mode are as
                follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_HoRT', operation='sum',
                                 raise_error=raise_error,
                                 raise_warning=raise_warning,
                                 verbose=verbose, use_references=use_references,
                                 **kwargs)

    def get_H(self, units, T=c.T0('K'), raise_error=True, raise_warning=True,
              verbose=False, use_references=True, **kwargs):
        """Calculate the enthalpy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            H : float or (N+5,) `numpy.ndarray`_
                Enthalpy. N represents the number of misc models. If
                verbose is True, contribution to each mode are as follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        units = '{}/K'.format(units)
        R_adj = _get_R_adj(units=units, elements=self.elements)
        return self.get_HoRT(verbose=verbose, raise_error=raise_error,
                             raise_warning=raise_warning, T=T,
                             use_references=use_references, **kwargs)*T*R_adj

    def get_SoR(self, verbose=False, raise_error=True, raise_warning=True,
                use_references=True, **kwargs):
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
        Returns
        -------
            SoR : float or (N+5,) `numpy.ndarray`_
                Dimensionless entropy. N represents the number of misc
                models. If verbose is True, contribution to each mode are as
                follows: [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_SoR', operation='sum',
                                 raise_error=raise_error,
                                 raise_warning=raise_warning,
                                 verbose=verbose, use_references=use_references,
                                 **kwargs)

    def get_S(self, units, verbose=False, raise_error=True,
              raise_warning=True, use_references=True, **kwargs):
        """Calculate the entropy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            S : float or (N+5,) `numpy.ndarray`_
                Entropy. N represents the number of misc models. If
                verbose is True, contribution to each mode are as follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        R_adj = _get_R_adj(units=units, elements=self.elements)
        return self.get_SoR(verbose=verbose, raise_error=raise_error,
                            raise_warning=raise_warning, 
                            use_references=use_references, **kwargs)*R_adj

    def get_FoRT(self, verbose=False, raise_error=True, raise_warning=True,
                 use_references=True, **kwargs):
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            FoRT : float or (N+5,) `numpy.ndarray`_
                Dimensionless Helmoltz energy. N represents the number of
                misc models. If verbose is True, contribution to each mode
                are as follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_FoRT', operation='sum',
                                 raise_error=raise_error,
                                 raise_warning=raise_warning,
                                 verbose=verbose, use_references=use_references,
                                 **kwargs)

    def get_F(self, units, T=c.T0('K'), verbose=False, raise_error=True,
              raise_warning=True, use_references=True, **kwargs):
        """Calculate the Helmholtz energy

        Parameters
        ----------
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            F : float or (N+5,) `numpy.ndarray`_
                Helmholtz energy. N represents the number of misc models. If
                verbose is True, contribution to each mode are as follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        units = '{}/K'.format(units)
        R_adj = _get_R_adj(units=units, elements=self.elements)
        return self.get_FoRT(verbose=verbose, T=T, raise_error=raise_error,
                             raise_warning=raise_warning,
                             use_references=use_references, **kwargs)*T*R_adj

    def get_GoRT(self, verbose=False, raise_error=True, raise_warning=True,
                 use_references=True, **kwargs):
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
            use_references : bool, optional
                If True, adds contribution from references. Default is True
            kwargs : key-word arguments
                Parameters passed to each mode
        Returns
        -------
            GoRT : float or (N+5,) `numpy.ndarray`_
                Dimensionless Gibbs energy. N represents the number of misc
                models. If verbose is True, contribution to each mode are as
                follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        return self.get_quantity(method_name='get_GoRT', operation='sum',
                                 raise_error=raise_error,
                                 raise_warning=raise_warning,
                                 verbose=verbose, use_references=use_references,
                                 **kwargs)

    def get_G(self, units, T=c.T0('K'), verbose=False, raise_error=True,
              raise_warning=True, use_references=True, **kwargs):
        """Calculate the Gibbs energy

        Parameters
        ----------
            verbose : bool, optional
                If False, returns the Gibbs energy. If True, returns
                contribution of each mode.
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
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
            G : float or (N+5,) `numpy.ndarray`_
                Gibbs energy. N represents the number of misc models. If
                verbose is True, contribution to each mode are as follows:
                [trans, vib, rot, elec, nucl, misc_models (if any)]

        .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
        """
        units = '{}/K'.format(units)
        R_adj = _get_R_adj(units=units, elements=self.elements)
        return self.get_GoRT(verbose=verbose, raise_error=raise_error,
                             raise_warning=raise_warning, T=T,
                             use_references=use_references, **kwargs)*T*R_adj

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = {
            'class': str(self.__class__),
            'name': self.name,
            'trans_model': self.trans_model.to_dict(),
            'vib_model': self.vib_model.to_dict(),
            'rot_model': self.rot_model.to_dict(),
            'elec_model': self.elec_model.to_dict(),
            'nucl_model': self.nucl_model.to_dict(),
            'smiles': self.smiles,
            'notes': self.notes}

        if _is_iterable(self.misc_models):
            obj_dict['misc_models'] = \
                    [mix_model.to_dict() for mix_model in self.misc_models]
        else:
            obj_dict['misc_models'] = self.misc_models

        try:
            obj_dict['references'] = self.references.to_dict()
        except AttributeError:
            obj_dict['references'] = self.references

        return obj_dict

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
        json_obj = json_pmutt.remove_class(json_obj)
        name = json_obj['name']
        trans_model = json_pmutt.json_to_pmutt(json_obj['trans_model'])
        vib_model = json_pmutt.json_to_pmutt(json_obj['vib_model'])
        rot_model = json_pmutt.json_to_pmutt(json_obj['rot_model'])
        elec_model = json_pmutt.json_to_pmutt(json_obj['elec_model'])
        nucl_model = json_pmutt.json_to_pmutt(json_obj['nucl_model'])
        references = json_pmutt.json_to_pmutt(json_obj['references'])
        notes = json_obj['notes']

        return cls(name=name,
                   trans_model=trans_model,
                   vib_model=vib_model,
                   rot_model=rot_model,
                   elec_model=elec_model,
                   nucl_model=nucl_model,
                   references=references,
                   notes=notes)


presets = {
    'idealgas': {
        'model': StatMech,
        'trans_model': trans.FreeTrans,
        'n_degrees': 3,
        'vib_model': vib.HarmonicVib,
        'elec_model': elec.GroundStateElec,
        'rot_model': rot.RigidRotor,
        'required': ('molecular_weight', 'vib_wavenumbers', 'potentialenergy',
                     'spin', 'geometry', 'rot_temperatures', 'symmetrynumber'),
        'optional': ('atoms')
        },
    'harmonic': {
        'model': StatMech,
        'vib_model': vib.HarmonicVib,
        'elec_model': elec.GroundStateElec,
        'required': ('vib_wavenumbers', 'potentialenergy', 'spin'),
    },
    'electronic': {
        'model': StatMech,
        'elec_model': elec.GroundStateElec,
        'required': ('potentialenergy', 'spin'),
    },
    'placeholder': {
        'model': StatMech,
        'trans_model': EmptyMode,
        'vib_model': EmptyMode,
        'elec_model': EmptyMode,
        'rot_model': EmptyMode,
        'nucl_model': EmptyMode,
        'required': (),
    },
    'constant': {
        'model': StatMech,
        'elec_model': ConstantMode,
        'optional': ('q', 'Cv', 'Cp', 'U', 'H', 'S', 'F', 'G')
    }
}
"""Commonly used models. The 'required' and 'optional' fields indicate
parameters that still need to be passed."""
