"""
pMuTT.mixture
Vlachos group code for handling the effect of mixing on thermodynamic
properties
Created on Fri Feb 8 14:30:00 2018
"""

import numpy as np
from pMuTT import _get_specie_kwargs, _apply_numpy_operation, _get_mode_quantity


def _get_mix_quantity(misc_models, method_name, raise_error=True,
                      raise_warning=True, default_value=0., **kwargs):
    """Calculate contribution from mixing models to desired quantity

    Parameters
    ----------
        misc_models : list (length N) of ``pMuTT.mixture`` objects
            Mix models to calculate the property
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
        default_value : float, optional
            Default value if the object does not contain the method. Default is
            0
        kwargs : key-word arguments
            Arguments to calculate mixture model properties, if any
    Returns
    -------
        mix_quantity : (N,) `numpy.ndarray`_
            Mixing quantity of interest. If verbose is True, each element
            corresponds to the contribution of each mix_model

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ndarray.html
    """
    # Return default value if no mixture models exist
    if misc_models is None:
        return np.array([default_value])

    # Calculate contribution from mixing models if any
    mix_quantity = np.full_like(a=misc_models, fill_value=default_value)
    for i, mix_model in enumerate(misc_models):
        if mix_model is None:
            continue

        specie_kwargs = _get_specie_kwargs(mix_model.name_j, **kwargs)
        mix_quantity[i] = _get_mode_quantity(mode=mix_model,
                                             method_name=method_name,
                                             raise_error=raise_error,
                                             raise_warning=raise_warning,
                                             default_value=default_value,
                                             **specie_kwargs)
    return mix_quantity
