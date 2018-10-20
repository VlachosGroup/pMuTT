# -*- coding: utf-8 -*-
import json

class pMuTTEncoder(json.JSONEncoder):
    """Encodes pMuTT objects to JSON format. The object (and complex subobjects
    ) must have the method: ``to_dict()``.
    """
    def default(self, o):
        try:
            o_dict = o.to_dict()
        except AttributeError:
            super().default(o)
        else:
            return o_dict

def json_to_pMuTT(json_obj):
    """Object hook to convert json to pMuTT objects. Any complex object should
    be in the ``type_to_class_dict`` in ``pMuTT.io_.jsonio.type_to_class``
    function.

    Parameters
    ----------
        json_obj : Type supported by JSON
            JSON object to be converted to pMuTT object. If this is a complex
            pMuTT object, must have the 'class' entry in the dictionary.
    Returns
    -------
        obj : pMuTT object
            Parsed pMuTT object
    """
    try:
        class_name = type_to_class(json_obj['class'])
    except (KeyError, TypeError):
        return json_obj
    else:
        obj = class_name.from_dict(json_obj)
        return obj

def type_to_class(class_str):
    """Converts between type of object and pMuTT classes.

    Parameters
    ----------
        class_str : str
            Output of str(pMuTT_obj.__class__)
    Returns
    -------
        class : class
            Class corresponding to class_str
    """

    # Although it is usually inadvisible to import functions within a function,
    # this was done purposefully. Importing outside the function caused circular
    # import errors. This way the imports are limited to the function.
    from pMuTT.models.eos import IdealGasEOS, vanDerWaalsEOS
    from pMuTT.models.reaction import Reaction
    from pMuTT.models.empirical import BaseThermo
    from pMuTT.models.empirical.nasa import Nasa
    from pMuTT.models.empirical.references import Reference, References
    from pMuTT.models.empirical.zacros import Zacros
    from pMuTT.models.statmech import StatMech, EmptyMode
    from pMuTT.models.statmech.trans import IdealTrans
    from pMuTT.models.statmech.vib import HarmonicVib, QRRHOVib, EinsteinVib
    from pMuTT.models.statmech.rot import RigidRotor
    from pMuTT.models.statmech.elec import IdealElec
    from pMuTT.models.statmech.nucl import IdealNucl

    type_to_class_dict = {
        "<class 'pMuTT.models.eos.IdealGasEOS'>": IdealGasEOS,
        "<class 'pMuTT.models.eos.vanDerWaalsEOS'>": vanDerWaalsEOS,
        "<class 'pMuTT.models.reaction.Reaction'>": Reaction,
        "<class 'pMuTT.models.empirical.BaseThermo'>": BaseThermo,
        "<class 'pMuTT.models.empirical.nasa.Nasa'>": Nasa,
        "<class 'pMuTT.models.empirical.references.Reference'>": Reference,
        "<class 'pMuTT.models.empirical.references.References'>": References,
        "<class 'pMuTT.models.empirical.zacros.Zacros'>": Zacros,
        "<class 'pMuTT.models.statmech.StatMech'>": StatMech,
        "<class 'pMuTT.models.statmech.EmptyMode'>": EmptyMode,
        "<class 'pMuTT.models.statmech.trans.IdealTrans'>": IdealTrans,
        "<class 'pMuTT.models.statmech.vib.HarmonicVib'>": HarmonicVib,
        "<class 'pMuTT.models.statmech.vib.QRRHOVib'>": QRRHOVib,
        "<class 'pMuTT.models.statmech.vib.EinsteinVib'>": EinsteinVib,
        "<class 'pMuTT.models.statmech.rot.RigidRotor'>": RigidRotor,
        "<class 'pMuTT.models.statmech.elec.IdealElec'>": IdealElec,
        "<class 'pMuTT.models.statmech.nucl.IdealNucl'>": IdealNucl,
    }
    return type_to_class_dict[class_str]

def remove_class(json_obj):
    """Removes the 'class' entry from the JSON object.

    Parameters
    ----------
        json_obj : dict
            JSON object with 'class' entry
    Returns
    -------
        json_obj : dict
            JSON object without 'class' entry
    """
    try:
        del json_obj['class']
    except KeyError:
        pass
    return json_obj