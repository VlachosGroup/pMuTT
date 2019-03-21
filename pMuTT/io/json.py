# -*- coding: utf-8 -*-
import json


class pMuTTEncoder(json.JSONEncoder):
    """Encodes pMuTT objects to JSON format. The object (and complex subobjects
    ) must have the method: to_dict().
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
    be in the :data:`~pMuTT.io.json.type_to_class.type_to_class_dict`
    dictionary.

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
    # this was done purposefully. Importing outside the function caused
    # circular import errors. This way the imports are limited to the function.
    from pMuTT.eos import IdealGasEOS, vanDerWaalsEOS
    from pMuTT.reaction import Reaction, Reactions
    from pMuTT.reaction.bep import BEP
    from pMuTT.empirical import EmpiricalBase
    from pMuTT.empirical.nasa import Nasa
    from pMuTT.empirical.shomate import Shomate
    from pMuTT.empirical.references import Reference, References
    from pMuTT.empirical.zacros import Zacros
    from pMuTT.statmech import StatMech, EmptyMode
    from pMuTT.statmech.trans import IdealTrans
    from pMuTT.statmech.vib import HarmonicVib, QRRHOVib, EinsteinVib
    from pMuTT.statmech.rot import RigidRotor
    from pMuTT.statmech.elec import IdealElec
    from pMuTT.statmech.nucl import IdealNucl
    from pMuTT.mixture.cov import CovEffect
    from pMuTT.chemkin import CatSite

    type_to_class_dict = {
        "<class 'pMuTT.eos.IdealGasEOS'>": IdealGasEOS,
        "<class 'pMuTT.eos.vanDerWaalsEOS'>": vanDerWaalsEOS,
        "<class 'pMuTT.reaction.Reaction'>": Reaction,
        "<class 'pMuTT.reaction.Reactions'>": Reactions,
        "<class 'pMuTT.reaction.bep.BEP'>": BEP,
        "<class 'pMuTT.empirical.EmpiricalBase'>": EmpiricalBase,
        "<class 'pMuTT.empirical.nasa.Nasa'>": Nasa,
        "<class 'pMuTT.empirical.shomate.Shomate'>": Shomate,
        "<class 'pMuTT.empirical.references.Reference'>": Reference,
        "<class 'pMuTT.empirical.references.References'>": References,
        "<class 'pMuTT.empirical.zacros.Zacros'>": Zacros,
        "<class 'pMuTT.statmech.StatMech'>": StatMech,
        "<class 'pMuTT.statmech.EmptyMode'>": EmptyMode,
        "<class 'pMuTT.statmech.trans.IdealTrans'>": IdealTrans,
        "<class 'pMuTT.statmech.vib.HarmonicVib'>": HarmonicVib,
        "<class 'pMuTT.statmech.vib.QRRHOVib'>": QRRHOVib,
        "<class 'pMuTT.statmech.vib.EinsteinVib'>": EinsteinVib,
        "<class 'pMuTT.statmech.rot.RigidRotor'>": RigidRotor,
        "<class 'pMuTT.statmech.elec.IdealElec'>": IdealElec,
        "<class 'pMuTT.statmech.nucl.IdealNucl'>": IdealNucl,
        "<class 'pMuTT.mixture.cov.CovEffect'>": CovEffect,
        "<class 'pMuTT.chemkin.CatSite'>": CatSite,
    }
    return type_to_class_dict[class_str]


def remove_class(json_obj):
    """Removes unnecessary entries from the JSON object when reinitializing the
    pMuTT object

    Parameters
    ----------
        json_obj : dict
            JSON object unnecessary entries
    Returns
    -------
        json_obj : dict
            JSON object without unnecessary entries
    """
    json_obj.pop('class', None)
    json_obj.pop('type', None)
    json_obj.pop('_id', None)
    return json_obj
