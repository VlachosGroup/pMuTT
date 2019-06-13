# -*- coding: utf-8 -*-
import json


class pmuttEncoder(json.JSONEncoder):
    """Encodes pmutt objects to JSON format. The object (and complex subobjects
    ) must have the method: to_dict().
    """
    def default(self, o):
        try:
            o_dict = o.to_dict()
        except AttributeError:
            super().default(o)
        else:
            return o_dict


def json_to_pmutt(json_obj):
    """Object hook to convert json to pmutt objects. Any complex object should
    be in the :data:`~pmutt.io.json.type_to_class.type_to_class_dict`
    dictionary.

    Parameters
    ----------
        json_obj : Type supported by JSON
            JSON object to be converted to pmutt object. If this is a complex
            pmutt object, must have the 'class' entry in the dictionary.
    Returns
    -------
        obj : pmutt object
            Parsed pmutt object
    """
    try:
        class_name = type_to_class(json_obj['class'])
    except (KeyError, TypeError):
        return json_obj
    else:
        obj = class_name.from_dict(json_obj)
        return obj


def type_to_class(class_str):
    """Converts between type of object and pmutt classes.

    Parameters
    ----------
        class_str : str
            Output of str(pmutt_obj.__class__)
    Returns
    -------
        class : class
            Class corresponding to class_str
    """

    # Although it is usually inadvisible to import functions within a function,
    # this was done purposefully. Importing outside the function caused
    # circular import errors. This way the imports are limited to the function.
    from pmutt.eos import IdealGasEOS, vanDerWaalsEOS
    from pmutt.reaction import Reaction, Reactions
    from pmutt.reaction.bep import BEP
    from pmutt.empirical import EmpiricalBase
    from pmutt.empirical.nasa import Nasa
    from pmutt.empirical.shomate import Shomate
    from pmutt.empirical.references import Reference, References
    from pmutt.empirical.zacros import Zacros
    from pmutt.statmech import StatMech, EmptyMode
    from pmutt.statmech.trans import FreeTrans
    from pmutt.statmech.vib import HarmonicVib, QRRHOVib, EinsteinVib, DebyeVib
    from pmutt.statmech.rot import RigidRotor
    from pmutt.statmech.elec import GroundStateElec
    from pmutt.statmech.nucl import EmptyNucl
    from pmutt.mixture.cov import PiecewiseCovEffect
    from pmutt.chemkin import CatSite

    type_to_class_dict = {
        "<class 'pmutt.eos.IdealGasEOS'>": IdealGasEOS,
        "<class 'pmutt.eos.vanDerWaalsEOS'>": vanDerWaalsEOS,
        "<class 'pmutt.reaction.Reaction'>": Reaction,
        "<class 'pmutt.reaction.Reactions'>": Reactions,
        "<class 'pmutt.reaction.bep.BEP'>": BEP,
        "<class 'pmutt.empirical.EmpiricalBase'>": EmpiricalBase,
        "<class 'pmutt.empirical.nasa.Nasa'>": Nasa,
        "<class 'pmutt.empirical.shomate.Shomate'>": Shomate,
        "<class 'pmutt.empirical.references.Reference'>": Reference,
        "<class 'pmutt.empirical.references.References'>": References,
        "<class 'pmutt.empirical.zacros.Zacros'>": Zacros,
        "<class 'pmutt.statmech.StatMech'>": StatMech,
        "<class 'pmutt.statmech.EmptyMode'>": EmptyMode,
        "<class 'pmutt.statmech.trans.FreeTrans'>": FreeTrans,
        "<class 'pmutt.statmech.vib.HarmonicVib'>": HarmonicVib,
        "<class 'pmutt.statmech.vib.QRRHOVib'>": QRRHOVib,
        "<class 'pmutt.statmech.vib.EinsteinVib'>": EinsteinVib,
        "<class 'pmutt.statmech.vib.DebyeVib'>": DebyeVib,
        "<class 'pmutt.statmech.rot.RigidRotor'>": RigidRotor,
        "<class 'pmutt.statmech.elec.GroundStateElec'>": GroundStateElec,
        "<class 'pmutt.statmech.nucl.EmptyNucl'>": EmptyNucl,
        "<class 'pmutt.mixture.cov.PiecewiseCovEffect'>": PiecewiseCovEffect,
        "<class 'pmutt.chemkin.CatSite'>": CatSite,
    }
    return type_to_class_dict[class_str]


def remove_class(json_obj):
    """Removes unnecessary entries from the JSON object when reinitializing the
    pmutt object

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
