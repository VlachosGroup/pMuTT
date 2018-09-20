# -*- coding: utf-8 -*-
import json

class PyMuTTEncoder(json.JSONEncoder):
    """Encodes PyMuTT objects to JSON format. The object (and complex subobjects
    ) must have the method: ``to_dict()``.
    """
    def default(self, obj):
        try:
            obj_dict = obj.to_dict()
        except AttributeError:
            super().default(obj)
        else:
            return obj_dict

def json_to_PyMuTT(json_obj):
    """Object hook to convert json to PyMuTT objects. Any complex object should
    be in the ``type_to_class_dict`` in ``PyMuTT.io_.jsonio.type_to_class``
    function.

    Parameters
    ----------
        json_obj : Type supported by JSON
            JSON object to be converted to PyMuTT object. If this is a complex
            PyMuTT object, must have the 'class' entry in the dictionary.
    Returns
    -------
        obj : PyMuTT object
            Parsed PyMuTT object
    """
    try:
        class_name = type_to_class(json_obj['class'])
    except (KeyError, TypeError):
        return json_obj
    else:
        obj = class_name.from_dict(json_obj)
        return obj

def type_to_class(class_str):
    """Converts between type of object and PyMuTT classes.

    Parameters
    ----------
        class_str : str
            Output of str(PyMuTT_obj.__class__)
    Returns
    -------
        class : class
            Class corresponding to class_str
    """

    # Although it is usually inadvisible to import functions within a function,
    # this was done purposefully. Using only a dictionary resulted in circular
    # import errors. This way the imports are limited to the function.
    from PyMuTT.models.empirical import BaseThermo
    from PyMuTT.models.empirical.nasa import Nasa
    from PyMuTT.models.empirical.references import References
    from PyMuTT.models.empirical.zacros import Zacros
    from PyMuTT.models.statmech import StatMech, EmptyMode
    from PyMuTT.models.statmech.trans import IdealTrans
    from PyMuTT.models.statmech.vib import HarmonicVib, QRRHOVib
    from PyMuTT.models.statmech.rot import RigidRotor
    from PyMuTT.models.statmech.elec import IdealElec
    from PyMuTT.models.statmech.nucl import IdealNucl

    type_to_class_dict = {
        "<class 'PyMuTT.models.empirical.BaseThermo'>": BaseThermo,
        "<class 'PyMuTT.models.empirical.nasa.Nasa'>": Nasa,
        "<class 'PyMuTT.models.empirical.references.References'>": References,
        "<class 'PyMuTT.models.empirical.zacros.Zacros'>": Zacros,
        "<class 'PyMuTT.models.statmech.StatMech'>": StatMech,
        "<class 'PyMuTT.models.statmech.EmptyMode'>": EmptyMode,
        "<class 'PyMuTT.models.statmech.trans.IdealTrans'>": IdealTrans,
        "<class 'PyMuTT.models.statmech.vib.HarmonicVib'>": HarmonicVib,
        "<class 'PyMuTT.models.statmech.vib.QRRHOVib'>": QRRHOVib,
        "<class 'PyMuTT.models.statmech.rot.RigidRotor'>": RigidRotor,
        "<class 'PyMuTT.models.statmech.elec.IdealElec'>": IdealElec,
        "<class 'PyMuTT.models.statmech.nucl.IdealNucl'>": IdealNucl,
    }
    return type_to_class_dict[class_str]