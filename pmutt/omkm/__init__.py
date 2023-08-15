from collections import namedtuple

_Param = namedtuple('_Param', 'label val units')
"""Parameters as a NamedTuple for easier unit passing for YAML file.

    Attributes
    ----------
        label : str
            Label name for attribute being assigned.
        val : obj
            Value that is assigned to header[label] if the key does not exist or
            val is None.
        val_units : str
            Units for ``val`` where quantities are proceeded by '_'.
            e.g. '_length3/_time' for the volumetric flow rate.
"""

def _assign_yaml_val(param, header, units=None):
    """Helper method to assign label to header
    
    Parameters
    ----------
        param : _Param namedtuple
            Parameter with three attributes: ``label``, ``val``, ``units``
        header : dict
            Upper level dictionary that ``label`` will be nested under.
        units : :class:`~pmutt.omkm.units.Unit` object, optional
            Units to write file.
    """
    # Do nothing if the label was previously assigned
    if param.label in header:
        return
    # Do nothing if value is not specified
    if param.val is None:
        return
    if isinstance(param.val, str):
        param = param._replace(val='\"{}\"'.format(param.val))
    # Assign the value as is if the unit type is None
    if param.units is None:
        header[param.label] = param.val
        return
    # Assume SI units if units is not specified
    if param.units is None:
        header[param.label] = param.val
        return
    # If the value is numerical and units were specified, add the units
    if isinstance(param.val, (int, float)):
        val_str = '\"{} {}\"'.format(param.val, param.units)
        # Replace activation energy first to avoid complications with energy
        val_str = val_str.replace('_act_energy', units.act_energy)
        # Add appropriate units
        for unit_type, unit in units.__dict__.items():
            val_str = val_str.replace('_{}'.format(unit_type), unit)
        header[param.label] = val_str
    # If the value is a list and units were specified, add units to each entry
    elif isinstance(param.val, list):
        vals_list = ['\"{} {}\"'.format(i, param.units) for i in param.val]
        for unit_type, unit in units.__dict__.items():
            old_str = '_{}'.format(unit_type)
            for i, val in enumerate(vals_list):
                vals_list[i] = val.replace(old_str, unit)
        header[param.label] = vals_list
