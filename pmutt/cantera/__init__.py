import more_itertools as mit


def _get_range_CTI(objs, parent_obj=None, delimiter='_'):
    """Returns the objs IDs
    
    Parameters
    ----------
        objs : list
            List of objects to be grouped as a range of CTI strings. Objects
            should be str, or have the ``name`` or ``id`` attribute.
        parent_obj : obj, optional
            Primarily used for more helpful error messages. Default is None.
        delimiter : str, optional
            Delimiter to separate header and footer of obj strings. Default is
            '_'
    Returns
    -------
        CTI_range_out : str
            Objects expressed in CTI range format.
    Raises
    ------
        TypeError:
            Raised when members of ``objs`` do not have an ``id`` or ``name``
            attribute and cannot be expressed as a string.
        ValueError:
            Raised when ``objs`` does not have an integer-compatible section
            that is required by the external library, more_itertools.
    """
    if isinstance(objs, str):
        CTI_out = objs
    elif len(objs) == 0:
        CTI_out = '[]'
    elif objs is None:
        CTI_out = '[]'
    else:
        CTI_out = '['

        # Get unique ids
        unique_headers = {}
        for obj in objs:
            # Try getting obj id from obj object
            try:
                obj_id = obj.id
            except AttributeError:
                try:
                    obj_id = obj.name
                except AttributeError:
                    obj_id = obj
            # Check that obj_id is a string type
            if not isinstance(obj_id, str):
                err_msg = (
                    'Expected obj ID for obj "{}" in {} to be str type. '
                    'Instead, received "{}" (type {}).'
                    ''.format(str(obj), parent_obj.__class__, obj_id,
                              type(obj_id)))
                raise TypeError(err_msg)

            # Separate header and footer
            i = obj_id.rfind(delimiter)
            # If delimiter does not exist, assign to generic list
            if i == -1:
                header = ''
                footer = obj_id
            else:
                header = obj_id[:i]
                footer = obj_id[i + 1:]

            # Convert footer to integer so more_itertools can process it
            try:
                footer = int(footer)
            except ValueError:
                err_msg = (
                    'Expected footer of obj ID for obj "{}" in {} to be '
                    'int-compatible. Instead received {}.'
                    ''.format(str(obj), parent_obj.__class__, footer))
                raise ValueError(err_msg)

            # Add to the dictionary of unique headers
            try:
                unique_headers[header].append(footer)
            except KeyError:
                unique_headers[header] = [footer]

        for header, footer_list in unique_headers.items():
            # Process header+delimiter
            if header == '':
                header_delim = header
            else:
                header_delim = '{}{}'.format(header, delimiter)
            # Sort the list and separate them into consecutive groups
            footer_list.sort()
            footer_ranges = mit.consecutive_groups(footer_list)
            # Write the appropriate ranges
            for footer_range in footer_ranges:
                footer_range = list(footer_range)
                if len(footer_range) == 1:
                    CTI_range = '"{}{:04d}", '.format(header_delim,
                                                      footer_range[0])
                else:
                    CTI_range = ('"{0}{1:04d} to {0}{2:04d}", '
                                 ''.format(header_delim, footer_range[0],
                                           footer_range[-1]))
                CTI_out += CTI_range
        CTI_out = '{}]'.format(CTI_out[:-2])
    return CTI_out
