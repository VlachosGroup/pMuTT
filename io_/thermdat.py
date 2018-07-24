def write_thermdat(filename, nasa_species):
    """
    Writes thermdats in the Chemkin format

    Parameters
        filename - str
            Output file name
        nasa_species - list 
            List of Thermochemistry.models.empirical.nasa.Nasa objects
    """
    with open(filename, 'w') as f_ptr:
        f_ptr.write('THERMO ALL\n       100       500      1500\n')
        
        float_string = '%.8E'    
        for nasa_specie in nasa_species:
            _write_line1(f_ptr, nasa_specie)
            _write_line2(f_ptr, nasa_specie, float_string)
            _write_line3(f_ptr, nasa_specie, float_string)
            _write_line4(f_ptr, nasa_specie, float_string)
        f_ptr.write('END')    
        
def _write_line1(thermdat_file, nasa_specie):
    """
    Writes the first line of the thermdat file, which contains information on the composition, phase, and temperature ranges

    Parameters
        thermdat_file - file object
            Thermdat file that is being written to
        nasa_specie - Thermochemistry.models.empirical.thermdat.Thermdat object
    """
    element_pos = [
    	24, #Element 1
    	28, #Element 1#
    	29, #Element 2
    	33, #Element 2#
    	34, #Element 3
    	38, #Element 3#
    	39, #Element 4
    	43] #Element 4#
    temperature_pos = [
        44, # Phase
        45, # T_low
        55, # T_high
        65, # T_mid
        79] # Line num

    #Adjusts the position based on the number of elements
    line1_pos = []
    for element, val in nasa_specie.elements.items():
    	if val > 0.:
    		line1_pos.append(element_pos.pop(0))
    		line1_pos.append(element_pos.pop(0))
    line1_pos.extend(temperature_pos)

    #Creating a list of the text to insert
    line1_fields = [nasa_specie.name]
    for element, val in nasa_specie.elements.items():
    	if val > 0.:
	        line1_fields.extend([element, '%d' % val])
    line1_fields.extend([nasa_specie.phase, '%.1f' % nasa_specie.T_low, '%.1f' % nasa_specie.T_high, '%.1f' % nasa_specie.T_mid])

    #Write the content with appropriate spacing
    line = ''
    for pos, field in zip(line1_pos, line1_fields):
        line += field
        line = _insert_space(pos, line)
    line += '1\n'
    thermdat_file.write(line)

def _write_line2(thermdat_file, nasa_specie, float_string):
    """
    Writes the second line of the thermdat file
    """
    line = ''
    for i in range(5):
        a = nasa_specie.a_high[i]
        if a >= 0:
            line += ' '            
        line += float_string % a
    line += '    2\n'
    thermdat_file.write(line)

def _write_line3(thermdat_file, nasa_specie, float_string):
    """
    Writes the third line of the thermdat file
    """
    line = ''
    for i in range(5):
        if i < 2:
            a = nasa_specie.a_high[i+5]
        else:
            a = nasa_specie.a_low[i-2]
        if a >= 0:
            line += ' '
        line += float_string % a
    line += '    3\n'
    thermdat_file.write(line)

def _write_line4(thermdat_file, nasa_specie, float_string):
    """
    Writes the fourth line of the thermdat file
    """
    line = ''
    for i in range(3,7):
        a = nasa_specie.a_low[i]
        if a >= 0:
            line += ' '
        line += float_string % a
    line += '                   4\n'
    thermdat_file.write(line)

def _insert_space(end_index, string):
    """
    Inserts the number of spaces required given the string and the position of
    the next non-blank field.
    """
    string += ' ' * (end_index - len(string))
    return string



def read_thermdat(filename):
    """
    Reads the thermdat file and returns an list of molecule objects with the
    symbols, temperature ranges and CHON fields.
    """
    raise NotImplementedError
    thermdats = []
    #Possible atoms 
    atoms = ['C', 'H', 'O', 'N']
    site_type = 1
    #Positions at which coefficients start in lines 2, 3, 4
    pos_index = [0, 15, 30, 45, 60]
    with open(thermdat_path, 'r') as thermdat_file:
        for line in thermdat_file:
            if len(line) > 2:
                #Finds the line number
                line_num = line[-3:]
                if '1' in line_num:
                    """
                    Encountered new species
					"""
                    #Parsing name
                    name_index = line.find(' ')
                    symbol = line[0:symbol_index]
                    buf= line[symbol_index:-1]

                    """
                    Parsing elements
                    """
                    CHON = []
                    for atom in atoms:
                        CHON.append(_get_CHON_value(buf, symbol, atom, verbose, warn))

                    """
                    Parsing phase
                    """
                    try:
                        phase = re.search("(G|S)", buf).group(0)
                    except AttributeError:
                        warnings.warn("Phase not found. Assuming surface phase.")
                        phase='S'

                    """
                    Parsing temperatures
                    """
                    T_range = _get_T_values(line, symbol)
                elif '2' in line_num:
                    a_high = []
                    for i in pos_index:
                        a_high.append(float(line[i:i+15]))
                elif '3' in line_num:
                    a_low  = []
                    for i in pos_index:
                        if len(a_high) < 7:
                            a_high.append(float(line[i:i+15]))
                        else:
                            a_low.append(float(line[i:i+15]))
                elif '4' in line_num:
                    for i in pos_index:
                        if len(a_low) < 7:
                            a_low.append(float(line[i:(i+15)]))
                    thermdats.append(
                    	Thermdat(
                    		name = name, 
                            elements = elements, 
                            phase = phase,
                            T_low = T_range[0],
                            T_mid = T_range[2],
                            T_high = T_range[1],
                            a_low = np.array(a_low),
                            a_high = np.array(a_high)))
    return thermdats

def _get_CHON_value(string, symbol, atom, verbose = True, warn = True):
    """Returns the number of C, H, O or N in an atom based on the string inputted."""
    #Looks for a pattern which starts with the atom type (e.g. C), has spaces
    #and then ends with a digit    
    pattern = '%s +\d+' % atom  
    try:
        buf = re.search(pattern, string).group(0)
    except AttributeError:
        if warn:
            warnings.warn("Unable to find %s atom in species %s. Returning 0." % (atom, symbol))
        return 0
    buf = re.search('\d+', buf).group(0)
    return int(float(buf))

def _get_T_values(string, symbol, verbose = True):
    """Returns T_low, T_mid and T_high given a string inputted."""
    #Looks for three floating numbers separated by a space
    pattern = "[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)? +[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)? +[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"
    try:
        buf = re.search(pattern, string).group(0)
    except AttributeError:
        if verbose:
            print(("Warning: Unable to find temperature limits in species %s. Returning 0s." % (symbol)))
        return [0, 0, 0]
    T_limits = re.split(" +", buf)
    T_out = []
    for T in T_limits:
        T_out.append(float(T))
    if len(T_out) < 3:
        if verbose:
            print(("Warning: Not all temperatures found for species %s. Returning 0s for unfound values." % symbol))
        while len(T_out) < 3:
            T_out.append(0)
    return T_out


def thermdat_to_json(filename, thermdats):
    pass

def json_to_thermdats(filename):
    pass