class Units:
    """Expresses units as Cantera CTI file.
    
    Attributes
    ----------
        length : str, optional
            Length unit. Default is 'cm'.
        time : str, optional
            Time unit. Default is 's'
        quantity : str, optional
            Quantity unit. Default is 'molec'
        energy : str, optional
            Energy unit. Default is 'cal/mol'
        act_energy : str, optional
            Activation energy unit. Default is 'cal/mol'
        pressure : str, optional
            Pressure unit. Default is 'bar'
        mass : str, optional
            Mass unit. Default is 'kg'
    """
    def __init__(self, length='cm', time='s', quantity='molec',
                 energy='cal/mol', act_energy='cal/mol', pressure='bar',
                 mass='kg'):
        self.length = length
        self.time = time
        self.quantity = quantity
        self.energy = energy
        self.act_energy = act_energy
        self.pressure = pressure
        self.mass = mass

    def to_CTI(self):
        """Writes the object in Cantera's CTI format.

        Parameters
        ----------
            max_line_len : int, optional
                Maximum number of characters in the line. Default is 80.
        Returns
        -------
            CTI_str : str
                Object represented as a CTI string.
        """
        cti_str = ('units(length="{}", time="{}", quantity="{}", energy="{}",\n'
                   '      act_energy="{}", pressure="{}", mass="{}")').format(
                       self.length, self.time, self.quantity, self.energy,
                       self.act_energy, self.pressure, self.mass)
        return cti_str

    def to_CTI_dict(self):
        """Returns a useful dictionary for CTI IO functions.
        
        Returns
        -------
            CTI_dict : dict
                Dictionary whose keys are the parameter names of CTI IO
                functions.
        """
        return {'length_unit': self.length,
                'time_unit': self.time,
                'quantity_unit': self.quantity,
                'energy_unit': self.energy,
                'act_energy_unit': self.act_energy,
                'pressure_unit': self.pressure,
                'mass_unit': self.mass}