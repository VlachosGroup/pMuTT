class Units:
    def __init__(self, length='cm', time='s', quantity='molecule',
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
        cti_str = ('units(length="{}", time="{}", quantity="{}", energy="{}",\n'
                   '      act_energy="{}", pressure="{}", mass="{}")').format(
                       self.length, self.time, self.quantity, self.energy,
                       self.act_energy, self.pressure, self.mass)
        return cti_str
