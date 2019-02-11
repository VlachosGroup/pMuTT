# -*- coding: utf-8 -*-

class CatSite:
    """Catalyst site for Chemkin

    Attributes
    ----------
        name : str
            Name of the catalyst site
        site_density : float
            Catalyst site density in mol/cm2
        density : float
            Catalyst density in g/cm3
        bulk_specie : str
            Name of the bulk specie
    """

    def __init__(self, name, site_density, density, bulk_specie):
        self.name = name
        self.site_density = site_density
        self.density = density
        self.bulk_specie = bulk_specie