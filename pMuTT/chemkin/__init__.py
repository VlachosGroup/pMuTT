# -*- coding: utf-8 -*-
from pMuTT.io_.jsonio import remove_class


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

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        return {'class': str(self.__class__),
                'name': self.name,
                'site_density': self.site_density,
                'density': self.density,
                'bulk_specie': self.bulk_specie, }

    @classmethod
    def from_dict(cls, json_obj):
        """Recreate an object from the JSON representation.

        Parameters
        ----------
            json_obj : dict
                JSON representation
        Returns
        -------
            BaseThermo : BaseThermo object
        """
        json_obj = remove_class(json_obj)
        return cls(**json_obj)
