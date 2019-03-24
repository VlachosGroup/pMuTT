# -*- coding: utf-8 -*-
"""
pMuTT.io.db

Read from/write to databases
"""

from pymongo import MongoClient
from pMuTT.io.json import json_to_pMuTT

class PMRester:
    """User-friendly REST interface with the database

    Attributes
    ----------
        client : `pymongo.mongo_client.MongoClient`_ object
        database : `pymongo.database.Database`_ object
        collection : `pymongo.collection.Collection`_ object

    .. _`pymongo.mongo_client.MongoClient`: http://api.mongodb.com/python/current/api/pymongo/mongo_client.html#pymongo.mongo_client.MongoClient
    .. _`pymongo.database.Database`: http://api.mongodb.com/python/current/api/pymongo/database.html#pymongo.database.Database
    .. _``pymongo.collection.Collection`: http://api.mongodb.com/python/current/api/pymongo/collection.html#pymongo.collection.Collection
    """
    
    def __init__(self, client, database, collection):
        self.client = client
        self.database = database
        self.collection = collection

    @classmethod
    def connect(cls, uri, database, collection):
        """Connect with the database

        Parameters
        ----------
            uri : str
                URI to connect to the database. See the
                `URI MongoDB documentation`_ for sample URIs
            database : str
                Name of the database to access from the `uri`
            collection : str
                Name of the collection within `database` to access
        Returns
        -------
            Nasa : Nasa object

        .. _`URI MongoDB documentation`: https://docs.mongodb.com/manual/reference/connection-string/#connections-connection-examples
        """
        client = MongoClient(uri)
        database = client[database]
        collection = database[collection]
        return cls(client=client, database=database, collection=collection)

    def get_specie_from_smiles(self, smiles, obj_type=None, max_records=1):
        """Searches database for object matching SMILES string

        Parameters
        ----------
            smiles : str
                Smiles representation of interested species
            obj_type : str, optional
                Type of object to return. Default will return the first
                encountered specie matching the smiles string regardless of
                type. Suggested types are:
                
                - 'nasa'
                - 'shomate'
                - 'statmech'

            max_records : int, optional
                Maximum number of records to find. Set this parameter to a
                -1 if all records are desired. Default is 1
        Returns
        -------
            specie : pMuTT model object or list of pMuTT model objects
                Object(s) representing the interested species
        """
        # Set up filter to find records
        db_filter = {'smiles': smiles}
        if obj_type is not None:
            db_filter['type'] = obj_type

        # Get species from database
        species = []
        for i, result in enumerate(self.collection.find(filter=db_filter)):
            specie = json_to_pMuTT(result)
            # Only return one record if desired
            if max_records == 1:
                species = specie
                break
            species.append(specie)
            # Determine if record limit reached
            if (i+1) == max_records:
                break
        return species

    def get_species_from_smiles(self, species_smiles, obj_type=None):
        """Searches database for object matching SMILES string

        Parameters
        ----------
            species_smiles : list of str
                Smiles representation of interested species
            obj_type : str, optional
                Type of object to return. Default will return the first
                encountered specie matching the smiles string regardless of
                type. Suggested types are:
                
                - 'nasa'
                - 'shomate'
                - 'statmech'

        Returns
        -------
            species : list of pMuTT model objects
                Objects representing the interested species
        """
        species = []
        for smiles in species_smiles:
            species.append(self.get_specie_from_smiles(smiles=smiles, 
                                                       obj_type=obj_type,
                                                       max_records=1))
        return species

    def find(self, **kwargs):
        """Find records in collection using desired criteria

        Parameters
        ----------
            max_records : int
                Maximum number of records to return
            kwargs : keyword arguments
                Arguments to use for `self.collection.find` method. See 
                `pymongo.collection.Collection.find`_ method
        Returns
        -------
            results : list of pMuTT objects
        .. _`pymongo.collection.Collection.find`: http://api.mongodb.com/python/current/api/pymongo/collection.html#pymongo.collection.Collection.find     
        """
        return [json_to_pMuTT(result)
                for result in self.collection.find(**kwargs)]

    def submit(self, pMuTT_obj):
        """Submit pMuTT object to `collection` encoded as JSON

        Parameters
        ----------
            pMuTT_obj : `pMuTT` object
                pMuTT object that can be encoded using `to_dict` method
        """
        self.collection.insert_one(pMuTT_obj.to_dict())