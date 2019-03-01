# -*- coding: utf-8 -*-
"""
pMuTT.io_.db

Read from/write to databases
"""
from pymongo import MongoClient

def get_collection(uri, database, collection):
    """Accesses a MongoDB and returns a specified collection

    Parameters
    ----------
        uri : str
            URI to connect to the database. See the
            `URI MongoDB documentation`_ for sample URIs
        database : str
            Name of the database to access
        collection : str
            Name of the collection to access
    Returns
    -------
        db_collection : `pymongo.collection.Collection`_ object
            Collection object where the user can insert or query records

    .. _`URI MongoDB documentation`: https://docs.mongodb.com/manual/reference/connection-string/#connections-connection-examples
    .. _`pymongo.collection.Collection`: http://api.mongodb.com/python/current/api/pymongo/collection.html#pymongo.collection.Collection
    """
    client = MongoClient(uri)
    db = client[database]
    return db[collection]