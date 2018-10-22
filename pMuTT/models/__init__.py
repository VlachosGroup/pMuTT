# -*- coding: utf-8 -*-
"""
pMuTT.models
Vlachos group code to models.
Created on Tues Jul 10 12:40:00 2018
"""

def pMuTT_list_to_dict(pMuTT_list, key='name'):
    """Converts a pMuTT list to a dictionary using a specified attribute. This 
    allows for quicker searching.
    
    Parameters
    ----------
        pMuTT_list : list of objects
            List of pMuTT objects to convert
        key : str
            Name of attribute used as the keys for the dictionary
    Returns
    -------
        pMuTT_dict : dict
            Dictionary of pMuTT objects
    """
    return {getattr(obj, key): obj for obj in pMuTT_list}