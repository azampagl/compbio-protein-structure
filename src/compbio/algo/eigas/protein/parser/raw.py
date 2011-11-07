"""
Protein parser based on raw data.

@author Aaron Zampaglione <azampagl@azampagl.com>
@package EIGAs
@copyright 2011 Aaron Zampaglione
@license MIT
"""
from core import ProteinParser

class RAWProteinParser(ProteinParser):
    
    def __init__(self, name, data):
        """
        Init.
        
        Key arguments:
        name -- name of the protein.
        data -- data dictionary containing the protein data.
        """
        super(self.__class__, self).__init__(name)
        
        self._coords = data['coords']