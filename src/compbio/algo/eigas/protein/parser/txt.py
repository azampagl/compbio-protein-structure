"""
Parses a text file containing protein data.

Currently only supports tab-delimited alpha carbon atom coordinates.

File format should be in the following:
        
        x     y     z
        
        23.0\t24.0\t25.0\n
        26.0\t27.0\t28.0\n
        ...

@author Aaron Zampaglione <azampagl@azampagl.com>
@package EIGAs
@copyright 2011 Aaron Zampaglione
@license MIT
"""
from core import ProteinParser

from decimal import Decimal

class TXTProteinParser(ProteinParser):
    
    def __init__(self, name, file_name):
        """
        Init.
        
        Key arguments:
        name      -- name of the protein
        file_name -- the location of the txt file.
        """
        super(self.__class__, self).__init__(name)
        
        self._coords = []
        
        for line in open(file_name, 'r'):
            self._coords.append(tuple(map(Decimal, line[:-1].split('\t'))))