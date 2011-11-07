"""
Base class for parsing protein data.

Currently only supports alpha carbon atom coordinates.

@author Aaron Zampaglione <azampagl@azampagl.com>
@package EIGAs
@copyright 2011 Aaron Zampaglione
@license MIT
"""
from ...exception import EIGAsException

from abc import ABCMeta

class ProteinParser(object):
    __metaclass__ = ABCMeta
    
    # Different implementations of the protein parser.
    TYPES = {'pdb': ('compbio.algo.eigas.protein.parser.pdb', 'PDBProteinParser'),
             'raw': ('compbio.algo.eigas.protein.parser.raw', 'RAWProteinParser'),
             'txt': ('compbio.algo.eigas.protein.parser.txt', 'TXTProteinParser'),
             }
    
    @classmethod
    def factory(cls, name, args):
        """
        Returns a specific protein parser implementation.
        
        Key arguments:
        name -- the type of parser (pdb, raw, txt).
        args -- the arguments to pass to the new object.
        """
        if name in ProteinParser.TYPES:
            meta = ProteinParser.TYPES[name]
            mod = __import__(meta[0], fromlist=[meta[1]])
            kls = getattr(mod, meta[1])
            return kls(*args)
        
        raise EIGAsException('Protein parser type not supported: ' + name)
        
    def __init__(self, name):
        """
        Init.
        
        Key arguments:
        name -- the name of the protein
        """
        self._name = name
    
    def coords(self):
        """
        Return the alpha carbon atom coordinates.
        """
        return self._coords
    
    def name(self):
        """
        Return the name of the protein.
        """
        return self._name