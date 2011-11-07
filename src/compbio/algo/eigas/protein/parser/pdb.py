"""
Parses a pdb file containing protein data.

@author Aaron Zampaglione <azampagl@azampagl.com>
@package EIGAs
@copyright 2011 Aaron Zampaglione
@license MIT
"""
from ...exception import EIGAsException
from core import ProteinParser

from Bio.PDB.PDBParser import PDBParser
from decimal import Decimal

class PDBProteinParser(ProteinParser):
    
    def __init__(self, name, file_name):
        """
        Init.
        
        Key arguments:
        name      -- name of the protein
        file_name -- the location of the pdb formatted file.
        """
        super(self.__class__, self).__init__(name)
        
        # Build parse object.
        parser = PDBParser(QUIET=True)
        
        structure = parser.get_structure(name, file_name)
        
        # Unsure how to process multiple structures.
        if len(structure.get_list()) < 1:
            raise EIGAsException("Error reading structure for '" + name + '"')
        
        model = structure.get_list()[0]
        
        # Unsure how to process multiple models.
        if len(model.get_list()) < 1:
            raise EIGAsException("Error reading model for '" + name + '"')
            
        # The chain contains a list of residues.
        chain = model.get_list()[0]
            
        # Keep track of the alpha carbon coordinates.
        self._coords = []
        
        for residue in chain:
            for atom in residue:
                # We're only looking at the primary carbon atom.
                if atom.get_name() == 'CA':
                    self._coords.append(map(Decimal, (map(str, atom.get_coord()))))