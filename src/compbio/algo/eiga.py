'''
Python implementation of the EIGA algorithm.

@author Aaron Zampaglione <azampagl@azapagl.com>
@copyright MIT
'''
from Bio.PDB.PDBParser import PDBParser
from scipy.spatial.distance import cdist

class Eiga(object):
    """
    EIGA.
    """
    
    # Distance matrix
    d = {}
    
    def __init__(self, structures):
        """
        Initializes the EIGA algorithm.
        
        Key arguments:
        structures -- structure tuple that contains the structure name and the .ent
                      file location.
        """
        for structure,file_name in structures:
            # Find the atomic coordinates for this structure.
            coords = self.parse_pdb(structure, file_name)
            # Calculate the distance matrices.
            self.d[structure] = cdist(coords, coords)
    
    def parse_pdb(self, structure, file_name):
        """
        Parses a pdb-format file and returns a list of coordinates.
        
        Key arguments:
        structure -- the structure's id 
        file_name -- the location of the .ent file.
        """
        # Build parse object.
        parser = PDBParser()
        
        structure = parser.get_structure(structure, file_name)
        
        # Unsure how to process multiple models.
        if len(structure.get_list()) > 1:
            raise
        
        model = structure.get_list()[0]
        
        # Unsure how to process multiple models.
        if len(model.get_list()) > 1:
            raise
        
        # The chain contains a list of residues.
        chain = model.get_list()[0]
        
        # Keep track of the alpha carbon coordinates.
        coords = []
        
        for residue in chain:
            for atom in residue:
                # We're only looking at the primary carbon atom.
                if atom.get_name() == "C":
                    coord = atom.get_coord()
                    coords.append((coord[0], coord[1], coord[2]))
        
        return coords
        