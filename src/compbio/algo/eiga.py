'''
Python implementation of the EIGA algorithm.

A Spectral Approach to Protein Structure Alignment
Yosi Shibberu and Allen Holder

@author Aaron Zampaglione <azampagl@azapagl.com>
@copyright MIT
'''
from Bio.PDB.PDBParser import PDBParser
from numpy import zeros, float
from scipy.spatial.distance import cdist

class Eiga(object):
    """
    EIGA.
    """
    
    @staticmethod
    def coords(structure, file_name):
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
    
    @staticmethod
    def dmatrix(coords):
        """
        Returns the distance matrix for a list of atomic coordinates.
        
        Key arguments:
        coords -- list of coordinates.
        """
        return cdist(coords, coords)
    
    # Options.
    opts = {
            # Cutoff parameter for the contact function.
            'k': 1.0,
            }
    
    # Distance matrices
    dmatrices = {}
    
    # Contact matrices
    cmatrices = {}
    
    def __init__(self, structures, opts=None):
        """
        Initializes the EIGA algorithm.
        
        Key arguments:
        structures -- structure tuple that contains the structure name and the .ent
                      file location.
        """
        # Set the options.
        if opts:
            for opt, value in opts:
                self.opts[opt] = value
        
        for structure,file_name in structures:
            # Find the atomic coordinates for this structure.
            coords = Eiga.coords(structure, file_name)
            # Calculate the distance matrices.
            dmatrix = Eiga.dmatrix(coords)
            self.dmatrices[structure] = dmatrix
            self.cmatrices[structure] = self.cmatrix(dmatrix)
    
    def cmatrix(self, dmatrix):
        """
        Returns the contact matrix.
        
        Key arguments:
        dmatrix -- distance matrix
        """
        l = len(dmatrix)
        rows = cols = range(l)
        
        # Find our cutoff parameter and it's inverse.
        k = self.opts['k']
        ik = 1 / k
        
        cmatrix = zeros((l, l), dtype=float)
        
        for i in rows:
            for j in cols:
                value = dmatrix[i][j]
                if value >= 0.0 and value <= k:
                    cmatrix[i][j] = 1 - ik * value
        
        return cmatrix