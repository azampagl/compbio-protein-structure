'''
Python implementation of the EIGA algorithm.

A Spectral Approach to Protein Structure Alignment
Yosi Shibberu and Allen Holder

@author Aaron Zampaglione <azampagl@azapagl.com>
@copyright MIT
'''
from Bio.PDB.PDBParser import PDBParser
from math import sqrt
from numpy import float, zeros
from numpy.linalg import eig
from scipy.spatial.distance import cdist

class Eiga(object):
    """
    EIGA.
    """
    
    class Protein(object):
        """
        Protein
        """
        
        name = ""
        
        #
        fingerprint = []
        
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
            
            # Unsure how to process multiple structures.
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
        
        def __init__(self, structure, file_name, opts=None):
            """
            Initializes a protein.
            
            Key arguments:
            structures -- structure tuple that contains the structure name and the .ent
                          file location.
            """
            # Set the protein name.
            self.name = structure
            
            # Find the atomic coordinates from the pdb-format file.
            coords = Eiga.Protein.coords(structure, file_name)
            
            # Calculate distance matrix
            dmatrix = cdist(coords, coords)
            
            # Create the contact matrix.
            cmatrix = self.cmatrix(dmatrix, opts['k'])
            
            # Find the eigvalues and eigvectors of the contact matrix.
            eigvalues, eigvectors = eig(cmatrix)
            
            # For each residue, we want to assign the "best"
            #  eigenvalue.
            #
            # @see section 3.1 of the report
            for j in range(len(coords)):
                
                # Init our max value and max angle
                max_value = None
                max_angle = None
                
                for i in range(len(eigvalues)):
                    angle = sqrt(eigvalues[i]) * abs(eigvectors[i][j])
                    if angle > max_angle:
                        max_value = eigvalues[i]
                        max_angle = angle
                
                self.fingerprint.append(max_value)
            
            print(self.fingerprint)
        
        def cmatrix(self, dmatrix, k):
            """
            Creates a contact matrix.
            
            Key arguments:
            dmatrix -- distance matrix
            k       -- cutoff threshold
            """
            l = len(dmatrix)
            rows = cols = range(l)
            
            ik = 1 / k
            
            cmatrix = zeros((l, l), dtype=float)
            
            for i in rows:
                for j in cols:
                    value = dmatrix[i][j]
                    if value >= 0.0 and value <= k:
                        cmatrix[i][j] = 1 - ik * value
            
            return cmatrix