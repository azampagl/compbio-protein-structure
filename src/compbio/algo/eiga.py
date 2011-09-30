'''
Python implementation of the EIGA algorithm.

A Spectral Approach to Protein Structure Alignment
Yosi Shibberu and Allen Holder

@author Aaron Zampaglione <azampagl@azapagl.com>
@copyright MIT
'''
from Bio.PDB.PDBParser import PDBParser
from math import sqrt
from numpy import diag, float, reshape, transpose, zeros
from numpy.linalg import det, eig, svd
from scipy.spatial.distance import cdist

import pprint
import sys

class Eiga(object):
    """
    EIGA.
    """
    
    # Contact matrix cutoff threshold.
    cutoff = 8.0
    
    @staticmethod
    def align(protein1, protein2):
        """
        Calculates the score between two proteins.
        
        Key arguments:
        protein1 -- the first protein.
        protein2 -- the second protein.
        """
        fingerprint1 = protein1.fingerprint
        fingerprint2 = protein2.fingerprint
        
        rows = len(fingerprint1)
        cols = len(fingerprint2)
        
        matrix = zeros((rows, cols), dtype=float)
        
        print(fingerprint1)
        print(fingerprint2)
        # Compute the score matrix of the vectors.
        for i in range(rows):
            for j in range(cols):
                #print(abs(fingerprint1[i] - fingerprint2[j]))
                matrix[i][j] = abs(fingerprint1[i] - fingerprint2[j])
        
        # Preprocess the matrix
        #pprint.pprint(matrix)
        
    
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
        
        def __init__(self, structure, file_name):
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
            cmatrix = self.cmatrix(dmatrix)
            
            #pprint.pprint(cmatrix)
            #sys.exit(0)
            
            # Find the eigvalues and eigvectors of the contact matrix.
            eigvalues, eigvectors = eig(cmatrix)
            eigvectors = reshape(eigvectors, order='F')
            eigvalues = sorted(eigvalues, reverse=True)
            print(eigvalues)
            sys.exit(0)
            # Find our intrinsic contact coords.
            r = (diag(eigvalues) ** 0.5) * transpose(eigvectors)
            #print(diag(eigvalues))
            #print("")
            #print(eigvectors)
            
            #print(det(cmatrix))
            #print("")
            #print("")
            #print(eigvectors * transpose(eigvectors))
            #print(eigvectors * diag(eigvalues) * transpose(eigvectors))
            #print(det(eigvectors) * det(diag(eigvalues)) * det(transpose(eigvectors)))
            #r = sqrt(eigvalues)
            #print(r)
            sys.exit(0)
            #print(eigvalues)
            
            #u, s, _ = svd(cmatrix)
            
            #pprint.pprint(eigvalues)
            #pprint.pprint(diag(u))
            #sys.exit(0)
            
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
            print("")
            print("")
            #print(self.fingerprint)
        
        def cmatrix(self, dmatrix):
            """
            Creates a contact matrix.
            
            Key arguments:
            dmatrix -- distance matrix
            """
            l = len(dmatrix)
            rows = cols = range(l)
            
            ik = 1 / Eiga.cutoff
            
            cmatrix = zeros((l, l), dtype=float)
            
            for i in rows:
                for j in cols:
                    value = dmatrix[i][j]
                    if value >= 0.0 and value <= Eiga.cutoff:
                        cmatrix[i][j] = 1 - ik * value
            
            return cmatrix