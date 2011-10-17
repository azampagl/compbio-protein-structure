'''
Python implementation of the EIGA algorithm.

A Spectral Approach to Protein Structure Alignment
Yosi Shibberu and Allen Holder

@author Aaron Zampaglione <azampagl@azapagl.com>
@copyright MIT
'''
from Bio.PDB.PDBParser import PDBParser
from numpy import array, arange, diag, dot, empty
from numpy import float, inner, transpose, vectorize, zeros
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
        
        class __Node():
            """
            Node of the DP matrix.
            """
            # Backpointer to the previous node.
            prev = None
            
            # Score at the current node.
            score = None
            
            # Value at this node.
            value = None
            
            def __init__(self, value):
                """
                Initializes the value during initialization.
                
                Key arguments:
                value -- the value for this node.
                """
                self.value = value
        
        fingerprint1 = protein1.fingerprint
        fingerprint2 = protein2.fingerprint
        
        rows = len(fingerprint1)
        cols = len(fingerprint2)
        
        #Build a matrix full of Prefix objects
        vnode = vectorize(__Node)
        matrix = empty((rows, cols), dtype=object)
        
        # Initialize the matrix.
        for i in range(rows):
            for j in range(cols):
                matrix[i][j] = vnode(abs(fingerprint1[i] - fingerprint2[j]))
        
        # Initialize the scroe for the first row and column.
        for i in range(rows):
            matrix[i][0].score = matrix[i][0].value
        for j in range(cols):
            matrix[0][j].score = matrix[0][j].value
            
        # Determine score using DP.
        for i in range(1, rows):
            for j in range(1, cols):
                # Current node.
                current = matrix[i][j]
                # Find the previous top, diag, and left scores
                top = matrix[i - 1][j]
                diag = matrix[i - 1][j - 1]
                left = matrix[i][j - 1]
        
        # Preprocess the matrix
        pprint.pprint(matrix)
        
    class Protein(object):
        """
        Protein
        """
        
        # Name of the protein.
        name = ""
        
        # Fingerprint of the protein
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
            # Initialize fingerprint.
            self.fingerprint = []
            
            # Set the protein name.
            self.name = structure
            
            # Find the atomic coordinates from the pdb-format file.
            coords = Eiga.Protein.coords(structure, file_name)
            
            # Calculate distance matrix
            dmatrix = cdist(coords, coords)
            
            # Create the contact matrix.
            cmatrix = self.cmatrix(dmatrix)
            
            # Test definition 1
            #ei = zeros((1, len(cmatrix[0:])), dtype=float)[0]
            #ei[0] = 1.0
            #ej = zeros((1, len(cmatrix[0:])), dtype=float)[0]
            #ej[1] = 1.0
            #print(inner(ei, inner(cmatrix, ej)) == cmatrix[0][1])
            
            # Find the eigvalues and eigvectors of the contact matrix.
            eigvectors, eigvalues, eigvectorsT = svd(cmatrix)
            
            # Check SVD decomposition
            #c = dot(eigvectors, dot(diag(eigvalues), eigvectorsT))
            #print(abs(cmatrix - c) < (1 ** -15))
            
            # Calculate r
            r = dot((diag(eigvalues) ** 0.5), eigvectorsT)
            
            # Test proof 2
            #i = 0
            #j = 0
            #ei = zeros((1, len(cmatrix[0:])), dtype=float)[0]
            #ei[i] = 1.0
            #ej = zeros((1, len(cmatrix[0:])), dtype=float)[0]
            #ej[j] = 1.0
            #ri = inner(r, ei).transpose()
            #rj = inner(r, ej)
            #print(abs(inner(ri, rj) - cmatrix[i][j]) < (1 ** -15))
                       
            # For each residue, we want to assign the "best"
            #  eigenvalue.
            #
            # @see section 3.1 of the report
            for j in range(len(coords)):
                
                max_index = None
                max_value = None
                
                for i in range(len(r[:][0])):
                    if (r[i][j] > max_value):
                        max_index = i
                        max_value = r[i][j]
                
                self.fingerprint.append(eigvalues[max_index])
            
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