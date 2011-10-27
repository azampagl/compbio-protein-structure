"""
Python implementation of the EIGAs algorithm.

A Spectral Approach to Protein Structure Alignment
Yosi Shibberu and Allen Holder

@author Aaron Zampaglione <azampagl@azampagl.com>
@package EIGAs
@copyright 2011 (c) Aaron Zampaglione
@license MIT
"""
from Bio.PDB.PDBParser import PDBParser
from numpy import diag, dot, empty, float, zeros
from numpy.linalg import svd
from scipy.spatial.distance import cdist

from exception import EIGAsException

class EIGAs(object):
    """
    EIGAs.
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
        
        class __Node(object):
            """
            Node of the DP matrix.
            """
            # Backpointer.
            prev = None
            
            # Score at the current node.
            score = 0.0
            
            # Indices of protein 1
            indices1 = None
            
            # Indices of protein 2
            indices2 = None
            
            def __init__(self):
                """
                Initializes the node.
                """
                self.prev = None
                self.score = 0.0
                self.indices1 = None
                self.indices2 = None
        
        # Quick reference the protein fingerprints.
        fingerprint1 = protein1.fingerprint
        fingerprint2 = protein2.fingerprint
        
        rows = len(fingerprint1)
        cols = len(fingerprint2)
        
        #Build a matrix full of Prefix objects
        matrix = empty((rows, cols), dtype=object)
        # Initialize the matrix.
        for i in range(rows):
            for j in range(cols):
                matrix[i][j] = __Node()
            
        # Determine score using DP.
        for i in range(1, rows):
            for j in range(1, cols):                
                # Find the previous top, diag, and left components.
                top = matrix[i - 1][j]
                diag = matrix[i - 1][j - 1]
                left = matrix[i][j - 1]
                
                # Find the scores.
                top_score = top.score + abs(fingerprint1[i - 1] - fingerprint2[j]) + 1.0
                diag_score = diag.score + abs(fingerprint1[i - 1] - fingerprint2[j - 1])
                left_score = left.score + abs(fingerprint1[i] - fingerprint2[j - 1]) + 1.0
                
                # Top
                if (top_score <= diag_score and top_score <= left_score):
                    matrix[i][j].prev = top
                    matrix[i][j].score = top_score
                    matrix[i][j].indices1 = i - 1
                # Diagonal
                elif (diag_score <= top_score and diag_score <= left_score):
                    matrix[i][j].prev = diag
                    matrix[i][j].score = diag_score
                    matrix[i][j].indices1 = i - 1
                    matrix[i][j].indices2 = j - 1
                # Left
                else:
                    matrix[i][j].prev = left
                    matrix[i][j].score = left_score
                    matrix[i][j].indices2 = j - 1
        
        # Follow the pointers backwards to rebuild the globally aligned sequences.
        s1 = []
        s2 = []
        node = matrix[-1][-1]
        while node.prev != None:
            s1.insert(0, node.indices1)
            s2.insert(0, node.indices2)
            node = node.prev
        
        return matrix[-1][-1].score, s1, s2
    
    @staticmethod
    def aligned(protein1=None, protein2=None, seq1=None, seq2=None):
        """
        Counts the number of aligned items.
        
        Either both proteins OR both sequences must be provided.
        
        Key arguments:
        protein1 -- first protein. [optional]
        protein2 -- second protein. [optional]
        seq1     -- first sequence. [optional]
        seq2     -- second sequence. [optional]
        """
        aligned = 0
        
        # Align the proteins if the sequences weren't provided.
        if not seq1 or not seq2:
            _, seq1, seq2 = EIGAs.align(protein1, protein2)
        
        for i in range(len(seq1)):
            if seq1[i] and seq2[i]:
                aligned += 1
        
        return aligned
    
    class Protein(object):
        """
        Protein
        """
        
        # Name of the protein.
        name = ""
        
        # Fingerprint of the protein
        fingerprint = []
        
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
            coords = self.coords(file_name)
            
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
            _, eigvalues, eigvectorsT = svd(cmatrix)
            
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
            
            ik = 1 / EIGAs.cutoff
            
            cmatrix = zeros((l, l), dtype=float)
            
            for i in rows:
                for j in cols:
                    value = dmatrix[i][j]
                    if value >= 0.0 and value <= EIGAs.cutoff:
                        cmatrix[i][j] = 1 - ik * value
            
            return cmatrix
        
        def coords(self, file_name):
            """
            Parses a pdb-format file and returns a list of coordinates.
            
            Key arguments:
            file_name -- the location of the .ent file.
            """
            # Build parse object.
            parser = PDBParser()
            
            structure = parser.get_structure(self.name, file_name)
            
            # Unsure how to process multiple structures.
            if len(structure.get_list()) != 1:
                raise EIGAsException("Error reading structure for '" + self.name + '"')
            
            model = structure.get_list()[0]
            
            # Unsure how to process multiple models.
            if len(model.get_list()) != 1:
                raise EIGAsException("Error reading model for '" + self.name + '"')
            
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