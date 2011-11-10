"""
Python implementation of the EIGAs algorithm.

A Spectral Approach to Protein Structure Alignment
Yosi Shibberu and Allen Holder

@author Aaron Zampaglione <azampagl@azampagl.com>
@package EIGAs
@copyright 2011 Aaron Zampaglione
@license MIT
"""
from numpy import empty

class EIGAs(object):
    
    # Contact matrix cutoff threshold.
    cutoff = 8.0
    
    class _Node(object):
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
            _, seq1, seq2 = EIGAs.global_align(protein1, protein2)
        
        for i in range(len(seq1)):
            if seq1[i] != None and seq2[i] != None:
                aligned += 1
        
        return aligned
    
    @classmethod
    def global_align(cls, protein1, protein2):
        """
        Globally aligns two proteins.
        
        Key arguments:
        protein1 -- the first protein.
        protein2 -- the second protein.
        """
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
                matrix[i][j] = EIGAs._Node()
        
        # Init first row.
        for i in range(rows):
            matrix[i][0].score = float(i)#abs(fingerprint1[i] - fingerprint2[0])
            matrix[i][0].indices1 = i
            matrix[i][0].indices2 = 0
    
        # Init first col.
        for i in range(cols):
            matrix[0][i].score = float(i)#abs(fingerprint1[0] - fingerprint2[i])
            matrix[0][i].indices1 = 0
            matrix[0][i].indices2 = i
            
        # Determine score using DP.
        for i in range(1, rows):
            for j in range(1, cols):                
                # Find the previous top, diag, and left components.
                top = matrix[i - 1][j]
                diag = matrix[i - 1][j - 1]
                left = matrix[i][j - 1]
                
                # Find the scores.
                top_score = top.score + 1.0
                diag_score = diag.score + abs(fingerprint1[i] - fingerprint2[j])
                left_score = left.score + 1.0
                
                # Top
                if (top_score <= diag_score and top_score <= left_score):
                    matrix[i][j].prev = top
                    matrix[i][j].score = top_score
                    matrix[i][j].indices1 = i
                # Diagonal
                elif (diag_score <= top_score and diag_score <= left_score):
                    matrix[i][j].prev = diag
                    matrix[i][j].score = diag_score
                    matrix[i][j].indices1 = i
                    matrix[i][j].indices2 = j
                # Left
                else:
                    matrix[i][j].prev = left
                    matrix[i][j].score = left_score
                    matrix[i][j].indices2 = j
        
        # Follow the pointers backwards to rebuild the globally aligned sequences.
        s1 = []
        s2 = []
        node = matrix[-1][-1]
        while node != None:
            s1.insert(0, node.indices1)
            s2.insert(0, node.indices2)
            node = node.prev
        
        return matrix[-1][-1].score, s1, s2
    
    @classmethod
    def local_align(cls, protein1, protein2):
        """
        Finds the best local alignment between two proteins.
        
        Key arguments:
        protein1 -- the first protein.
        protein2 -- the second protein.
        """
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
                matrix[i][j] = EIGAs._Node()
        
        # Init first row.
        for i in range(rows):
            matrix[i][0].score = max(fingerprint1[i], fingerprint2[0])
            matrix[i][0].indices1 = i
            matrix[i][0].indices2 = 0
    
        # Init first col.
        for i in range(cols):
            matrix[0][i].score = max(fingerprint1[0], fingerprint2[i])
            matrix[0][i].indices1 = 0
            matrix[0][i].indices2 = i
        
        optimal = None
        
        # Determine score using DP.
        for i in range(1, rows):
            for j in range(1, cols):                
                # Find the previous top, diag, and left components.
                top = matrix[i - 1][j]
                diag = matrix[i - 1][j - 1]
                left = matrix[i][j - 1]
                
                # Find the scores.
                score = abs(fingerprint1[i] - fingerprint2[j])
                top_score = top.score + score + 1.0
                diag_score = diag.score + score
                left_score = left.score + score + 1.0
                
                # Top
                if (top_score <= diag_score and top_score <= left_score):
                    matrix[i][j].prev = top
                    matrix[i][j].score = top_score
                    matrix[i][j].indices1 = i
                # Diagonal
                elif (diag_score <= top_score and diag_score <= left_score):
                    matrix[i][j].prev = diag
                    matrix[i][j].score = diag_score
                    matrix[i][j].indices1 = i
                    matrix[i][j].indices2 = j
                # Left
                else:
                    matrix[i][j].prev = left
                    matrix[i][j].score = left_score
                    matrix[i][j].indices2 = j
                
                #if matrix[i][j].score < optimal.score:
                #    optimal = matrix[i][j]
        
        # Follow the pointers backwards to rebuild the globally aligned sequences.
        s1 = []
        s2 = []
        node = optimal
        while node != None:
            s1.insert(0, node.indices1)
            s2.insert(0, node.indices2)
            node = node.prev
        
        return optimal.score, s1, s2