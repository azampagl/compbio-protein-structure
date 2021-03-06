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
    
    # Gap penalty.
    GAP_PENALTY = 1
    
    class Node(object):
        """
        Node of the DP matrix.
        """
        
        def __init__(self, score=0):
            """
            Initializes the node.
            
            Key arguments:
            score -- initial score
            """
            self.prev = None
            self.score = score
            self.value = 0
            self.gaps = 0
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
                # Initialize matrix with difference in fingerprints.
                matrix[i][j] = EIGAs._Node(abs(fingerprint1[i] - fingerprint2[j]))
        
        # Init first row.
        for i in range(rows):
            matrix[i][0].value = float(i)
            matrix[i][0].indices1 = i
            matrix[i][0].indices2 = 0
    
        # Init first col.
        for j in range(cols):
            matrix[0][j].value = float(j)
            matrix[0][j].indices1 = 0
            matrix[0][j].indices2 = j
            
        # Determine score using DP.
        for i in range(1, rows):
            for j in range(1, cols):                
                # Find the previous top, diag, and left components.
                top = matrix[i - 1][j]
                diag = matrix[i - 1][j - 1]
                left = matrix[i][j - 1]
                
                # Find the scores.
                top_score = top.value + matrix[i - 1][j].score + EIGAs.GAP_PENALTY
                diag_score = diag.value + matrix[i - 1][j - 1].score
                left_score = left.value + matrix[i][j - 1].score + EIGAs.GAP_PENALTY
                
                # Top
                if (top_score <= diag_score and top_score <= left_score):
                    matrix[i][j].prev = top
                    matrix[i][j].value = top.value + EIGAs.GAP_PENALTY
                    matrix[i][j].indices1 = i
                # Diagonal
                elif (diag_score <= top_score and diag_score <= left_score):
                    matrix[i][j].prev = diag
                    matrix[i][j].value = diag.value
                    matrix[i][j].indices1 = i
                    matrix[i][j].indices2 = j
                # Left
                else:
                    matrix[i][j].prev = left
                    matrix[i][j].value = left.value + EIGAs.GAP_PENALTY
                    matrix[i][j].indices2 = j
        
        # Follow the pointers backwards to rebuild the globally aligned sequences.
        s1 = []
        s2 = []
        node = matrix[-1][-1]
        while node != None:
            s1.insert(0, node.indices1)
            s2.insert(0, node.indices2)
            node = node.prev
        
        # Fill in missing indices in the beginning.
        while s1[0]:
            s1.insert(0, s1[0] - 1)
            s2.insert(0, None)
        while s2[0]:
            s1.insert(0, None)
            s2.insert(0, s2[0] - 1)
                
        return matrix, s1, s2
    
    @classmethod
    def local_align(cls, protein1, protein2, max_gaps=0):
        """
        Local alignment.
        
        Key arguments:
        protein1  -- the first protein.
        protein2  -- the second protein.
        min_value -- only return alignments that meet a minimum value (think
                     of this as a threshold for the minimum number of structurally
                     aligned residues). [optional]
        """
        # Quick reference the protein fingerprints.
        fingerprint1 = protein1.fingerprint
        fingerprint2 = protein2.fingerprint
        
        rows = len(fingerprint1) + 1
        cols = len(fingerprint2) + 1
        
        #Build a matrix full of Prefix objects
        matrix = empty((rows, cols), dtype=object)
        # Initialize the matrix.
        for i in range(1, rows):
            for j in range(1, cols):
                # Initialize matrix with difference in fingerprints.
                matrix[i][j] = EIGAs.Node(abs(fingerprint1[i - 1] - fingerprint2[j - 1]))
        
        # Init first row.
        for i in range(1, rows):
            matrix[i][0] = EIGAs.Node()
            matrix[i][0].gaps = i
            matrix[i][0].indices1 = i
            matrix[i][0].indices2 = None
    
        # Init first col.
        for j in range(1, cols):
            matrix[0][j] = EIGAs.Node()
            matrix[0][j].gaps = j
            matrix[0][j].indices1 = None
            matrix[0][j].indices2 = j
        
        # Top left corner.
        matrix[0][0] = EIGAs.Node()
        
        optimals = []
            
        # Determine score using DP.
        for i in range(1, rows):
            for j in range(1, cols):                
                # Find the previous top, diag, and left components.
                top = matrix[i - 1][j]
                diag = matrix[i - 1][j - 1]
                left = matrix[i][j - 1]
                
                # See if the fingeprints at the current position 'match'.
                if matrix[i][j].score < EIGAs.GAP_PENALTY:
                    value = 2 * EIGAs.GAP_PENALTY
                else:
                    value = -1 * EIGAs.GAP_PENALTY
                
                # Find the values for each direction.
                top_score = top.value - EIGAs.GAP_PENALTY
                diag_score = diag.value + value
                left_score = left.value - EIGAs.GAP_PENALTY
                
                # Diagonal
                if (diag_score >= top_score and diag_score >= left_score and diag_score >= 0):
                    matrix[i][j].prev = diag
                    matrix[i][j].gaps = diag.gaps
                    matrix[i][j].value = diag_score
                    matrix[i][j].indices1 = i
                    matrix[i][j].indices2 = j 
                # Top
                elif (top_score >= diag_score and top_score >= left_score and top_score >= 0):
                    matrix[i][j].prev = top
                    matrix[i][j].gaps = top.gaps + 1
                    matrix[i][j].value = top_score
                    matrix[i][j].indices1 = i
                # Left
                elif (left_score >= 0):
                    matrix[i][j].prev = left
                    matrix[i][j].gaps = left.gaps + 1
                    matrix[i][j].value = left_score
                    matrix[i][j].indices2 = j                               
                # Fall through not necessary, it's handled during init.
                #else:
                #    matrix[i][j].prev = None
                #    matrix[i][j].value = 0
                
                # If the value is greater than 0, we'll consider it for now.
                #  We need to preserve the order, so add to the front of the list.
                if matrix[i][j].gaps <= max_gaps:
                    optimals.append(matrix[i][j])
        
        # Sort the optimal list, order is preserved so the top valued node deepest into
        #  the matrix is considered first.
        optimals.sort(key=lambda node: node.value)
        
        # We could probably do some post filtering here as well...
        #  for now, lets just return the top node
        #optimals = optimals[:1]
        
        seqs = []
        while optimals:
            optimal = optimals.pop()
            
            seq1 = []
            seq2 = []
            node = optimal
            
            while node != None:
                # If this node is in the list, remove it so we don't process this path again!
                try:
                    optimals.remove(node)
                except ValueError:
                    pass
                
                # The true indices need to be shifted back one due to DP.
                index1 = node.indices1
                if index1 != None:
                    index1 -= 1
                index2 = node.indices2
                if index2 != None:
                    index2 -= 1
                
                seq1.insert(0, index1)
                seq2.insert(0, index2)
                
                node = node.prev
            
            # Remove unnecessary gap due to first node at 0,0.
            if seq1[0] == None and seq2[0] == None:
                seq1.pop(0)
                seq2.pop(0)
            
            # This will the trim a lot of unecessary alignments that might have passed the
            #  max gap threshold.
            if len(seq1) > 0 and len(seq2) > 0:
                seqs.append((optimal, seq1, seq2))
        
        return matrix, seqs