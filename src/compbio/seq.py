"""
Sequence methods for global and local alignment/comparison.

@author  azampagl@azampagl.com (Aaron Zampaglione)
@copyright MIT

@see Setubal, Medianis. Introduction to Computational Biology.
@see http://en.wikipedia.org/wiki/Sequence_alignment
@see http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
@see http://www.clcbio.com/index.php?id=1046
"""
import numpy as np

def global_align(seq1, seq2):
    """Globally aligns two sequences and 
    returns the alignment score and two
    globally aligned sequences.
    
    Keyword arguments:
    seq1  -- the first sequence.
    seq2  -- the second sequence.
    """ 
    # Init score matrix.
    (rows, cols), score = __cmp_init(seq1, seq2, -2)
    
    # Build a matrix full of Prefix objects
    vprefix = np.vectorize(__Prefix)
    prefixes = np.empty((rows, cols), dtype=object)
    prefixes[:,:] = vprefix(score)
    
    # Determine score using DP.
    for i in range(1, rows):
        for j in range(1, cols):
            
            # Find the previous top, diag, and left scores
            top = prefixes[i - 1][j].score - 2
            diag = prefixes[i - 1][j - 1].score + __p(seq1[i - 1], seq2[j - 1])
            left = prefixes[i][j - 1].score - 2
            
            # Top
            if (top > diag and top > left):
                prefixes[i][j].prev = prefixes[i - 1][j]
                prefixes[i][j].score = top
                prefixes[i][j].seq1 = seq1[i - 1]
            # Diagonal
            elif (diag > top and diag > left):
                prefixes[i][j].prev = prefixes[i - 1][j - 1]
                prefixes[i][j].score = diag
                prefixes[i][j].seq1 = seq1[i - 1]
                prefixes[i][j].seq2 = seq2[j - 1]
            # Left
            else:
                prefixes[i][j].prev = prefixes[i][j - 1]
                prefixes[i][j].score = left
                prefixes[i][j].seq2 = seq2[j - 1]
    
    # Follow the pointers backwards to rebuild the globally aligned sequences.
    s1 = ""
    s2 = ""
    node = prefixes[i][j]
    while node.prev != None:
        s1 = node.seq1 + s1
        s2 = node.seq2 + s2
        node = node.prev
    
    # Return the score
    return s1, s2, prefixes[-1][-1].score

def local_align(seq1, seq2):
    """Returns the optimal (longest) local alignment (substring) of two sequences.
    
    Keyword arguments:
    seq1 -- the first sequence.
    seq2 -- the second sequence.
    """
    # Init score matrix.
    (rows, cols), score = __cmp_init(seq1, seq2, 0)
    
    # Build a matrix full of Prefix objects
    vprefix = np.vectorize(__Prefix)
    prefixes = np.empty((rows, cols), dtype=object)
    prefixes[:,:] = vprefix(score)
    
    # Initialize max/optimal
    optimal = prefixes[0][0]
    
    # Determine score using DP.
    for i in range(1, rows):
        for j in range(1, cols):
            
            # Find the previous top, diag, and left scores
            top = prefixes[i - 1][j].score - 2
            diag = prefixes[i - 1][j - 1].score + __p(seq1[i - 1], seq2[j - 1])
            left = prefixes[i][j - 1].score - 2
            
            # Top
            if (top > diag and top > left and top > 0):
                prefixes[i][j].prev = prefixes[i - 1][j]
                prefixes[i][j].score = top
                prefixes[i][j].seq1 = seq1[i - 1]
            # Diagonal
            elif (diag > top and diag > left and diag > 0):
                prefixes[i][j].prev = prefixes[i - 1][j - 1]
                prefixes[i][j].score = diag
                prefixes[i][j].seq1 = seq1[i - 1]
                prefixes[i][j].seq2 = seq2[j - 1]
            # Left
            elif (left > top and left > diag and left > 0):
                prefixes[i][j].prev = prefixes[i][j - 1]
                prefixes[i][j].score = left
                prefixes[i][j].seq2 = seq2[j - 1]
    
            if prefixes[i][j].score > optimal.score:
                optimal = prefixes[i][j]
            
    # Follow the pointers backwards to rebuild the optimal aligned sequence.
    s = ""
    node = optimal
    while node.prev != None:
        s = node.seq1 + s
        node = node.prev
    
    # Return the optimal sequence
    return s

class __Prefix():
    """
    Prefix object for a similarity matrix.
    """
    
    # Backpointer to the previous prefix.
    prev = None
    
    # Score at the current prefix
    score = 0
    
    # Sequence 1 modification at this prefix.
    seq1 = " "
    
    # Sequence 2 modification at this prefix.
    seq2 = " "
        
    def __init__(self, score):
        """
        Initializes the score during vectorization.
        
        Key arguments:
        score -- the score for this prefix.
        """
        self.score = score
            
def __cmp_init(seq1, seq2, v):
    """Initializes a similarity matrix.
    
    Keyword arguments:
    seq1 -- the first sequence.
    seq2 -- the second sequence.
    v    -- multiple for the first row and column.
    """
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    
    # Build our score matrix.
    score = np.zeros([rows, cols], dtype=np.int)
    
    # Set the first column of score values.
    for i in range(1, rows):
        score[i][0] = i * v
    
    # Set the first row of our score values.
    for i in range(1, cols):
        score[0][i] = i * v
    
    return (rows, cols), score

def __p(a, b):
    """Determines if two characters in a sequence are equivalent.
    
    Keyword arguments:
    a -- the character in the first sequence.
    b -- the character in the second sequence.
    """
    if (a == b):
        return 1
    
    return -1