"""
Sequence methods for global and local alignment/comparison.

@author  azampagl@azampagl.com (Aaron Zampaglione)
@copyright MIT

@see Setubal, Medianis. Introduction to Computational Biology.
@see http://en.wikipedia.org/wiki/Sequence_alignment
@see http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
@see http://www.clcbio.com/index.php?id=1046
"""
from numpy import int, zeros

def global_align(seq1, seq2):
    """Globally aligns two sequences.
    
    Keyword arguments:
    seq1 -- the first sequence.
    seq2 -- the second sequence.
    """
    return __global_align(seq1, seq2, len(seq1), len(seq2), __global_cmp(seq1, seq2))

def global_cmp(seq1, seq2):
    """Globally similarity between two sequences.
    
    Keyword arguments:
    seq1 -- the first sequence.
    seq2 -- the second sequence.
    """
    return __global_cmp(seq1, seq2)[-1][-1]

def local_align(seq1, seq2, threshold = 2):
    """Returns the optimal (longest) local alignment (substring) of two sequences.
    
    This method can be easily modified to return other local
    alignment strings found between the two sequences.  Currently,
    it only returns the optimal (longest similar substring).
    
    Keyword arguments:
    seq1      -- the first sequence.
    seq2      -- the second sequence.
    threshold -- the minimum substring length (default 2).
    """
    indices, score = __local_cmp(seq1, seq2, threshold)
    return __local_align(seq1, seq2, indices[0][0], indices[0][1], score)


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
    score = zeros([rows, cols], dtype=int)
    
    # Set the first column of score values.
    for i in range(1, rows):
        score[i][0] = i * v
    
    # Set the first row of our score values.
    for i in range(1, cols):
        score[0][i] = i * v
    
    return (rows, cols), score

def __global_align(seq1, seq2, i, j, score):
    """
    @see global_align
    
    Keyword arguments:
    seq1  -- the first sequence.
    seq2  -- the second sequence.
    i     -- the index for the first sequence (row).
    j     -- the index for the second sequence (col).
    score -- the score matrix.
    """
    s = ""
    t = ""
    
    # Base case.
    if (i == 0 and j == 0):
        s = ""
        t = ""
        
    # Left operation.
    elif (i > 0 and score[i][j] == (score[i - 1][j] - 2)):
        s, t = __global_align(seq1, seq2, i - 1, j, score)
        
        s += seq1[i - 1]
        t += " "
        
    # Diagonal operation.
    elif (i > 0 and j > 0 and score[i][j] == (score[i - 1][j - 1] + __p(seq1[i - 1], seq2[j - 1]))):
        s, t = __global_align(seq1, seq2, i - 1, j - 1, score)
        
        s += seq1[i - 1]
        t += seq2[j - 1]
        
    # Top operation.
    else:
        s, t = __global_align(seq1, seq2, i, j - 1, score)
        
        s += " "
        t += seq2[j - 1]
        
    return s, t
        
def __global_cmp(seq1, seq2):
    """Builds the global similarity matrix for two sequences.
    
    Keyword arguments:
    seq1  -- the first sequence.
    seq2  -- the second sequence.
    """
    # Init score matrix.
    (rows, cols), score = __cmp_init(seq1, seq2, -2)
    
    # Determine score using DP.
    for i in range(1, rows):
        for j in range(1, cols):
            # Calculate score
            score[i][j] = max(score[i - 1][j] - 2,
                              score[i - 1][j - 1] + __p(seq1[i - 1], seq2[j - 1]),
                              score[i][j - 1] - 2)
    
    # Return the score
    return score

def __local_align(seq1, seq2, i, j, score):
    """
    @see local_align
    
    Keyword arguments:
    seq1  -- the first sequence.
    seq2  -- the second sequence.
    i     -- the index for the first sequence (row).
    j     -- the index for the second sequence (col).
    score -- the score matrix.
    """
    
    # Base case.
    if (i == 0 and j == 0):
        return ""
    
    # Found a zero value
    elif (score[i - 1][j - 1] == 0):
        return seq1[i - 1]
    
    # Left operation.
    elif (i > 0 and score[i][j] == (score[i - 1][j] - 2)):
        return __local_align(seq1, seq2, i - 1, j, score) + seq1[i - 1]
    
    # Diagonal operation.
    elif (i > 0 and j > 0 and score[i][j] == (score[i - 1][j - 1] + __p(seq1[i - 1], seq2[j - 1]))):
        return __local_align(seq1, seq2, i - 1, j - 1, score) + seq1[i - 1]
    
    # Top operation.
    return __local_align(seq1, seq2, i, j - 1, score) + " "

def __local_cmp(seq1, seq2, threshold):
    """Builds the local similarity matrix for two sequences.
    
    Keyword arguments:
    seq1      -- the first sequence.
    seq2      -- the second sequence.
    threshold -- the minimum substring length (default 2).
    """
    # Init score matrix.
    (rows, cols), score = __cmp_init(seq1, seq2, 0)
    
    # Dictionary of the max scores
    max_scores = {}
    
    # Determine score using DP.
    for i in range(1, rows):
        for j in range(1, cols):
            # Calculate score
            score[i][j] = max(score[i][j - 1] - 2,
                              score[i - 1][j - 1] + __p(seq1[i - 1], seq2[j - 1]),
                              score[i - 1][j] - 2,
                              0)
            
            # Add to max scores if the value was above the given threshold.
            if (score[i][j] >= threshold):
                if (score[i][j] in max_scores):
                    max_scores[score[i][j]].append((i, j))
                else:
                    max_scores[score[i][j]] = [(i, j)]
    
    # Sort the indices by the greatest (optimal) values in the matrix.
    indices = []
    for key in sorted(max_scores.keys(), reverse=True):
        indices.extend(max_scores[key])
        
    # Return the sorted indices and the score matrix
    return indices, score

def __p(a, b):
    """Determines if two characters in a sequence are equivalent.
    
    Keyword arguments:
    a -- the character in the first sequence.
    b -- the character in the second sequence.
    """
    if (a == b):
        return 1
    
    return -1