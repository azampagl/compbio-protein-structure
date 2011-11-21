"""
Basic protein for the EIGAs algorithm.

Determining a fingerprint for a protein is based on the report:
A Spectral Approach to Protein Structure
Alignment
Yosi Shibberu, Allen Holder
Mathematics Department Rose-Hulman Institute of Technology
Terre Haute, IN 47803

@author Aaron Zampaglione <azampagl@azampagl.com>
@package EIGAs
@copyright 2011 Aaron Zampaglione
@license MIT
"""
from numpy import diag, dot, float, zeros
from numpy.linalg import svd
from scipy.spatial.distance import cdist

class Protein(object):
    
    # Contact matrix cutoff.
    cutoff = 8
    
    @classmethod
    def cmatrix(cls, dmatrix):
        """
        Creates a contact matrix.
        
        Key arguments:
        dmatrix -- distance matrix
        """
        l = len(dmatrix)
        rows = range(l)
        cols = range(l)
        
        k = 1 / Protein.cutoff
        
        cmatrix = zeros((l, l), dtype=float)
        
        for i in rows:
            for j in cols:
                value = dmatrix[i][j]
                if value >= 0.0 and value <= Protein.cutoff:
                    cmatrix[i][j] = 1 - k * value
        
        return cmatrix
    
    def __init__(self, parser):
        """
        Initializes a protein.
        
        Key arguments:
        name   -- the name of this protein.
        parser -- the protein parser.
        """
        # Initialize fingerprint.
        self.fingerprint = []
        
        # Set the protein name.
        self.name = parser.name()
        
        # Find the atomic coordinates from the parser.
        self.coords = parser.coords()
        
        # Calculate distance matrix
        dmatrix = cdist(self.coords, self.coords)
        
        # Create the contact matrix.
        cmatrix = Protein.cmatrix(dmatrix)
        
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
        r = dot(diag(eigvalues) ** 0.5, eigvectorsT)
        
        # Test proof 2
        #from numpy import inner
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
        for j in range(len(self.coords)):
            max_index = None
            max_value = None
                
            for i in range(len(r[:][0])):
                if (r[i][j] > max_value):
                    max_index = i
                    max_value = r[i][j]
                
            self.fingerprint.append(eigvalues[max_index])