"""
Unit tests for the EIGA class.

@see http://astral.berkeley.edu/pdbstyle-1.75.html

@author Aaron Zampaglione <azapagl@azampagl.com>
@package EIGAs
@copyright 2011 (c) Aaron Zampaglione
@license MIT
"""
import unittest

from compbio.algo.eigas.core import EIGAs
from compbio.algo.eigas.protein.core import Protein
from compbio.algo.eigas.protein.parser.core import ProteinParser
from compbio.tests.common.data import HARD

class TestCompbioAlgoEIGAs(unittest.TestCase):
    
    # Protein pairs to test.
    #   Protein 1, Protein 2, Alignment Score
    '''PROTEIN_PAIRS = [(('1FXI', DATA_DIR + '/1FXI.pdb'), ('1UBQ', DATA_DIR + '/1UBQ.pdb'), 74), # 1.
                     (('1TEN', DATA_DIR + '/1TEN.pdb'), ('3HHR', DATA_DIR + '/3HHR.pdb'), 88), # 2.
                     (('3HLA', DATA_DIR + '/3HLA.pdb'), ('2RHE', DATA_DIR + '/2RHE.pdb'), 95), # 3.
                     (('2AZA', DATA_DIR + '/2AZA.pdb'), ('1PAZ', DATA_DIR + '/1PAZ.pdb'), 109), # 4.
                     (('1CEW', DATA_DIR + '/1CEW.pdb'), ('1MOL', DATA_DIR + '/1MOL.pdb'), 88), # 5.
                     (('1CID', DATA_DIR + '/1CID.pdb'), ('2RHE', DATA_DIR + '/2RHE.pdb'), 113), # 6.
                     (('1CRL', DATA_DIR + '/1CRL.pdb'), ('1EDE', DATA_DIR + '/1EDE.pdb'), 304), # 7.
                     (('2SIM', DATA_DIR + '/2SIM.pdb'), ('1NSB', DATA_DIR + '/1NSB.pdb'), 339), # 8.
                     (('1BGE', DATA_DIR + '/1BGE.pdb'), ('2GMF', DATA_DIR + '/2GMF.pdb'), 121), # 9.
                     (('1TIE', DATA_DIR + '/1TIE.pdb'), ('4FGF', DATA_DIR + '/4FGF.pdb'), 120), # 10.
                     ]
    '''
    # Store Protein objects.
    protein_pairs = []
    
    def setUp(self):
        """
        Build the initial protein pairs.
        """
        """
        self.protein_pairs = []
        
        # Loop through all the proteins and get protein objects.
        for protein_pair in self.PROTEIN_PAIRS:
            protein1 = Protein(ProteinParser.factory('pdb', protein_pair[0]))
            protein2 = Protein(ProteinParser.factory('pdb', protein_pair[1]))
            self.protein_pairs.append((protein1,
                                       protein2,
                                       protein_pair[2]
                                       ))
        """
    
    def tearDown(self):
        """
        """
        pass
    
    def testAlignments(self):
        """
        """
        '''
        for protein_pair in self.protein_pairs:
            protein1 = protein_pair[0]
            protein2 = protein_pair[1]
            _, seq1, seq2 = EIGAs.global_align(protein1, protein2)
            aligned = EIGAs.aligned(seq1=seq1, seq2=seq2)
            
            msg = protein1.name + ' (' + str(len(protein1.fingerprint)) + ')' +  '\tvs.\t' + \
                  protein2.name + ' (' + str(len(protein2.fingerprint)) + ')' + '\t=>\t' + \
                  str(aligned) + ' <> ' + str(protein_pair[2])# + '\n' + \
                  #', '.join([str(value) for value in seq1]) + '\n' + \
                  #', '.join([str(value) for value in seq2]) + '\n'
            
            print(msg)
            #self.assertTrue(aligned == protein_pair[2], msg)
        '''
        print(SKOLNICK)

if __name__ == "__main__":
    unittest.main()