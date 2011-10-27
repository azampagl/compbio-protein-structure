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
from compbio.algo.eigas.exception import EIGAsException

class TestCompbioAlgoEIGAs(unittest.TestCase):
    
    # Root directory.
    ROOT = '../../../../../'
    
    # Data directory.
    DATA_DIR = ROOT + '/data/hard'
    
    # Folds to test.
    PROTEIN_PAIRS = [(('d1fxia_', DATA_DIR + '/d1fxia_.ent'), ('d1ubqa_', DATA_DIR + '/d1ubqa_.ent')), # 1.
                     (('d1crla_', DATA_DIR + '/d1crla_.ent'), ('d1edea_', DATA_DIR + '/d1edea_.ent'))] # 7.
    
    # Store Protein objects.
    protein_pairs = []
    
    def setUp(self):
        """
        Build the initial protein pairs.
        """
        self.protein_pairs = []
        
        # Loop through all the proteins and get protein objects.
        for protein_pair in self.PROTEIN_PAIRS:
            self.protein_pairs.append((EIGAs.Protein(protein_pair[0][0], protein_pair[0][1]),
                                       EIGAs.Protein(protein_pair[1][0], protein_pair[1][1]),
                                       ))
        
    def tearDown(self):
        """
        """
        pass
    
    def testAlignments(self):
        """
        """
        for protein_pair in self.protein_pairs:
            print(protein_pair[0].name + ' vs. ' + protein_pair[1].name + ': ' + \
                  str(EIGAs.aligned(protein1=protein_pair[0], protein2=protein_pair[1])))

if __name__ == "__main__":
    unittest.main()