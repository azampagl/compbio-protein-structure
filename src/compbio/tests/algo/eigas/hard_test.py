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

class TestCompbioAlgoEIGAs(unittest.TestCase):
    
    # Root directory.
    ROOT = '../../../../../'
    
    # Data directory.
    DATA_DIR = ROOT + '/data/hard/ent'
    
    # Folds to test.
    
    # .pdb style obtained from pdb.org
    PROTEIN_PAIRS = [(('d1fxia_', DATA_DIR + '/d1fxia_.ent'), ('d1ubqa_', DATA_DIR + '/d1ubqa_.ent')), # 1.
                     ]
    
    # .ent style obtained from berkley
    '''
    PROTEIN_PAIRS = [(('d1fxia_', DATA_DIR + '/d1fxia_.ent'), ('d1ubqa_', DATA_DIR + '/d1ubqa_.ent')), # 1.
                     (('d1tena_', DATA_DIR + '/d1tena_.ent'), ('d3hhrb1', DATA_DIR + '/d3hhrb1.ent')), # 2.
                     (('d3hlab_', DATA_DIR + '/d3hlab_.ent'), ('d2rhea_', DATA_DIR + '/d2rhea_.ent')), # 3.
                     (('d2azaa_', DATA_DIR + '/d2azaa_.ent'), ('d1paza_', DATA_DIR + '/d1paza_.ent')), # 4.
                     (('d1cewi_', DATA_DIR + '/d1cewi_.ent'), ('d1mola_', DATA_DIR + '/d1mola_.ent')), # 5.
                     (('d1cida1', DATA_DIR + '/d1cida1.ent'), ('d2rhea_', DATA_DIR + '/d2rhea_.ent')), # 6.
                     (('d1crla_', DATA_DIR + '/d1crla_.ent'), ('d1edea_', DATA_DIR + '/d1edea_.ent')), # 7.
                     (('d2sima_', DATA_DIR + '/d2sima_.ent'), ('d1nsba_', DATA_DIR + '/d1nsba_.ent')), # 8.
                     (('d1bgeb_', DATA_DIR + '/d1bgeb_.ent'), ('d2gmfa_', DATA_DIR + '/d2gmfa_.ent')), # 9..
                     (('d1tiea_', DATA_DIR + '/d1tiea_.ent'), ('d4fgfa_', DATA_DIR + '/d4fgfa_.ent')), # 10.
                     ]
    '''
    
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
            protein1 = protein_pair[0]
            protein2 = protein_pair[1]
            score, seq1, seq2 = EIGAs.global_align(protein1, protein2)
            print(protein1.name + ' (' + str(len(protein1.fingerprint)) + ')' +  '\tvs.\t' + \
                  protein2.name + ' (' + str(len(protein2.fingerprint)) + ')' + ':\t' + \
                  str(EIGAs.aligned(seq1=seq1, seq2=seq2)))
            #print(seq1)
            #print(seq2)
            #print('')

if __name__ == "__main__":
    unittest.main()