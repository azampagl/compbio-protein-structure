"""
Unit tests for the EIGA class.

@author Aaron Zampaglione <azapagl@azampagl.com>
@package EIGAs
@copyright 2011 (c) Aaron Zampaglione
@license MIT
"""
import unittest

from compbio.algo.eigas.core import EIGAs
from compbio.algo.eigas.protein.core import Protein
from compbio.algo.eigas.protein.parser.core import ProteinParser
from compbio.common.data import HARD

class TestCompbioAlgoEIGAs(unittest.TestCase):
    
    def testSelfAlign(self):
        """
        Tests self aligning a protein to itself.
        """
        protein = Protein(ProteinParser.factory('pdb', ('1BGE', HARD['1BGE'][0], HARD['1BGE'][1])))
        
        score, seq1, seq2 = EIGAs.global_align(protein, protein)
        self.assertEqual(score, 0.0)
        self.assertEqual(seq1, seq2)

if __name__ == "__main__":
    unittest.main()