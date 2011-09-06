"""
Unit tests for the sequence modules of compbio.

@author  azampagl@azampagl.com (Aaron Zampaglione)
@copyright MIT
"""
from compbio.seq import global_align, global_cmp, local_align
import unittest


class TestCompbioSeq(unittest.TestCase):
    """
    """

    def setUp(self):
        """Build global and local alignment fixtures."""
        
        # Global fixtures.
        self.global_seqs = [
            ("GACGGATTAG", "GATCGGAATAG", 6, "GA CGGATTAG", "GATCGGAATAG"),
        ]
        
        # Local alignment fixtures.
        self.local_seqs = [
            ("CAGCCUCGCUUAG", "AAUGCCAUUGACGG", "GCC"),
            ("ACACACTA", "AGCACACA", "ACACA"),
        ]

    def tearDown(self):
        """Do Nothing."""
        pass

    def test_global_align(self):
        """Test global alignment."""
        for seq in self.global_seqs:
            self.assertTrue(global_align(seq[0], seq[1]) == (seq[3], seq[4]))
    
    def test_global_cmp(self):
        """Test global similarity."""
        for seq in self.global_seqs:
            self.assertTrue(global_cmp(seq[0], seq[1]) == seq[2])

    def test_local_align(self):
        """Test local alignment."""
        for seq in self.local_seqs:
            self.assertTrue(local_align(seq[0], seq[1]) == seq[2])
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testGlobalCmp01']
    unittest.main()