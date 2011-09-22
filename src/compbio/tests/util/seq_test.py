"""
Unit tests for the sequence modules of compbio.

@author  azampagl@azampagl.com (Aaron Zampaglione)
@copyright MIT
"""
import unittest

from compbio.util.seq import global_align, local_align

class TestCompbioUtilSeq(unittest.TestCase):
    """
    """

    def setUp(self):
        """Build global and local alignment fixtures."""
        
        # Global fixtures.
        self.global_seqs = [
            ("GACGGATTAG", "GATCGGAATAG", "GA CGGATTAG", "GATCGGAATAG", 6),
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
            self.assertTrue(global_align(seq[0], seq[1]) == (seq[2], seq[3], seq[4]))

    def test_local_align(self):
        """Test local alignment."""
        for seq in self.local_seqs:
            self.assertTrue(local_align(seq[0], seq[1]) == seq[2])
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testGlobalCmp01']
    unittest.main()