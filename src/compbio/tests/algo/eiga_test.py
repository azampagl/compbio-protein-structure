'''
Unit tests for the EIGA class.

@author Aaron Zampaglione <azapagl@azampagl.com>
'''
import unittest

from compbio.algo.eiga import Eiga

class TestCompbioAlgoEiga(unittest.TestCase):
    
    # Eiga class.
    eiga = None
    
    # Structures to test.
    structures = [("d1a04a2", "../../../../data/ent/d1a04a2.ent"),
                  ("d1axib2", "../../../../data/ent/d1axib2.ent")]
    
    def setUp(self):
        """
        """
        pass
    
    def tearDown(self):
        """
        """
        pass
    
    def testCmatrix(self):
        """
        """
        eiga = Eiga([self.structures[0]])
        eiga.cmatrices
    
    def testCoords(self):
        """
        """
        Eiga.coords(self.structures[0][0], self.structures[0][1])
    
    def testDmatrix(self):
        """
        """
        coords = Eiga.coords(self.structures[0][0], self.structures[0][1])
        Eiga.dmatrix(coords)


if __name__ == "__main__":
    unittest.main()