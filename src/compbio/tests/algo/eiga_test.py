'''
Created on Sep 15, 2011

@author: sysadmin
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
        self.eiga = Eiga(self.structures)
    
    def tearDown(self):
        pass
    
    def testName(self):
        print(self.eiga.d)
    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()