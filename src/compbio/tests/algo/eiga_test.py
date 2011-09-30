'''
Unit tests for the EIGA class.

@author Aaron Zampaglione <azapagl@azampagl.com>
'''
import unittest

from compbio.algo.eiga import Eiga

class TestCompbioAlgoEiga(unittest.TestCase):
    
    # Proteins to test
    skolnick = {
                'flavodxin-like': [("d1b00a_", "../../../../data/ent/d1b00a_.ent"),
                                   ("d1dbwa_", "../../../../data/ent/d1dbwa_.ent")],
                
                'cupredoxin-like': [("d1bawa_", "../../../../data/ent/d1bawa_.ent")],
                }
    
    def setUp(self):
        """
        """
        pass
    
    def tearDown(self):
        """
        """
        pass
    
    def testAlign(self):
        """
        """
        #for fold, proteins in self.skolnick.items():
        #    for i in range(proteins):
        proteins1 = self.skolnick['flavodxin-like']
        proteins2 = self.skolnick['cupredoxin-like']
        
        protein1 = Eiga.Protein(proteins1[0][0], proteins1[0][1])
        protein2 = Eiga.Protein(proteins2[0][0], proteins2[0][1])
        
        Eiga.align(protein1, protein2)
        
    def testCmatrix(self):
        """
        """
        pass
        #eiga = Eiga.Protein(self.structures[0][0], self.structures[0][1], {'k': 8.0})
        #print(eiga.cmatrices)
        #for structure, eigvector in eiga.eigvectors.items():
        #    print(eigvector.diagonal())
    
    def testCoords(self):
        """
        """
        #Eiga.coords(self.structures[0][0], self.structures[0][1])
    
    def testDmatrix(self):
        """
        """
        #coords = Eiga.coords(self.structures[0][0], self.structures[0][1])
        #Eiga.dmatrix(coords)

if __name__ == "__main__":
    unittest.main()