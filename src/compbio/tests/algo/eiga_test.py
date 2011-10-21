'''
Unit tests for the EIGA class.

@see http://astral.berkeley.edu/pdbstyle-1.75.html

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
                
                'tim-beta': [("d1amka_", "../../../../data/ent/d1amka_.ent")],
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
        family1 = self.skolnick['flavodxin-like']
        family2 = self.skolnick['cupredoxin-like']
        
        protein1_1 = Eiga.Protein(family1[0][0], family1[0][1])
        protein1_2 = Eiga.Protein(family1[1][0], family1[1][1])
        
        protein2_1 = Eiga.Protein(family2[0][0], family2[0][1])
        
        
        print(str(family1[0][0]) + " " + str(family1[1][0]))
        print(Eiga.align(protein1_1, protein1_2))
        print("")
        print(str(family1[0][0]) + " " + str(family2[0][0]))
        print(Eiga.align(protein1_1, protein2_1))
        
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