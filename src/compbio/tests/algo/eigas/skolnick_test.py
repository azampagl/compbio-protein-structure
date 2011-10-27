'''
Unit tests for the EIGA class.

@see http://astral.berkeley.edu/pdbstyle-1.75.html

@author Aaron Zampaglione <azapagl@azampagl.com>
@package EIGAs
@copyright 2011 (c) Aaron Zampaglione
@license MIT
'''
import os
import unittest

from compbio.algo.eiga.core import EIGA
from compbio.algo.eiga.exception import EIGAException

class TestCompbioAlgoEIGA(unittest.TestCase):
    
    # Folds to test.
    FOLDS = ['microbial_ribonuclease', 'cupredoxin-like']
    
    ROOT = '../../../../../'
    
    # Proteins to test
    skolnick = {}
    
    def setUp(self):
        """
        """
        self.skolnick = {}
        
        data_dir = self.ROOT + '/data/skolnick'
        
        #for path in os.listdir(data_dir):
        #    if os.path.isdir(data_dir + '/' + path):
        #        self.skolnick[path] = []
        
        for fold in self.FOLDS:
            self.skolnick[fold] = []
            for ent in os.listdir(data_dir + '/' + fold):
                try:
                    self.skolnick[fold].append(EIGA.Protein(ent[:-4], data_dir + '/' + fold + '/' + ent))
                except EIGAException as e:
                    print(e)
        
    def tearDown(self):
        """
        """
        pass
    
    def testAlign(self):
        """
        """
        for fold1 in self.FOLDS:
            for fold2 in self.FOLDS:
                for protein1 in self.skolnick[fold1]:
                    for protein2 in self.skolnick[fold2]:
                        if protein1 != protein2:
                            print(protein1.name + ' (' + str(fold1) + ') vs. ' + protein2.name + ' (' + str(fold2) + ')')
                            score, seq1, seq2, = EIGA.align(protein1, protein2)
                            aligned = 0
                            #print(seq1)
                            #print(seq2)
                            for i in range(len(seq1)):
                                if seq1[i] != None and seq2[i] != None:
                                    aligned += 1
                            print(str(aligned) + ' / ' + str(len(seq1)))
                            print('')
        '''family1 = self.skolnick['flavodxin-like']
        family2 = self.skolnick['cupredoxin-like']
        
        protein1_1 = Eiga.Protein(family1[0][0], family1[0][1])
        protein1_2 = Eiga.Protein(family1[1][0], family1[1][1])
        
        protein2_1 = Eiga.Protein(family2[0][0], family2[0][1])
        
        
        print(str(family1[0][0]) + " " + str(family1[1][0]))
        print(Eiga.align(protein1_1, protein1_2))
        print("")
        print(str(family1[0][0]) + " " + str(family2[0][0]))
        print(Eiga.align(protein1_1, protein2_1))'''
        
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