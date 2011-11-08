"""
@author Aaron Zampaglione <azapagl@azampagl.com>
@package EIGAs
@copyright 2011 Aaron Zampaglione
@license MIT
"""
import unittest

from itertools import product, combinations

from compbio.algo.eigas.protein.core import Protein
from compbio.algo.eigas.protein.parser.core import ProteinParser

class TestCompbioAlgoEIGAsHard(unittest.TestCase):
    
    # Root directory.
    ROOT = '../../../../../'
    
    # Data directory.
    DATA_DIR = ROOT + 'data/hard/pdb'
    
    # .pdb style obtained from pdb.org
    PROTEINS = [('1FXI', DATA_DIR + '/1FXI.pdb'),
                ('1UBQ', DATA_DIR + '/1UBQ.pdb'),
                ('1TEN', DATA_DIR + '/1TEN.pdb'),
                ('3HHR', DATA_DIR + '/3HHR.pdb'),
                ('3HLA', DATA_DIR + '/3HLA.pdb'),
                ('2RHE', DATA_DIR + '/2RHE.pdb'),
                ('2AZA', DATA_DIR + '/2AZA.pdb'),
                ('1PAZ', DATA_DIR + '/1PAZ.pdb'),
                ('1CEW', DATA_DIR + '/1CEW.pdb'),
                ('1MOL', DATA_DIR + '/1MOL.pdb'),
                ('1CID', DATA_DIR + '/1CID.pdb'),
                ('1CRL', DATA_DIR + '/1CRL.pdb'),
                ('1EDE', DATA_DIR + '/1EDE.pdb'),
                ('2SIM', DATA_DIR + '/2SIM.pdb'),
                ('1NSB', DATA_DIR + '/1NSB.pdb'),
                ('1BGE', DATA_DIR + '/1BGE.pdb'),
                ('2GMF', DATA_DIR + '/2GMF.pdb'),
                ('1TIE', DATA_DIR + '/1TIE.pdb'),
                ('4FGF', DATA_DIR + '/4FGF.pdb'),
                ]
    
    def setUp(self):
        """
        Build the initial protein pairs.
        """
        self.proteins = []
        
        # Loop through all the proteins and get protein objects.
        for protein in self.PROTEINS:
            self.proteins.append(Protein(ProteinParser.factory('pdb', protein)))
    
    def tearDown(self):
        """
        """
        pass
    
    def testFingerprintMax(self):
        """
        """
        for protein1, protein2 in combinations(self.proteins, 2):
            max_diff = max([abs(fingerprinti - fingerprintj) for fingerprinti, fingerprintj in product(protein1.fingerprint, protein2.fingerprint)])
            print(protein1.name + ' (' + str(len(protein1.fingerprint)) + ')\tvs.\t' + protein2.name + ' (' + str(len(protein2.fingerprint)) + ')\t' + str(max_diff))
            #print(protein.name + ': ' + str(min(protein.finger + print)) + ' (Min)')
            #print('\t' + protein.name + '(' + str(len(protein.fingerprint)) + '):\t' + str(max(protein.fingerprint)))

if __name__ == "__main__":
    unittest.main()