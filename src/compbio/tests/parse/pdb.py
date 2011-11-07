"""
Unit tests for parsing pdb files.

This is to verify the Biopython is installed correctly
and we know how to use the basic parsing method.

@author Aaron Zampaglione <azapagl@azampagl.com>
@package compbio
@copyright 2011 Aaron Zampaglione
@license MIT
"""
import unittest

from decimal import Decimal

from Bio.PDB.PDBParser import PDBParser

class TestPDBParse(unittest.TestCase):
    
    # Root directory.
    ROOT = '../../../../'
    
    # Data directory.
    DATA_DIR = ROOT + 'data/hard/pdb'
    
    # Folds to test.
    
    # .pdb style obtained from pdb.org
    # Format: Protein Name, PDB File Location, Number of Alpha Carbons.
    PROTEINS = [('1FXI', DATA_DIR + '/1FXI.pdb', 96),
                ]
    
    def setUp(self):
        """
        """
        self.parser = PDBParser()
    
    def tearDown(self):
        """
        """
        pass
    
    def testParse(self):
        """
        Tests if Biopython's PDBParser.
        """
        for protein in self.PROTEINS:
            structure = self.parser.get_structure(protein[0], protein[1])
            
            # Unsure how to process multiple structures.
            if len(structure.get_list()) < 1:
                raise Exception("Error reading structure for '" + protein[0] + '"')
            
            model = structure.get_list()[0]
            
            # Unsure how to process multiple models.
            if len(model.get_list()) < 1:
                raise Exception("Error reading model for '" + protein[0] + '"')
            
            # The chain contains a list of residues.
            chain = model.get_list()[0]
            
            coords = []
            
            # Look at all the alpha carbon atoms.
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == 'CA':
                        coords.append(map(Decimal, (map(str, atom.get_coord()))))
            
            self.assertTrue(len(coords) == protein[2], "Invalid number of alpha carbons.")

if __name__ == "__main__":
    unittest.main()