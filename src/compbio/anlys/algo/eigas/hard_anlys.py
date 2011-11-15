"""
Analysis of hard alignment proteins
@author Aaron Zampaglione <azapagl@azampagl.com>
@package EIGAs
@copyright 2011 (c) Aaron Zampaglione
@license MIT
"""
from compbio.algo.eigas.core import EIGAs
from compbio.algo.eigas.protein.core import Protein
from compbio.algo.eigas.protein.parser.core import ProteinParser
from compbio.common.data import HARD

PROTEIN_PAIRS = [(HARD['1FXIa'], HARD['1UBQ'], 74),
                 (HARD['1TEN'], HARD['3HHRb'], 88),
                 (HARD['3HLAb'], HARD['2RHE'], 95),
                 (HARD['2AZAa'], HARD['1PAZ'], 109),
                 (HARD['1CEWi'], HARD['1MOLa'], 88),
                 (HARD['1CID'], HARD['2RHE'], 133),
                 (HARD['1CRL'], HARD['1EDE'], 304),
                 (HARD['2SIM'], HARD['1NSBa'], 339),
                 (HARD['1BGEb'], HARD['2GMFa'], 121),
                 (HARD['1TIE'], HARD['4FGF'], 120),
                 ]

# Title
print('Protein1\t|\tProtein2\t|\tScore\t|\tReport\'s Score\n')

for pair in PROTEIN_PAIRS:
    protein1 = Protein(ProteinParser.factory('pdb', (pair[0])))
    protein2 = Protein(ProteinParser.factory('pdb', (pair[1])))
    aligned = EIGAs.aligned(protein1=protein1, protein2=protein2)
        
    print(protein1.name + ' (' + str(len(protein1.fingerprint)) + ')' + \
          '\t|\t' + \
          protein2.name + ' (' + str(len(protein2.fingerprint)) + ')' + \
          '\t|\t' + str(aligned) + \
          '\t|\t' + str(pair[2])
          )