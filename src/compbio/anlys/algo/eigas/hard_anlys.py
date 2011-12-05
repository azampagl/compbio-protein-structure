"""
Analysis of hard alignment proteins
@author Aaron Zampaglione <azapagl@azampagl.com>
@package EIGAs
@copyright 2011 Aaron Zampaglione
@license MIT
"""
from compbio.algo.eigas.core import EIGAs
from compbio.algo.eigas.protein.core import Protein
from compbio.algo.eigas.protein.parser.core import ProteinParser
from compbio.common.data import HARD

import getopt
import sys

# Get command line args.
opts, args = getopt.getopt(sys.argv[1:], ':o')
if not len(args):
    raise Exception('Missing output file.')

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
html = """
<table>
    <tr>
        <th>Protein 1</th>
        <th>Protein 2</th>
        <th>Alignment</th>
        <th>Report's Alignment</th>
    </tr>
"""

for pair in PROTEIN_PAIRS:
    protein1 = Protein(ProteinParser.factory('pdb', (pair[0])))
    protein2 = Protein(ProteinParser.factory('pdb', (pair[1])))
    _, s1, s2 =  EIGAs.global_align(protein1, protein2)
    aligned = EIGAs.aligned(protein1=protein1, protein2=protein2)
    
    s = """
    <tr>
        <td>{0} ({1})</td>
        <td>{2} ({3})</td>
        <td>{4}</td>
        <td>{5}</td>
    </tr>
    """
    html += s.format(protein1.name, len(protein1.fingerprint),
               protein2.name, len(protein2.fingerprint),
               aligned,
               pair[2])

html += """
</table>
"""

open(args[0], 'w').write(html)
print('Complete.')