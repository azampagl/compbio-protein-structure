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
from compbio.common.data import SKOLNICK

import getopt
import itertools
import sys

# Get command line args.
opts, args = getopt.getopt(sys.argv[1:], ':o')
if not len(args):
    raise Exception('Missing output file.')

html = """
"""

for family in SKOLNICK:
    title = """
    <h1>{0}</h1>
    """.format(family)
    
    # Table Headers
    html += title + """
    <table>
        <tr>
            <th>Protein 1</th>
            <th>Protein 2</th>
            <th>Aligned</th>
        </tr>
    """

    for protein1, protein2 in itertools.combinations(SKOLNICK[family].keys(), 2):
        protein1 = Protein(ProteinParser.factory('pdb', (SKOLNICK[family][protein1])))
        protein2 = Protein(ProteinParser.factory('pdb', (SKOLNICK[family][protein2])))
        aligned = EIGAs.aligned(protein1=protein1, protein2=protein2)
    
        html += """
        <tr>
            <td>{0} ({1})</td>
            <td>{2} ({3})</td>
            <td>{4}</td>
        </tr>
        """.format(protein1.name, len(protein1.fingerprint),
                   protein2.name, len(protein2.fingerprint),
                   aligned)

    html += """
    </table>
    """

open(args[0], 'w').write(html)
print('Complete.')