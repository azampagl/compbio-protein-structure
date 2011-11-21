"""
Analysis of skolnick alignment proteins.

@author Aaron Zampaglione <azapagl@azampagl.com>
@package EIGAs
@copyright 2011 Aaron Zampaglione
@license MIT
"""
from compbio.algo.eigas.core import EIGAs
from compbio.algo.eigas.protein.core import Protein
from compbio.algo.eigas.protein.parser.core import ProteinParser
from compbio.common.data import SKOLNICK

from numpy import average, std

from getopt import getopt
from itertools import product
from sys import argv

def combinations_with_replacement(iterable, r):
    # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
    pool = tuple(iterable)
    n = len(pool)
    for indices in product(range(n), repeat=r):
        if sorted(indices) == list(indices):
            yield tuple(pool[i] for i in indices)

# Get command line args.
opts, args = getopt(argv[1:], ':o')
if not len(args):
    raise Exception('Missing output file.')

html = """
"""

for family1, family2 in combinations_with_replacement(SKOLNICK.keys(), 2):
    title = """
    <h1>{0} vs. {1}</h1>
    """.format(family1, family2,)
    
    # Table Headers
    html += title + """
    <table>
        <tr>
            <th>Protein 1</th>
            <th>Protein 2</th>
            <th>Aligned</th>
            <th>Alignment Length</th>
            <th>Percentage</th>
        </tr>
    """
    alignments = []
    percentages = []
    for protein1, protein2 in product(SKOLNICK[family1].keys(), SKOLNICK[family2].keys()):
        protein1 = Protein(ProteinParser.factory('pdb', (SKOLNICK[family1][protein1])))
        protein2 = Protein(ProteinParser.factory('pdb', (SKOLNICK[family2][protein2])))
        matrix, seq1, seq2 = EIGAs.global_align(protein1, protein2)
        aligned = EIGAs.aligned(seq1=seq1, seq2=seq2)
    
        html += """
        <tr>
            <td>{0} ({1})</td>
            <td>{2} ({3})</td>
            <td>{4}</td>
            <td>{5}</td>
            <td>{6}</td>
        </tr>
        """.format(protein1.name, len(protein1.fingerprint),
                   protein2.name, len(protein2.fingerprint),
                   aligned,
                   len(seq1),
                   aligned / float(len(seq1)))
        
        alignments.append(aligned)
        percentages.append(aligned / float(len(seq1)))
        
    html += """
        <tr>
            <td colspan="3">&nbsp;</td>
        </tr>
        <tr>
            <td colspan="2">Aligned Average:</td>
            <td>{0}</td>
        </tr>
        <tr>
            <td colspan="2">Aligned STDEV:</td>
            <td>{1}</td>
        </tr>
        <tr>
            <td colspan="2">Percentage Average:</td>
            <td>{2}</td>
        </tr>
        <tr>
            <td colspan="2">Percentage STDEV:</td>
            <td>{3}</td>
        </tr>
    </table>
    """.format(average(alignments),
               std(alignments),
               average(percentages),
               std(percentages))

open(args[0], 'w').write(html)
print('Complete.')