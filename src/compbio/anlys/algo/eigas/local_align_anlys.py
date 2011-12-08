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
from compbio.common.data import MOTIF

import getopt
import sys
import pprint

# Get command line args.
#opts, args = getopt.getopt(sys.argv[1:], ':o')
#if not len(args):
#    raise Exception('Missing output file.')

protein1 = Protein(ProteinParser.factory('pdb', MOTIF['MICROBIAL-RIBONUCLEASE']['1RDS']))
protein2 = Protein(ProteinParser.factory('pdb', MOTIF['MICROBIAL-RIBONUCLEASE']['1BU4']))

matrix, seqs = EIGAs.local_align(protein1, protein2, 6.0)

print(len(seqs))
for seq in seqs:
    print(seq[0])
    print(seq[1])
    print(seq[2])
    print([abs(protein1.fingerprint[i] - protein2.fingerprint[j]) for i, j in zip(seq[1], seq[2])])
    print('')
#pprint.pprint(seqs[0])
#pprint.pprint(seqs[2])
#pprint.pprint(seqs[30])
    