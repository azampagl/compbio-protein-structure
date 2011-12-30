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
from compbio.common.data import LOCAL

import getopt
import csv
import sys

# Get command line args.
opts, args = getopt.getopt(sys.argv[1:], ':o')
if not len(args):
    raise Exception('Missing output directory.')

# A specific order is necessary due to the GoTERM analysis by another group (mostly for visual effects).
ORDER= [
    '1ADF',
    '2JHF',
    '1MGO',
    '1EE2',
    '1QV6',
    '1A71',
    '1N8K',
    '1N92',
    '1P1R',
    '1YE3',
    '3BTO',
    '1QV7',
    '1KUV',
    '1KUY',
    '1KUX',
    '1L0C',
    '1IB1e',
    '1IB1f',
    '1IB1g',
    '1IB1h',
    '1B6B',
    '1B6Bb',
    '1B73',
    '2DWU',
    '1ZUW',
    '1B74',
    '3IST',
    '2JFU',
    '2GZM',
    '3ISV',
    '3HFR',
    '2JFV',
    '2OHO',
    '1EP0',
    '1PM7',
    '1NXM',
    '2IXK',
    '1WLT',
    '1DZR',
    '1RTV',
    '2B9U',
    '1NZC',
    '1NYW',
    '1CLK',
    '1A0C',
    '1MUW',
    '1BXB',
    '2GLK',
    '1QT1',
    '1XIM',
    '1XLA',
    '1DXI',
    '2GYI',
    '1XYL',
    '1S5N',
    '1XYA',
    '3RUBl',
    '8RUC',
    '2VDI',
    '1UZH',
    '1UWA',
    '1UPMs',
    '1SVDm',
    '1RXOc',
    '1IWAb',
    '1IR1s',
    '1GK8i'
]

writer1 = csv.writer(open(args[0] + '/local-align-len.csv', 'wb'))
writer2 = csv.writer(open(args[0] + '/local-align-len-norm.csv', 'wb'))
writer3 = csv.writer(open(args[0] + '/local-align-score.csv', 'wb'))
writer4 = csv.writer(open(args[0] + '/local-align-score-norm.csv', 'wb'))

# write our first header
header = [None]
header.extend(ORDER)
writer1.writerow(header)
writer2.writerow(header)
writer3.writerow(header)
writer4.writerow(header)

# Build a result dictionary.
results1 = {}
results2 = {}
results3 = {}
results4 = {}
for name1 in ORDER:
    results1[name1] = {}
    results2[name1] = {}
    results3[name1] = {}
    results4[name1] = {}
    for name2 in ORDER:
        results1[name1][name2] = None
        results2[name1][name2] = None
        results3[name1][name2] = None
        results4[name1][name2] = None

# Let's spit out the length of the optimal motif.
for name1 in ORDER:
    row1 = [name1]
    row2 = [name1]
    row3 = [name1]
    row4 = [name1]
    
    for name2 in ORDER:
        # See if we processed this combination already
        #if results1[name1][name2] != None:
        #    continue
        print(name1 + ' ' + name2)
        
        protein1 = Protein(ProteinParser.factory('pdb', LOCAL[name1]))
        protein2 = Protein(ProteinParser.factory('pdb', LOCAL[name2]))
        
        # Norm factor will be the smallest protein fingerprint (which is also
        #  the largest possible alignment).
        norm1 = min(len(protein1.fingerprint), len(protein2.fingerprint))
        # The norm for the second set is the value given when two items 
        #  match during local alignment (2) and largest possible alignment.
        norm2 = norm1 * 2
        
        # Local align the two proteins.
        matrix, seqs = EIGAs.local_align(protein1, protein2)
        
        results1[name1][name2] = results1[name2][name1] = len(seqs[0][1])
        results2[name1][name2] = results2[name2][name1] = float(len(seqs[0][1])) / norm1
        results3[name1][name2] = results3[name2][name1] = seqs[0][0].value
        results4[name1][name2] = results4[name2][name1] = float(seqs[0][0].value) / norm2
        
        row1.append(results1[name1][name2])
        row2.append(results2[name1][name2])
        row3.append(results3[name1][name2])
        row4.append(results4[name1][name2])
    
    writer1.writerow(row1)
    writer2.writerow(row2)
    writer3.writerow(row3)
    writer4.writerow(row4)