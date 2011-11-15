"""
Common data sets for all tests.

@author Aaron Zampaglione <azampagl@my.fit.edu>
@package EIGAs
@copyright 2011 Aaron Zampaglione
@license MIT
"""
import compbio
import os

# Data dir should always be two levels up.
DIR = os.sep.join((os.path.dirname(compbio.__file__).split(os.sep))[:-2]) + os.sep + 'data'

# Hard alignment proteins.
HARD_EXT = 'pdb'
HARD_DIR = os.path.abspath(DIR + os.sep + 'hard' + os.sep + HARD_EXT + os.sep)
HARD = {}
for f in os.listdir(HARD_DIR):
    HARD[f.rsplit('.', 1)[0]] = (f.rsplit('.', 1)[0], os.path.abspath(HARD_DIR + os.sep + f), None)

# Special cases for hard proteins.
HARD['1FXIa'] = ('1FXIa', HARD['1FXI'][1], 'A')
del HARD['1FXI']

HARD['3HHRb'] = ('3HHRb', HARD['3HHR'][1], 'B')
del HARD['3HHR']

HARD['3HLAb'] = ('3HLAb', HARD['3HLA'][1], 'B')
del HARD['3HLA']

HARD['2AZAa'] = ('2AZAa', HARD['2AZA'][1], 'A')
del HARD['2AZA']

HARD['1CEWi'] = ('1CEWi', HARD['1CEW'][1], 'I')
del HARD['1CEW']

HARD['1MOLa'] = ('1MOLa', HARD['1MOL'][1], 'A')
del HARD['1MOL']

HARD['1NSBa'] = ('1NSBa', HARD['1NSB'][1], 'A')
del HARD['1NSB']

HARD['1BGEb'] = ('1BGEb', HARD['1BGE'][1], 'B')
del HARD['1BGE']

HARD['2GMFa'] = ('2GMFa', HARD['2GMF'][1], 'A')
del HARD['2GMF']

#
#
#

# Skolnick data set.
SKOLNICK_EXT = 'pdb'
SKOLNICK_DIR = os.path.abspath(DIR + os.sep + 'skolnick' + os.sep + HARD_EXT + os.sep)
SKOLNICK = {}
for d in os.listdir(SKOLNICK_DIR):
    SKOLNICK[d.upper()] = {}
    for f in os.listdir(SKOLNICK_DIR + os.sep + d):
        SKOLNICK[d.upper()][f.rsplit('.', 1)[0]] = (f.rsplit('.', 1)[0], os.path.abspath(SKOLNICK_DIR + os.sep + d + os.sep +f), None)

# Special cases for skolnick.
SKOLNICK['FLAVODXIN-LIKE']['1QMPA'] = ('1QMPA', SKOLNICK['FLAVODXIN-LIKE']['1QMP'][1], 'A')
SKOLNICK['FLAVODXIN-LIKE']['1QMPB'] = ('1QMPB', SKOLNICK['FLAVODXIN-LIKE']['1QMP'][1], 'B')
SKOLNICK['FLAVODXIN-LIKE']['1QMPC'] = ('1QMPC', SKOLNICK['FLAVODXIN-LIKE']['1QMP'][1], 'C')
SKOLNICK['FLAVODXIN-LIKE']['1QMPD'] = ('1QMPD', SKOLNICK['FLAVODXIN-LIKE']['1QMP'][1], 'D')
del SKOLNICK['FLAVODXIN-LIKE']['1QMP']

SKOLNICK['FLAVODXIN-LIKE']['4TMYA'] = ('4TMYA', SKOLNICK['FLAVODXIN-LIKE']['4TMY'][1], 'A')
SKOLNICK['FLAVODXIN-LIKE']['4TMYB'] = ('4TMYB', SKOLNICK['FLAVODXIN-LIKE']['4TMY'][1], 'B')
del SKOLNICK['FLAVODXIN-LIKE']['4TMY']

SKOLNICK['CUPREDOXIN-LIKE']['1BYOA'] = ('1BYOA', SKOLNICK['CUPREDOXIN-LIKE']['1BYO'][1], 'A')
SKOLNICK['CUPREDOXIN-LIKE']['1BYOB'] = ('1BYOB', SKOLNICK['CUPREDOXIN-LIKE']['1BYO'][1], 'B')
del SKOLNICK['CUPREDOXIN-LIKE']['1BYO']

SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1A'] = ('1RN1A', SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1'][1], 'A')
SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1B'] = ('1RN1B', SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1'][1], 'B')
SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1C'] = ('1RN1C', SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1'][1], 'C')
del SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1']