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
    HARD[f.rsplit('.', 1)[0]] = (os.path.abspath(HARD_DIR + os.sep + f), None)

# Skolnick data set.
SKOLNICK_EXT = 'pdb'
SKOLNICK_DIR = os.path.abspath(DIR + os.sep + 'skolnick' + os.sep + HARD_EXT + os.sep)
SKOLNICK = {}
for d in os.listdir(SKOLNICK_DIR):
    SKOLNICK[d.upper()] = {}
    for f in os.listdir(SKOLNICK_DIR + os.sep + d):
        SKOLNICK[d.upper()][f.rsplit('.', 1)[0]] = (os.path.abspath(SKOLNICK_DIR + os.sep + d + os.sep +f), None)

# Special cases for skolnick.
SKOLNICK['FLAVODXIN-LIKE']['1QMPA'] = (SKOLNICK['FLAVODXIN-LIKE']['1QMP'][0], 'A')
SKOLNICK['FLAVODXIN-LIKE']['1QMPB'] = (SKOLNICK['FLAVODXIN-LIKE']['1QMP'][0], 'B')
SKOLNICK['FLAVODXIN-LIKE']['1QMPC'] = (SKOLNICK['FLAVODXIN-LIKE']['1QMP'][0], 'C')
SKOLNICK['FLAVODXIN-LIKE']['1QMPD'] = (SKOLNICK['FLAVODXIN-LIKE']['1QMP'][0], 'D')
del SKOLNICK['FLAVODXIN-LIKE']['1QMP']

SKOLNICK['FLAVODXIN-LIKE']['4TMYA'] = (SKOLNICK['FLAVODXIN-LIKE']['4TMY'][0], 'A')
SKOLNICK['FLAVODXIN-LIKE']['4TMYB'] = (SKOLNICK['FLAVODXIN-LIKE']['4TMY'][0], 'B')
del SKOLNICK['FLAVODXIN-LIKE']['4TMY']

SKOLNICK['CUPREDOXIN-LIKE']['1BYOA'] = (SKOLNICK['CUPREDOXIN-LIKE']['1BYO'][0], 'A')
SKOLNICK['CUPREDOXIN-LIKE']['1BYOB'] = (SKOLNICK['CUPREDOXIN-LIKE']['1BYO'][0], 'B')
del SKOLNICK['CUPREDOXIN-LIKE']['1BYO']

SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1A'] = (SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1'][0], 'A')
SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1B'] = (SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1'][0], 'B')
SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1C'] = (SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1'][0], 'C')
del SKOLNICK['MICROBIAL_RIBONUCLEASE']['1RN1']