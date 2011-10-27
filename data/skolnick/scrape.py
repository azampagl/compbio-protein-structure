"""
Scapes PDB database for skolnick data set.

@see http://astral.berkeley.edu/pdbstyle-1.75.html
@see ./README.txt

@author Aaron Zampaglione <azampagl@my.fit.edu>
@copyright 2011 (c) Aaron Zampaglione
@license MIT
"""
import os
import sys
import time
import urllib

BASE_URL = 'http://astral.berkeley.edu/pdbstyle.cgi?id='
EXT = 'ent'

# Skolnick dictionary.
skolnick = {}

#
# Parse the readme to find fold families and proteins.
#
last_fold = ''

readme = open('README.txt', 'r')
for line in readme:
    if line == "\n":
        continue
    
    line = line[:-1].split(" ")
    
    # Line contains the fold name.
    if len(line) == 1:
        last_fold = line[0]
        skolnick[line[0]] = []
        continue
    
    # Append the protein to the last fold discovered.
    skolnick[last_fold].append(line[1])
readme.close()

# Create directories.
for fold in skolnick.keys():
    if os.path.exists(EXT + '/' + fold):
        os.removedirs(EXT + '/' + fold)
    os.mkdir(EXT + '/' + fold)

for fold, proteins in skolnick.items():
    for protein in proteins:
        result = urllib.urlopen(BASE_URL + protein)
        ent = open(EXT + '/' + fold + '/' + protein + '.' + EXT, 'w')
        ent.write(result.read())
        ent.close()
        time.sleep(2)