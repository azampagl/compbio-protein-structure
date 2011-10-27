"""
Scapes PDB database for proteus 300 data set.

@see http://www.irisa.fr/symbiose/old/softwares/resources/p_300_set
@see http://astral.berkeley.edu/pdbstyle-1.75.html
@see ./README.txt

@author Aaron Zampaglione <azampagl@my.fit.edu>
@copyright 2011 (c) Aaron Zampaglione
@license MIT
"""
import os
import re
import sys
import time
import urllib

BASE_URL = 'http://astral.berkeley.edu/pdbstyle.cgi?id='
EXT = 'ent'

readme = open('README.txt', 'r')
line = readme.readline()
for match in re.finditer(r"<td>(?P<domains>.*?)(\s<br />)?</td>.*?<td>(?P<fold>.*?)(\s<br />)?</td>.*?<td>(?P<family>.*?)(\s<br />)?</td>", line, re.I | re.S):
    # Clean the fold.
    fold = match.group('fold').lower().strip() \
        .replace(' ', '_') \
        .replace('/', '_') \
        .replace('(', '_') \
        .replace(')', '_')
    
    # Create the fold directory.
    if not os.path.exists(EXT + '/' + fold):
        os.mkdir(EXT + '/' + fold)
    
     # Clean the domains.
    domains = match.group('domains').lower().strip().split(', ')
    for domain in domains:
        result = urllib.urlopen(BASE_URL + domain)
        ent = open(EXT + '/' + fold + '/' + domain + '.' + EXT, 'w')
        ent.write(result.read())
        ent.close()
        time.sleep(2)

readme.close()