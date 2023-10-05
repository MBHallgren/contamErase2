import os
import sys
import hashlib

with open('rmlst_hashes.txt', 'w') as out:
    with open('rmlst.fsa', 'r') as f:
        sequence = ''
        for line in f:
            if line.startswith('>'):
                if sequence != '':
                    out.write('{}\t{}\n'.format(hashlib.md5(sequence.encode('utf-8')).hexdigest(), header))
                sequence = ''
                header = line.strip()[1:]
            else:
                sequence += line.strip()
        out.write('{}\t{}\n'.format(hashlib.md5(sequence.encode('utf-8')).hexdigest(), header))
