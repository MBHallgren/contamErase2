import os
import sys

infile = '/home/people/malhal/test/test_10_5/rmlst_alignment.res'

top_hits = {}
with open(infile, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            allele = line[0]
            gene = allele.split('_')[0]
            if gene not in top_hits:
                top_hits[gene] = [allele, line[1], 0]
            else:
                if line[1] > top_hits[gene][1]:
                    top_hits[gene] = [allele, line[1], 0]

print (top_hits)

name_file = '/home/people/malhal/test/test_10_5/specie_db.name'
t = 1
with open(name_file, 'r') as f:
    for line in f:
        line = line.strip()
        allele = line
        gene = allele.split('_')[0]
        if gene in top_hits:
            if top_hits[gene][0] == allele:
                top_hits[gene][-1] = t
        t += 1

cmd = ''
for item in top_hits:
    cmd += str(top_hits[item][-1]) + ','
print (cmd[:-1])
