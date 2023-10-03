import sys
import json

infile = 'rmlst_scheme.txt'

rmlst_dict = {}

with open(infile, 'r') as f:
    headers = f.readline().strip().split('\t')
    for line in f:
        line = line.strip().split('\t')
        if line != ['']:
            species = line[55]
            if species not in rmlst_dict:
                rmlst_dict[species] = set()
                for item in line[1:54]:
                    if item != '' or item != 'N':
                        rmlst_dict[species].add(item)
            else:
                for item in line[1:54]:
                    if item != '' or item != 'N':
                        rmlst_dict[species].add(item)
for item in rmlst_dict:
    if item == 'Escherichia coli':
        print (item, len(rmlst_dict[item]))

print (len(rmlst_dict))