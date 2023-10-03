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
                rmlst_dict[species] = dict()
                rmlst_dict[species][line[0]] = line[1:54]
            else:
                rmlst_dict[species][line[0]] = line[1:54]

with open ('rmlst_headers.txt', 'w') as f:
    f.write('\n'.join(headers))
with open('rmlst_scheme.json', 'w') as f:
    json.dump(rmlst_dict, f)