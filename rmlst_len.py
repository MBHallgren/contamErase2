import os
import sys

rmlst_dict = {}

with open('rmlst.fsa', 'r') as f:
    sequence = ''
    for line in f:
        if line.startswith('>'):
            if sequence != '':
                if str(len(sequence)) not in rmlst_dict[gene]:
                    rmlst_dict[gene][str(len(sequence))] = 0
                rmlst_dict[gene][str(len(sequence))] += 1
                sequence = ''
            line = line.strip()
            gene = line.split('_')[0][1:]
            if gene not in rmlst_dict:
                rmlst_dict[gene] = dict()
            sequence = ''
        else:
            sequence += line.strip()
    if sequence != '':
        if str(len(sequence)) not in rmlst_dict[gene]:
            rmlst_dict[gene][str(len(sequence))] = 0
        rmlst_dict[gene][str(len(sequence))] += 1


for gene in rmlst_dict:
    print (gene)
    if len(rmlst_dict[gene]) > 1:
        for hit in rmlst_dict[gene]:
            print (hit, rmlst_dict[gene][hit])