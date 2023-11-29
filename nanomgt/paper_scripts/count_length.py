import os
import sys

path = '/home/people/malhal/contamErase/data/illumina/intra/paper/paper_rmlst_profiles/'
files = os.listdir(path)

rmlst_dict = {}

for item in files:
    if item.endswith('.res'):
        with open(path + item, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    gene = line.strip().split('\t')[0].split('_')[0]
                    if gene not in rmlst_dict:
                        rmlst_dict[gene] = set()
                        rmlst_dict[gene].add(line.strip().split('\t')[3])
                    else:
                        rmlst_dict[gene].add(line.strip().split('\t')[3])

for item in rmlst_dict:
    print (item, rmlst_dict[item])
