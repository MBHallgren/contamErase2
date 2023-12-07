import os
import sys

#ecoli
prim = '/home/people/malhal/papers/rmlst/test/gridion_0011/DTU_2023_1025625_100_SI_SW_887_89_2_93_SE_5/rmlst_alignment.fsa'
#sec = '/home/people/malhal/papers/rmlst/test/gridion_0011/DTU_2023_1025624_100_SI_371_SE_5/rmlst_alignment.fsa'
#sec = '/home/people/malhal/data/test/test/rmlst_alignment.fsa'
#sec = '/home/people/malhal/papers/rmlst/test/44_test/barcode03_5/rmlst_alignment.fsa'

#sarues
#prim = '/home/people/malhal/papers/rmlst/benchmarking/sample_benchmarks/1028222.fsa'
#sec = '/home/people/malhal/papers/rmlst/benchmarking/sample_benchmarks/1028247.fsa'

gene_dict_1 = {}
gene_dict_2 = {}

with open(prim, 'r') as f:
    for line in f:
        if line.startswith('>'):
            header = line[1:].strip()
            gene = header.split('_')[0]
            gene_dict_1[gene] = ''
        else:
            gene_dict_1[gene] += line.strip()

with open(sec, 'r') as f:
    for line in f:
        if line.startswith('>'):
            header = line[1:].strip()
            gene = header.split('_')[0]
            gene_dict_2[gene] = ''
        else:
            gene_dict_2[gene] += line.strip()



for key in gene_dict_1:
    seq_1 = gene_dict_1[key]
    seq_2 = gene_dict_2[key]
    print (key)
    #print (key, len(seq_1), len(seq_2))
    #if len(seq_1) != len(seq_2):
    #    print ('lengths differ')
    for i in range(len(seq_1)):
        if seq_1[i].upper() != seq_2[i].upper():
            print ('{}_{}_{}'.format(i+1, seq_1[i], seq_2[i]))