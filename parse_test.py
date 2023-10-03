import os
import sys

"""

file = 'rmlst_escherichia_coli.mat'
count = 0
with open(file, 'r') as f:
    flag = False
    for line in f:
        if line.startswith('#'):
            if line.rstrip() == '#BACT000003_17':
                flag = True
            else:
                flag = False
        if flag:
            print (line.strip())
            if line.startswith('A') or line.startswith('C') or line.startswith('G') or line.startswith('T'):
                count += 1
            elif line.startswith('N') or line.startswith('-'):
                print (count)
print (count)



file = 'rmlst_escherichia_coli.aln'
with open(file, 'r') as f:
    flag = False
    for line in f:
        if line.startswith('# BACT000003_17'):
            flag = True
        elif line.startswith('#'):
            flag = False
        if flag:
            print (line.strip())
"""

s = 'ATGGGTCAGAAAGTACATCCTAATGGTATTCGCCTGGGTATTGTAAAACCATGGAACTCTACCTGGTTTGCGAACACCAAAGAATTCGCTGACAACCTGGACAGCGATTTTAAAGTACGTCAGTACCTGACTAAGGAACTGGCTAAAGCGTCCGTATCTCGTATCGTTATCGAGCGTCCGGCTAAGAGCATCCGTGTAACCATTCACACTGCTCGCCCGGGTATCGTTATCGGTAAAAAAGGTGAAGACGTAGAAAAACTGCGTAAGGTCGTAGCGGACATCGCTGGCGTTCCTGCACAGATCAACATCGCCGAAGTTCGTAAGCCTGAACTGGACGCAAAACTGGTTGCTGACAGCATCACTTCTCAGCTGGAACGTCGCGTTATGTTCCGTCGTGCTATGAAGCGTGCTGTACAGAACGCAATGCGTCTGGGCGCTAAAGGTATTAAAGTTGAAGTTAGCGGCCGTCTGGGCGGCGCGGAAATCGCACGTACCGAATGGTACCGCGAAGGTCGCGTACCGCTGCACACTCTGCGTGCTGACATCGACTACAACACCTCTGAAGCGCACACCACTTACGGTGTAATCGGCGTTAAAGTGTGGATCTTCAAAGGCGAGATCCTGGGTGGTATGGCTGCTGTTGAACAACCGGAAAAACCGGCTGCTCAGCCTAAAAAGCAGCAGCGTAAAGGCCGTAAATAA'

print (s[245:255])