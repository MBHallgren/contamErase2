import os
import sys

mat_file = 'intra_spe.mat'
gene = '#BACT000036_19'

print_flag = False
with open (mat_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            if line.startswith(gene):
                print_flag = True
            else:
                print_flag = False
        if print_flag:
            print (line.strip())
