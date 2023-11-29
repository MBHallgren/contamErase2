import os
import sys

file_1 = 'DTU_2023_1025625_100_SI_SW_887_89_2_93_SE_5/rmlst_alignment.res'
file_2 = 'DTU_2023_1025620_100_SI_1_SE_5/rmlst_alignment.res'

list_1 = list()
list_2 = list()

with open (file_1, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            list_1.append(line[0])

with open (file_2, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            list_2.append(line[0])

for item in list_1:
    if item not in list_2:
        print (item)

#for item in list_2:
#    if item not in list_1:
#        print (item)