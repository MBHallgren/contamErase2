import os
import sys

file_1 = '/home/people/malhal/data/ecoli_results/DTU_2022_1013323_R1/DTU_2022_1013323_R1.res'
file_2 = '/home/people/malhal/data/ecoli_results/DTU_2022_1023536_R1/DTU_2022_1023536_R1.res'

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

for item in list_2:
    if item not in list_1:
        print (item)