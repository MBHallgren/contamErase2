import os
import sys

def find_mutations(string_1, string_2):
    for i in range(len(string_1)):
        if string_1[i].upper() != string_2[i].upper():
            print ('{}_{}_{}'.format(i+1, string_1[i], string_2[i]))

profile_1 = 'DTU_2023_1028222.fsa'
profile_2 = 'DTU_2023_1028247.fsa'

profile_1_dict = {}
profile_2_dict = {}

with open(profile_1, 'r') as f:
    for line in f:
        if line.startswith('>'):
            header = line[1:].strip()
            profile_1_dict[header] = ''
        else:
            profile_1_dict[header] += line.strip()

with open(profile_2, 'r') as f:
    for line in f:
        if line.startswith('>'):
            header = line[1:].strip()
            profile_2_dict[header] = ''
        else:
            profile_2_dict[header] += line.strip()

for item in profile_1_dict:
    gene = item.split('_')[0]
    for item2 in profile_2_dict:
        gene2 = item2.split('_')[0]
        if gene == gene2:
            print ('{}>{}'.format(item, item2))
            find_mutations(profile_1_dict[item], profile_2_dict[item2])
