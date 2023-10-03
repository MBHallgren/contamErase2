import os
import sys

bac_dict = {}
with open('Escherichia_coli_read_pool_R1.frag', 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        if line[-2] == 'BACT000065_1182':
            if 'BACT000065_1182' not in bac_dict:
                bac_dict['BACT000065_1182'] = [line[-1]]
            else:
                bac_dict['BACT000065_1182'].append(line[-1])
        if line[-2] == 'BACT000065_10':
            if 'BACT000065_10' not in bac_dict:
                bac_dict['BACT000065_10'] = [line[-1]]
            else:
                bac_dict['BACT000065_10'].append(line[-1])

print (len(bac_dict['BACT000065_1182']))
print (len(bac_dict['BACT000065_10']))

print (len(set(bac_dict['BACT000065_1182']).intersection(set(bac_dict['BACT000065_10']))))
