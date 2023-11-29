import os
import sys

read_file = '/home/people/malhal/data/ecoli_results/test/36.txt'
read_set = set()

with open(read_file, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        read_set.add(line[-1])

#frag_file_list = ['/home/people/malhal/contamErase/test/trial_5_10/bacteria_alignment.frag']
                  #'/home/people/malhal/contamErase/test/trial_5_10/viral_alignment.frag',
                  #'/home/people/malhal/contamErase/test/trial_5_10/human_alignment.frag',
                  #'/home/people/malhal/contamErase/test/trial_5_10/plasmid_alignment.frag']

frag_file_list = ['/home/people/malhal/contamErase/test/trial_5_10/tmp_bacteria_alignment.frag']


#frag_file_list = ['/home/people/malhal/contamErase/test/trial_5_10/bacteria_alignment.frag',
#                  '/home/people/malhal/contamErase/test/trial_5_10/viral_alignment.frag',
#                  '/home/people/malhal/contamErase/test/trial_5_10/human_alignment.frag',
#                  '/home/people/malhal/contamErase/test/trial_5_10/plasmid_alignment.frag']

for frag_file in frag_file_list:
    print (frag_file)
    all_frags = dict()

    with open(frag_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            all_frags[line[-1]] = line[-2]

    for item in read_set:
        if item in all_frags:
            print (item, all_frags[item])
        else:
            print (item, 'NOT FOUND IN FRAG FILE')

