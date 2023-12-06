import os
import sys

ecoli_file = '/home/people/malhal/gridion_and_coli.txt'
ecoli_list = list()
for line in open(ecoli_file, 'r'):
    line = line.strip().split('\t')
    name = line[0].split('_')[0] + '_' + line[0].split('_')[1] + '_' + line[0].split('_')[2]
    ecoli_list.append(name)

path = '/home/people/malhal/data/GridION_0010'
files = os.listdir(path)
for file in files:
    name = file.split('.')[0]
    name = name.split('_')[0] + '_' + name.split('_')[1] + '_' + name.split('_')[2]
    if name in ecoli_list:
        print (name + ' IS ECOLI in path: ' + path)

path = '/home/people/malhal/data/GridION_0011'
files = os.listdir(path)
for file in files:
    name = file.split('.')[0]
    name = name.split('_')[0] + '_' + name.split('_')[1] + '_' + name.split('_')[2]
    if name in ecoli_list:
        print (name + ' IS ECOLI in path: ' + path)

path = '/home/people/malhal/data/GridION_0012'
files = os.listdir(path)
for file in files:
    name = file.split('.')[0]
    name = name.split('_')[0] + '_' + name.split('_')[1] + '_' + name.split('_')[2]
    if name in ecoli_list:
        print (name + ' IS ECOLI in path: ' + path)


path = '/home/people/malhal/data/GridION_0013/Raw_data'
files = os.listdir(path)
for file in files:
    name = file.split('.')[0]

    name = name.split('_')[0] + '_' + name.split('_')[1] + '_' + name.split('_')[2]
    if name in ecoli_list:
        print (name + ' IS ECOLI in path: ' + path)

