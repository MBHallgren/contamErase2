import os
import sys
"""
path = '/home/people/malhal/contamErase/benchmarking/confindr/data/'
files = os.listdir(path)

for item in files:
    for i in range(1, 11, 1):
        os.system('confindr -i {}/{}/{}/ -o output_001/{}_{} -d rmlst_db/ --rmlst -b 1 -bf 0.01'.format(path, item, i, item, i))
        
"""
path = '/home/people/malhal/contamErase/benchmarking/confindr/data/nanopore/'
for i in range(1, 11, 1):
    os.system('confindr -i {}/{}/ -o nanopore_output_005/{} -d rmlst_db/ --rmlst -b 1 -bf 0.05 -dt Nanopore -q 14'.format(path, i, i))