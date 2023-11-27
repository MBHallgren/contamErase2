import os
import sys

"""
path = '/home/people/malhal/contamErase/benchmarking/confindr/data/'
files = os.listdir(path)

for item in files:
    for i in range(1, 11, 1):
        os.system('confindr -i {}/{}/{}/ -o output_001/{}_{} -d rmlst_db/ --rmlst -b 1 -bf 0.01'.format(path, item, i, item, i))

"""

bf_rates = [1, 2, 3, 4, 5]

path = '/home/people/malhal/contamErase/benchmarking/confindr/data/nanopore/'
for rate in bf_rates:
    for i in range(1, 11, 1):
        os.system(
            'confindr -i {}/{}/ -o {}/{} -d /home/people/malhal/contamErase/benchmarking/confindr/rmlst_db/ --rmlst -b 3 -bf 0.0{} -dt Nanopore -q 14'
            .format(path, i, rate, i, rate))