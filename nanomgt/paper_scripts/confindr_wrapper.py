import os
import sys


bf_rates = [1, 2, 3, 4, 5]

path = '/home/people/malhal/papers/rmlst/data/saures/subsets'
#path = '/home/people/malhal/papers/rmlst/data/ecoli/subsets'
for rate in bf_rates:
    for i in range(1, 11, 1):
        os.system(
            'confindr -i {}/{}/ -o {}/{} -d /home/people/malhal/contamErase/benchmarking/confindr/rmlst_db/ --rmlst -b 3 -bf 0.0{} -dt Nanopore -q 20'
            .format(path, i, rate, i, rate))