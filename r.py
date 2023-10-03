import os
import sys

path = '/home/people/malhal/data/GridION_0012/'

files = os.listdir(path)

for item in files:
    if not '101_SI' in item:
        name = item.split('.')[0]
        cmd = 'kmergenetyper -nanopore {} -md 1.9 -t_db ~/nanodecon_paper/tmp_tests/nanodecon_db/rmlst_db -o output/{} -t 32'.format(path + item, name)
        os.system(cmd)