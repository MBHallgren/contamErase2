import os
import sys

path = '/home/people/malhal/contamErase/data/nanopore/intra/paper/subsets/'
files = os.listdir(path)

for item in files:
    if item != 't.py':
        os.system('mkdir {}'.format(item))
        read_file = "{}/{}/{}_intra.fastq".format(path, item, item)
        cmd = '/home/people/malhal/contamErase2/bin/nanodecon --nanopore {}' \
            ' --o {} --db_dir /home/people/malhal/contamErase/test/nanodecon_db/ --threads 16 --isd 0.01'.format(read_file, item)
        os.system(cmd)