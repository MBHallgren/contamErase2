import os
import sys

path = '/home/people/malhal/contamErase/benchmarking/confindr/data/nanopore/'
files = os.listdir(path)

for item in files:
    if item != 't.py':
        os.system('mkdir {}'.format(item))
        read_file = "{}/{}/intra.fastq".format(path, item)
        cmd = '/home/people/malhal/contamErase2/bin/nanomgt --nanopore {}' \
            ' --o {} --db_dir /home/people/malhal/contamErase/test/nanodecon_db/ --threads 16 --isd 0.01'.format(read_file, item)
        os.system(cmd)