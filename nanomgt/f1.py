import os
import sys

path = '/home/people/malhal/contamErase/data/illumina/intra/paper/paper_rmlst_profiles/genomes'
files = os.listdir(path)

for item in files:
    if item != 't.py':
        os.system('mkdir {}'.format(item))
        cmd = 'kma -i {}' \
                ' -o {}/{} -t_db /home/people/malhal/contamErase/test/nanodecon_db/rmlst_db -t 16 -ID 10 -ont -t 16 -md 1.5 -matrix'.format(read_file, item, second_item)
        os.system(cmd)