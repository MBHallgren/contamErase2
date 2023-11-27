import os
import sys

path = '/home/people/plan/kma_nanomlst/seqs/run44/guppy_5.0.11_sup'
files = os.listdir(path)

for item in files:
    name = item.split('_')[0]
    cmd = '~/contamErase2/bin/nanodecon --nanopore {}/{} --o {}_5 --db_dir /home/people/malhal/contamErase_db/ --q_score 14 --mrd 0.05 --threads 16 > {}_5.csv'.format(path, item, name, name)
    print (cmd)
    #os.system(cmd)