import os
import sys

path = '/home/people/malhal/data/nanodecon_data/'

black_list = set()

for item in os.listdir(path):
    if not item.startswith('.'):
        id = item[0:16]
        if id not in black_list:
            os.system('cat {}*_R1* >> ecoli/{}_R1.fastq'.format(path+id, id))
            os.system('cat {}*_R2* >> ecoli/{}_R2.fastq'.format(path+id, id))
            black_list.add(id)

path = '/home/people/malhal/data/ecoli/'

files = os.listdir(path)

for i in range(0, len(files), 2):
    name = files[i].split('.')[0]
    cmd = 'kmergenetyper -illumina {} {} -md 5 -t 16 -o ecoli_results/{} -t_db ~/contamErase/test/nanodecon_db/rmlst_db -keep'.format(path+files[i], path+files[i+1], name)
    os.system(cmd)