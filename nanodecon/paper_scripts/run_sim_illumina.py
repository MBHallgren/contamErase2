import os
import sys

path = '/home/people/malhal/contamErase/data/illumina/intra/paper/sim_subsets/'
files = os.listdir(path)

os.system('mkdir results/sim/isd_5')

for id in files:
    total_path = path + '/' + id
    print (total_path)
    cmd = '/home/people/malhal/nanodecon_paper/tmp_tests/NanoDecon/bin/nanodecon --illumina {}/{}_R1.fastq {}/{}_R2.fastq' \
          ' --o results/sim/isd_5/{} --db_dir /home/people/malhal/contamErase/test/nanodecon_db/ --threads 8 --isd 0.01'.format(total_path, id, total_path, id, id)
    #print (cmd)
    os.system(cmd)