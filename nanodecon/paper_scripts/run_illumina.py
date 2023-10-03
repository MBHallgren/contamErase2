import os
import sys

path = '/home/people/malhal/contamErase/data/illumina/intra/paper/subsets/isd_1/'
files = os.listdir(path)

for item in files:
    trial_path = path + item
    for id in os.listdir(trial_path):
        total_path = trial_path + '/' + id
        print (total_path)
        cmd = '/home/people/malhal/nanodecon_paper/tmp_tests/NanoDecon/bin/nanodecon --illumina {}/{}_R1.fastq {}/{}_R2.fastq' \
              ' --o results/{}/isd_5/{} --db_dir /home/people/malhal/contamErase/test/nanodecon_db/ --threads 8 --isd 0.05'.format(total_path, id, total_path, id, item, id)
        #print (cmd)
        os.system(cmd)