import os
import sys

path = '/home/people/malhal/data/nanodecon/illumina/'
file_list = os.listdir(path)
for i in range(0, len(file_list), 2):
    name = file_list[i].split('_')[0]
    cmd = '/home/people/malhal/nanodecon_paper/tmp_tests/NanoDecon/bin/nanodecon' \
          ' --illumina {} {} --db_dir /home/people/malhal/nanodecon_paper/tmp_tests/nanodecon_db' \
          ' --o {} --threads 16 >> nanodecon_output.txt'\
        .format(path + file_list[i], path + file_list[i+1], name)
    os.system(cmd)