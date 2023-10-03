import os
import sys

path = '/home/people/malhal/nano_gene_finder/paper/data/nanopore/'
file_list = os.listdir(path)
for item in file_list:
    if item.endswith('.fastq.gz'):
        name = item.split('.')[0]
        cmd = '/home/people/malhal/nanodecon_paper/tmp_tests/NanoDecon/bin/nanodecon' \
              ' --nanopore {} --db_dir /home/people/malhal/nanodecon_paper/tmp_tests/nanodecon_db' \
              ' --ugt 5 --o {} --threads 32'\
            .format(path + item, name)
        os.system(cmd)
