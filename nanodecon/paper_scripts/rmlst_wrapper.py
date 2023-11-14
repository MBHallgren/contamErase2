import os
import sys

bf_rates = [1, 2, 3, 4, 5]

path = '/home/people/malhal/contamErase/benchmarking/confindr/data/nanopore/'
for rate in bf_rates:
    for i in range(1, 11, 1):
        cmd ='/home/people/malhal/contamErase2/bin/nanodecon --nanopore {}/{}/intra.fastq --o {}/{} --db_dir /home/people/malhal/contamErase_db/ --threads 16 --mrd 0.0{} > {}/{}_{}.csv'\
            .format(path, i, rate, i, rate, rate, i, 'results')
        print (cmd)
        os.system(cmd)