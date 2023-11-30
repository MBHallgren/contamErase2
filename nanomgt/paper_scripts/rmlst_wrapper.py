import os
import sys

bf_rates = [1, 2, 3, 4, 5]

#path = '/home/people/malhal/contamErase/benchmarking/confindr/data/nanopore/'
path = '/home/people/malhal/papers/rmlst/data/mix/subsets'
#path = '/home/people/malhal/contamErase/benchmarking/confindr/data/nanopore/'
for rate in bf_rates:
    for i in range(1, 11, 1):
        cmd ='/home/people/malhal/contamErase2/bin/nanomgt --nanopore {}/{}/{}_intra.fastq --o {}/{} --db_dir /home/people/malhal/contamErase_db/ --threads 4 --mrd 0.0{} > {}/{}_{}.csv'\
            .format(path, i, i, rate, i, rate, rate, i, 'results')
        #cmd ='/home/people/malhal/contamErase2/bin/nanomgt --nanopore {}/{}/intra.fastq --o {}/{} --db_dir /home/people/malhal/contamErase_db/ --threads 4 --mrd 0.0{} > {}/{}_{}.csv'\
        #    .format(path, i, rate, i, rate, rate, i, 'results')
        os.system(cmd)