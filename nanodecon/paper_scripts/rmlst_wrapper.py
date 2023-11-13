import os
import sys

bf_rates = [1, 2, 3, 4, 5]

path = '/home/people/malhal/contamErase/benchmarking/confindr/data/nanopore/'
for rate in bf_rates:
    for i in range(1, 11, 1):
        os.system(
            '/home/people/malhal/contamErase2/bin/nanodecon --nanopore {}/{}/ --o {}/{} --db_dir /home/people/malhal/contamErase_db/ --rmlst --threads 16 --mrd 0.05 '
            .format(path, i, rate, i, rate))