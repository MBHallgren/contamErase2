import os
import sys
"""
path = '/home/people/malhal/contamErase/benchmarking/confindr/data/'
files = os.listdir(path)

for item in files:
    for i in range(1, 11, 1):
        os.system('mkdir {}/{}/{}'.format(path, item, i))
        os.system('cp /home/people/malhal/contamErase/data/illumina/intra/paper/subsets/isd_1/{}/{}/primary_R1.fastq {}/{}/{}/.'.format(item, i, path, item, i))
        os.system('cp /home/people/malhal/contamErase/data/illumina/intra/paper/subsets/isd_1/{}/{}/primary_R2.fastq {}/{}/{}/.'.format(item, i, path, item, i))
        os.system('cp /home/people/malhal/contamErase/data/illumina/intra/paper/subsets/isd_1/{}/{}/secondary_R1.fastq {}/{}/{}/.'.format(item, i, path, item, i))
        os.system('cp /home/people/malhal/contamErase/data/illumina/intra/paper/subsets/isd_1/{}/{}/secondary_R2.fastq {}/{}/{}/.'.format(item, i, path, item, i))
        os.system('cat {}/{}/{}/primary_R1.fastq {}/{}/{}/secondary_R1.fastq > {}/{}/{}/intra_R1.fastq'.format(path, item, i, path, item, i, path, item, i))
        os.system('cat {}/{}/{}/primary_R2.fastq {}/{}/{}/secondary_R2.fastq > {}/{}/{}/intra_R2.fastq'.format(path, item, i, path, item, i, path, item, i))
        os.system('rm {}/{}/{}/primary_R1.fastq'.format(path, item, i))
        os.system('rm {}/{}/{}/primary_R2.fastq'.format(path, item, i))
        os.system('rm {}/{}/{}/secondary_R1.fastq'.format(path, item, i))
        os.system('rm {}/{}/{}/secondary_R2.fastq'.format(path, item, i))
        
"""

path = '/home/people/malhal/contamErase/benchmarking/confindr/data/nanopore'

for i in range(1, 11, 1):
    os.system('mkdir {}/{}'.format(path, i))
    os.system('cp /home/people/malhal/contamErase/data/nanopore/intra/paper/subsets/{}/primary.fastq {}/{}/.'.format(i, path, i))
    os.system('cp /home/people/malhal/contamErase/data/nanopore/intra/paper/subsets/{}/secondary.fastq {}/{}/.'.format(i, path, i))
    os.system('cat {}/{}/primary.fastq {}/{}/secondary.fastq > {}/{}/intra.fastq'.format(path, i, path, i, path, i))
    os.system('rm {}/{}/primary.fastq'.format(path, i))
    os.system('rm {}/{}/secondary.fastq'.format(path, i))