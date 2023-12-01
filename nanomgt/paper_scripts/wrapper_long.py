import os
import sys

bf_rates = [1, 2, 3, 4, 5]

#path = '/home/people/malhal/contamErase/benchmarking/confindr/data/nanopore/'
#path = '/home/people/malhal/papers/rmlst/data/mix/subsets'
path = '/home/people/malhal/contamErase/benchmarking/confindr/data/nanopore/'
for rate in bf_rates:
    for i in range(1, 11, 1):
        #cmd ='/home/people/malhal/contamErase2/bin/nanomgt --nanopore {}/{}/{}_intra.fastq --o {}/{} --db_dir /home/people/malhal/contamErase_db/ --threads 4 --mrd 0.0{} --bp 2 --pp 2 > {}/{}_{}.csv'\
        #    .format(path, i, i, rate, i, rate, rate, i, 'results')
        #cmd ='/home/people/malhal/contamErase2/bin/nanomgt --nanopore {}/{}/intra.fastq --o {}/{} --db_dir /home/people/malhal/contamErase_db/ --threads 4 --mrd 0.0{} --bp 2 --pp 2 > {}/{}_{}.csv'\
        #    .format(path, i, rate, i, rate, rate, i, 'results')
        #os.system(cmd)
        os.system('mkdir {}/{}'.format(rate, i))
        #os.system('cp /home/people/malhal/papers/rmlst/test/longshot/majority_seqs.fasta {}/{}/majority_seqs.fasta'.format(rate, i))
        os.system('cp /home/people/malhal/test/sarues_test/majority_seqs.fasta {}/{}/majority_seqs.fasta'.format(rate, i))

        os.system('minimap2 -ax map-ont {}/{}/majority_seqs.fasta {}/{}/intra.fastq > {}/{}/output.sam'.format(rate, i, path, i, rate, i))
        os.system('samtools view -bS {}/{}/output.sam > {}/{}/output.bam'.format(rate, i, rate, i))
        os.system('samtools sort {}/{}/output.bam -o {}/{}/output.sorted.bam'.format(rate, i, rate, i))
        os.system('samtools index {}/{}/output.sorted.bam'.format(rate, i))
        os.system('samtools faidx {}/{}/majority_seqs.fasta'.format(rate, i))
        os.system('longshot --bam {}/{}/output.sorted.bam --ref {}/{}/majority_seqs.fasta --out {}/{}/output.vcf --min_mapq 14 --min_alt_frac 0.0{}'.format(rate, i, rate, i, rate, i, rate))