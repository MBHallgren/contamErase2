import os
import sys

p_r1 = '/home/people/malhal/data/008253_R1.fastq'
p_r2 = '/home/people/malhal/data/008253_R2.fastq'
p_s1 = '/home/people/malhal/data/026027_R1.fastq'
p_s2 = '/home/people/malhal/data/026027_R2.fastq'
#extra_s1 = '/home/people/malhal/contamErase/data/illumina/intra/paper/1000973_R1.fastq'
#extra_s2 = '/home/people/malhal/contamErase/data/illumina/intra/paper/1000973_R2.fastq'

total = 2000000

for i in range(1, 11, 1):
    os.system('mkdir sim_subsets/{}'.format(i))
    percentage = i / 100
    extra_size = int(0.1*total)
    print (percentage)
    s_size = int(total * percentage)
    r_size = total - s_size

    os.system('mkdir subsets/{}'.format(i))
    os.system('seqtk sample -s100 {} {} > sim_subsets/{}/primary_R1.fastq'.format(p_r1, r_size, i))
    os.system('seqtk sample -s100 {} {} > sim_subsets/{}/primary_R2.fastq'.format(p_r2, r_size, i))
    os.system('seqtk sample -s100 {} {} > sim_subsets/{}/secondary_R1.fastq'.format(p_s1, s_size, i))
    os.system('seqtk sample -s100 {} {} > sim_subsets/{}/secondary_R2.fastq'.format(p_s2, s_size, i))
    #os.system('seqtk sample -s100 {} {} > subsets/{}/extra_R1.fastq'.format(extra_s1, extra_size, i))
    #os.system('seqtk sample -s100 {} {} > subsets/{}/extra_R2.fastq'.format(extra_s2, extra_size, i))
    os.system('cat sim_subsets/{0}/*_R1.fastq > sim_subsets/{0}/{0}_R1.fastq'.format(i))
    os.system('cat sim_subsets/{0}/*_R2.fastq > sim_subsets/{0}/{0}_R2.fastq'.format(i))


