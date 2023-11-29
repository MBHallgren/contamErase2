import os
import sys

p_r = '/home/people/malhal/data/GridION_0011/DTU_2023_1025625_100_SI_SW_887_89_2_93_SE.fq.gz'
p_s = '/home/people/plan/kma_nanomlst/seqs/run44/guppy_5.0.11_sup/barcode03_q10sup.fq.gz'
#extra_s = 'DTU_2023_1026326_100_SI_26_SE.fq'

total = 140000
p_s_limit = 60000

for i in range(1, 11, 1):
    percentage = i / 100
    print (percentage)
    s_size = int(total * percentage)
    if s_size > p_s_limit:
        s_size = p_s_limit
        r_size = (1/percentage) * s_size
    else:
        r_size = total - s_size
    r_size = int(r_size)
    s_size = int(s_size)
    #extra_size = (r_size + s_size)*0.05
    print (s_size, r_size)


    os.system('mkdir subsets/{}'.format(i))
    os.system('seqtk sample -s100 {} {} > subsets/{}/primary.fastq'.format(p_r, r_size, i))
    os.system('seqtk sample -s100 {} {} > subsets/{}/secondary.fastq'.format(p_s, s_size, i))
    #os.system('seqtk sample -s100 {} {} > subsets/{}/extra.fastq'.format(extra_s, extra_size, i))
    os.system('cat subsets/{}/*fastq > subsets/{}/{}_intra.fastq'.format(i, i, i))

