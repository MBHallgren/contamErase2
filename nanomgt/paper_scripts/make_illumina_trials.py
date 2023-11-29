"""
DTU_2023_1017265_R1 and DTU_2022_1013180_R1: 48/53, Different: BACT000060, BACT000064, BACT000049, BACT000036, BACT000006
From Pakistan and Romania - FINDS 51 Correct mlst genes,

DTU_2020_1000585_R1 and DTU_2022_1013180_R1: 50/53, Different: BACT000060, BACT000006, BACT000001
From Saudi Arabia and Romania NO RMLST FOUND

DTU_2020_1000585_R1 and DTU_2022_1013323_R1: 38/53 1, 2, 5, 11, 19, 30, 35, 40, 44, 48, 49, 53, 59, 62, 64
From Saudi Arabia and Portugal - RERUN WITH ID 0.03

DTU_2022_1013323_R1 and DTU_2023_1018673_R1: 27/53 1, 2, 5, 9, 10, 11, 12, 20, 30, 32, 33, 35, 36, 38, 40, 45, 46, 48, 49, 50, 53, 57, 59, 60, 62, 64
From portugal and Ghana, Does not work well

DTU_2022_1013323_R1 and DTU_2022_1023536: 46/53 1, 9, 30, 40, 48, 53, 59, 60, 64
From Portugal and Belgium

DTU_2020_1000822 and DTU_2022_1023536: 32/53 1, 2, 5, 11, 12, 16, 20 ,30, 32, 33, 34, 35, 36, 38, 40, 45, 46, 49, 53, 60, 62
From Hong Kong and Belgium

TRIAL 7 WITH US/EU DATA



Trial 8 with simulated data - i.e. extreme bias towards database!




LIST:
DTU_2020_1000395_1_SI_GHA_SEK_018	1000395	1	GHA Ghana	Nextseq141
DTU_2020_1000585_1_SI_SAU_A1Q_016	1000585	1	SAU Saudi Arabia	NextSeq_144
DTU_2020_1000780_1_SI_HKG_SHA_013	1000780	1	HKG Hong Kong	NextSeq_146
DTU_2020_1000789_1_SI_HKG_SHA_022	1000789	1	HKG Hong Kong	NextSeq_146
DTU_2020_1000821_1_SI_HKG_SHA_054	1000821	1	HKG Hong Kong	NextSeq_146
DTU_2020_1000822_1_SI_HKG_SHA_055	1000822	1	HKG Hong Kong	NextSeq_146
DTU_2022_1013147_1_SI_M2020_10063258	1013147	1	HUN Hungary	NextSeq_197
DTU_2022_1013169_1_SI_12800_2IS	1013169	1	ROU Romania	NextSeq_197
DTU_2022_1013176_1_SI_13226_2DB	1013176	1	ROU Romania	NextSeq_197
DTU_2022_1013180_1_SI_13446_2PH	1013180	1	ROU Romania	NextSeq_197
DTU_2022_1013321_1_SI_PAT_20_24568_ECX	1013321	1	PRT Portugal	NextSeq_205
DTU_2022_1013323_1_SI_PAT_20_25369_ECX	1013323	1	PRT Portugal	NextSeq_205
DTU_2022_1016570_1_SI_Ec_Matrix_EQAsia_22_A	1016570	1	NULL	NextSeq_213
DTU_2022_1017218_1_SI_TWIW_02_PAK_PES_AA_068	1017218	1	PAK Pakistan	NextSeq_234
DTU_2022_1023536_1_SI_S21FP01602	1023536	1	BEL Belgium	NextSeq_230
DTU_2023_1017265_1_SI_TWIW_02_PAK_MUL_AC_045	1017265	1	PAK Pakistan	NextSeq_242
DTU_2023_1017300_1_SI_TWIW_02_HRV_ZAG_BL_010	1017300	1	HRV Croatia	NextSeq_242
DTU_2023_1018673_1_SI_TWIW_02_GHA_SEK_BJ_053	1018673	1	GHA Ghana	NextSeq_261
DTU_2023_1018913_1_SI_TWIW_02_NGA_LAF_AN_013	1018913	1	NGA Nigeria	NextSeq_262


"""

import os
import sys

#Trial 1
p_r1 = '/home/people/malhal/contamErase/data/illumina/intra/paper/DTU_2020_1000822_R1.fastq'
p_r2 = '/home/people/malhal/contamErase/data/illumina/intra/paper/DTU_2020_1000822_R2.fastq'
p_s1 = '/home/people/malhal/contamErase/data/illumina/intra/paper/DTU_2022_1023536_R1.fastq'
p_s2 = '/home/people/malhal/contamErase/data/illumina/intra/paper/DTU_2022_1023536_R2.fastq'
extra_s1 = '/home/people/malhal/contamErase/data/illumina/intra/paper/1000973_R1.fastq'
extra_s2 = '/home/people/malhal/contamErase/data/illumina/intra/paper/1000973_R2.fastq'

total = 2000000

os.system('mkdir subsets/isd_1/trial_6')
for i in range(1, 11, 1):
    percentage = i / 100
    extra_size = int(0.1*total)
    print (percentage)
    s_size = int(total * percentage)
    r_size = total - s_size

    os.system('mkdir subsets/isd_1/trial_6/{}'.format(i))
    os.system('seqtk sample -s100 {} {} > subsets/isd_1/trial_6/{}/primary_R1.fastq'.format(p_r1, r_size, i))
    os.system('seqtk sample -s100 {} {} > subsets/isd_1/trial_6/{}/primary_R2.fastq'.format(p_r2, r_size, i))
    os.system('seqtk sample -s100 {} {} > subsets/isd_1/trial_6/{}/secondary_R1.fastq'.format(p_s1, s_size, i))
    os.system('seqtk sample -s100 {} {} > subsets/isd_1/trial_6/{}/secondary_R2.fastq'.format(p_s2, s_size, i))
    os.system('seqtk sample -s100 {} {} > subsets/isd_1/trial_6/{}/extra_R1.fastq'.format(extra_s1, extra_size, i))
    os.system('seqtk sample -s100 {} {} > subsets/isd_1/trial_6/{}/extra_R2.fastq'.format(extra_s2, extra_size, i))
    os.system('cat subsets/isd_1/trial_6/{0}/*_R1.fastq > subsets/isd_1/trial_6/{0}/{0}_R1.fastq'.format(i))
    os.system('cat subsets/isd_1/trial_6/{0}/*_R2.fastq > subsets/isd_1/trial_6/{0}/{0}_R2.fastq'.format(i))
