import os
import sys

blacklist_set = set()

with open('DTU_2023_1028247_100_SI_DADST035_SE_q20.fq', 'r') as f:
    for line in f:
        if line.startswith('@'):
            blacklist_set.add(line.split(' ')[0])

with open('DTU_2023_1028247_100_SI_DADST035_SE_q18.fq', 'r') as f:
    for line in f:
        if line.startswith('@'):
            blacklist_set.add(line.split(' ')[0])

with open('DTU_2023_1028247_100_SI_DADST035_SE_q16.fq', 'r') as f:
    for line in f:
        if line.startswith('@'):
            blacklist_set.add(line.split(' ')[0])

with open('DTU_2023_1028247_100_SI_DADST035_SE_q14.fq', 'r') as f:
    for line in f:
        if line.startswith('@'):
            blacklist_set.add(line.split(' ')[0])

with open('DTU_2023_1028247_100_SI_DADST035_SE_min_q_12.fq', 'w') as w:
    with open('DTU_2023_1028247_100_SI_DADST035_SE_q12.fq', 'r') as f:
        for line in f:
            if line.startswith('@'):
                if line.split(' ')[0] not in blacklist_set:
                    w.write(line)
                    w.write(f.readline())
                    w.write(f.readline())
                    w.write(f.readline())
