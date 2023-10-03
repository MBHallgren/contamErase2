import os
import sys

res_dict_008253_5M = dict()
with open('008253_5M/008253_5M_R1.res', 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split()
            print (line)
            number = line[0].split('_')[0][-2:]
            print (number)
            depth = float(line[8])
            if depth > 10:
                if number not in res_dict_008253_5M:
                    res_dict_008253_5M[number] = [line[0].split('_')[-1]]
                else:
                    res_dict_008253_5M[number].append(line[0].split('_')[-1])

res_dict_026027 = dict()
with open('026027/026027_1M_R1.res', 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split()
            print (line)
            number = line[0].split('_')[0][-2:]
            print (number)
            depth = float(line[8])
            if depth > 3:
                if number not in res_dict_026027:
                    res_dict_026027[number] = [line[0].split('_')[-1]]
                else:
                    res_dict_026027[number].append(line[0].split('_')[-1])

res_dict_mix = dict()
with open('mix/mix_R1.res', 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split()
            print (line)
            number = line[0].split('_')[0][-2:]
            print (number)
            depth = float(line[8])
            if depth > 0.1:
                if number not in res_dict_mix:
                    res_dict_mix[number] = [line[0].split('_')[-1]]
                else:
                    res_dict_mix[number].append(line[0].split('_')[-1])

print (res_dict_008253_5M)
print (res_dict_026027)
print (res_dict_mix)

print ('008253_5M\t026027\tmix')

same = 0
same_but_one_typed = 0
same_both_typed = 0
for i in range(1, 66, 1):
    if i < 10:
        i = '0{}'.format(i)
    cmd = ''
    if str(i) in res_dict_008253_5M:
        cmd += str(res_dict_008253_5M[str(i)])
    else:
        cmd += 'None'
    cmd += '\t'
    if str(i) in res_dict_026027:
        cmd += str(res_dict_026027[str(i)])
    else:
        cmd += 'None'
    cmd += '\t'
    if str(i) in res_dict_mix:
        cmd += str(res_dict_mix[str(i)])
    else:
        cmd += 'None'
    split = cmd.split('\t')
    same_flag = False
    if split[0] == split[1] and 'None' not in split:
        cmd += '\tSAME'
        same += 1
        same_flag = True
    if 'None' not in cmd and not same_flag:
        print (cmd)

print ('SAME', same)