import os
import sys

def findreferenceNumber(name, database):
    if os.path.exists(database + ".name"):
        with open(database + ".name") as f:
            t = 1
            for line in f:
                if line.rstrip() == name:
                    return t
                else:
                    t += 1
        sys.exit('reference name not found in database name file.')
    else:
        print ("database name file not found. Please run kma index with -name option.")
        sys.exit(1)

rmlst_database = '/home/people/malhal/nanodecon_paper/tmp_tests/nanodecon_db/rmlst_db'

file_1 = 'NZ_CP051744/NZ_CP051744.res'
file_2 = 'NZ_CP091756/NZ_CP091756.res'

total_file = 'test_rmlst.res'

file_set_1 = set()
file_set_2 = set()
total_set = set()
with open(file_1, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            file_set_1.add(line[0])

with open(file_2, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            file_set_2.add(line[0])

with open(total_file, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            total_set.add(line[0])

print (file_1)
hit_count = 0
miss_count = 0
for item in file_set_1:
    if item in total_set:
        hit_count += 1
    else:
        miss_count += 1
        print (item)
print ('hit', hit_count)
print ('miss', miss_count)

print (file_2)
hit_unique_count = 0
hit_shared_count = 0
miss_count = 0
for item in file_set_2:
    if item in total_set:
        if item in file_set_1:
            hit_shared_count += 1
        else:
            hit_unique_count += 1
    else:
        miss_count += 1
        print (item)
print ('hit unique', hit_unique_count)
print ('hit shared', hit_shared_count)
print ('miss', miss_count)
sys.exit()

"""
count = 0
for item in file_set_1:
    if item in file_set_2:
        print (item)
    else:
        print ('not in file 2', item)
print (count)

        count += 1
"""

file_1_list = list()
file_2_list = list()

print ('In both:')
with open('test_CP091756.res', 'r') as f:
    for line in f:
        if not line.startswith('#'):
            if line.strip().split('\t')[0] in file_set_1:
                print (line.strip())

print ('Only in file 2:')
with open('test_CP091756.res', 'r') as f:
    for line in f:
        if not line.startswith('#'):
            if line.strip().split('\t')[0] not in file_set_1:
                print (line.strip())

"""
for item in file_set_1:
    print (item)
    file_1_list.append(str(findreferenceNumber(item, rmlst_database)))

for item in file_set_2:
    print (item)
    file_2_list.append(str(findreferenceNumber(item, rmlst_database)))



file_1_cmd = 'kma seq2fasta -seqs ' + ','.join(file_1_list) + ' -t_db ' + rmlst_database + ' > file_1.fasta'
file_2_cmd = 'kma seq2fasta -seqs ' + ','.join(file_2_list) + ' -t_db ' + rmlst_database + ' > file_2.fasta'

print (file_1_cmd)      
print (file_2_cmd)
os.system(file_1_cmd)
os.system(file_2_cmd)
"""