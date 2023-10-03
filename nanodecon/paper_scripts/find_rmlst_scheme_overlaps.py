import os
import sys

file_1 = '1028222.res'
file_2 = '1028247.res'

list_1 = list()
list_2 = list()

with open (file_1, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            if not int(line[1]) < 1000 and not float(line[-3]) < 4:
                list_1.append(line[0])

with open (file_2, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            if not int(line[1]) < 1000 and not float(line[-3]) < 4:
                list_2.append(line[0])

rmlst_file = '/home/people/malhal/contamErase/test/rmlst_scheme.txt'

rmlst_scores = dict()

with open(rmlst_file, 'r') as f:
    for line in f:
        if not line.startswith('rST'):
            line = line.strip().split('\t')
            rmlst_genes = line[1:54]
            rmlst_scores[line[0]] = set()
            for i in range(len(rmlst_genes)):
                full_name = headers[i] + '_' + rmlst_genes[i]
                rmlst_scores[line[0]].add(full_name)
        else:
            line = line.strip().split('\t')
            headers = line[1:54]

print ('scoring dict complete')
top_1 = 0
top_2 = 0
top_set_1 = set()
top_set_2 = set()
for item in rmlst_scores:
    hits_1 = 0
    hits_2 = 0
    for gene in rmlst_scores[item]:
        if gene in list_1:
            hits_1 += 1
        if gene in list_2:
            hits_2 += 1
    if hits_1 >= 46:
        print (item, hits_1, '1')
        os.system('grep \'{} \' /home/people/malhal/contamErase/test/rmlst_scheme.txt'.format(item))
    if hits_2 >= 46:
        print (item, hits_2, '2')
        os.system('grep \'{} \' /home/people/malhal/contamErase/test/rmlst_scheme.txt'.format(item))
    if hits_1 > top_1:
        top_1 = hits_1
        top_set_1 = rmlst_scores[item]
    if hits_2 > top_2:
        top_2 = hits_2
        top_set_2 = rmlst_scores[item]

print (top_1, top_set_1)
print (top_2, top_set_2)
print ('')

non_shared_genes = top_set_1.symmetric_difference(top_set_2)

unique_genes = set()

for item in non_shared_genes:
    unique_genes.add(item.split('_')[0])

print (len(unique_genes), unique_genes)