import os
import sys

data_path = '/home/projects/cge/backup/fastq/food_seqdata/'

all_data = set()

output_ids = set()
with open('/home/people/malhal/ecoli.txt') as f:
    for line in f:
        line = line.strip().split(' ')
        output_ids.add(line[0])


data_size_dict = {}
t = 0
with open ('/home/projects/cge/people/malhal/malte_request.tsv', 'r') as f:
    line = f.readline()
    header = line.strip().split('\t')
    for line in f:
        t +=1
        if t % 100 == 0:
            print (t)
        line = line.strip().split('\t')
        data_size_dict[line[0]] = [0, '']
        if line[2] == '1':
            file_path = '{}/{}/Raw_data'.format(data_path, line[-1])
            files = os.listdir(file_path)
            for item in files:
                if line[0] in item:
                    data_size_dict[line[0]][0] += os.path.getsize('{}/{}'.format(file_path, item))
                    data_size_dict[line[0]][1] = line
                    if line[0] in output_ids:
                        #print (line[0], line[-1], item)
                        os.system('cp {} /home/projects/cge/people/malhal/nanodecon_data/{}'.format('{}/{}'.format(file_path, item), item))
sys.exit()
for item in data_size_dict:
    if int(data_size_dict[item][0]/1000000) > 500:
        print (item, int(data_size_dict[item][0]/1000000), data_size_dict[item][1])