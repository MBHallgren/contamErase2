import os
import sys

def auto_identifiy_species(input, output, database):
    cmd = "kma -i {} -o {} -t_db {} -mem_mode -1t1 -t 8 -Sparse".format(input, output, database)
    os.system(cmd)
    with open(output + '.spa', 'r') as f:
        best_score = 0
        genome_size = 0
        for line in f:
            if not line.startswith('#'):
                score = float(line.split('\t')[2])
                if score > best_score:
                    best_score = score
                    best_line = line.split('\t')[0]
                    genome_size = int(line.split('\t')[4])

    return best_line, genome_size

path = '/home/people/malhal/data/test_intra/'
files = os.listdir(path)

for item in files:
    name = item.split('.')[0]
    best_line = auto_identifiy_species(path + item, name, '/home/people/malhal/databasenew/bacteria.ATG')
    print (name, best_line)




