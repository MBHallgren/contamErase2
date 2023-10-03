import sys
import os

def compare_lists(lists, names):
    matrix = [[""] + names]  # Adding an empty cell in the top-left corner

    for i, list1 in enumerate(lists):
        row = [names[i]]  # First element in each row is the name
        for j, list2 in enumerate(lists):
            common_strings = len(set(list1) & set(list2))
            row.append(common_strings)
        matrix.append(row)

    return matrix

path = '/home/people/malhal/data/ecoli_results/'
file_list = os.listdir(path)
names = list()
lists = list()
for item in file_list:
    if not item.endswith('.py'):
        with open (path + item + '/' + item + '.res', 'r') as f:
            names.append(item)
            lists.append(list())
            for line in f:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    lists[-1].append(line[0])


matrix = compare_lists(lists, names)

for row in matrix:
    print(row)


"""
DTU_2023_1017265_R1 and DTU_2022_1013180_R1: 48/53, Different: BACT000060, BACT000064, BACT000049, BACT000036, BACT000006

DTU_2020_1000585_R1 and DTU_2022_1013180_R1: 50/53, Different: BACT000060, BACT000006, BACT000001

DTU_2020_1000585_R1 and DTU_2022_1013323_R1: 38/53 1, 2, 5, 11, 19, 30, 35, 40, 44, 48, 49, 53, 59, 62, 64

DTU_2022_1013323_R1 and DTU_2023_1018673_R1: 27/53 1, 2, 5, 9, 10, 11, 12, 20, 30, 32, 33, 35, 36, 38, 40, 45, 46, 48, 49, 50, 53, 57, 59, 60, 62, 64



"""