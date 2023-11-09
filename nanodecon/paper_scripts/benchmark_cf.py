import os
import sys
import numpy as np


path_1028222 = '/home/people/malhal/papers/rmlst/benchmarking/sample_benchmarks/1028222'

p1_1028222 = dict()
p3_1028222 = dict()
p5_1028222 = dict()


with open('{}/1p_1028222.txt'.format(path_1028222), 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        gene = line[0]
        mutations = line[1][1:-1]
        depths = line[2][1:-1]
        if mutations.split(',') != ['']:
            for i in range(len(mutations.split(','))):
                current_mutation = mutations.split(',')[i].strip().replace('\'', '')
                current_depth = depths.split(',')[i].strip().replace('\'', '')
                if gene in p1_1028222:
                    p1_1028222[gene][0].append(current_mutation)
                    p1_1028222[gene][1].append(current_depth)
                else:
                    p1_1028222[gene] = [[current_mutation], [current_depth]]
        else:
            p1_1028222[gene] = [[], []]

with open('{}/3p_1028222.txt'.format(path_1028222), 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        gene = line[0]
        mutations = line[1][1:-1]
        depths = line[2][1:-1]
        if mutations.split(',') != ['']:
            for i in range(len(mutations.split(','))):
                current_mutation = mutations.split(',')[i].strip().replace('\'', '')
                current_depth = depths.split(',')[i].strip().replace('\'', '')
                if gene in p3_1028222:
                    p3_1028222[gene][0].append(current_mutation)
                    p3_1028222[gene][1].append(current_depth)
                else:
                    p3_1028222[gene] = [[current_mutation], [current_depth]]
        else:
            p3_1028222[gene] = [[], []]

with open('{}/5p_1028222.txt'.format(path_1028222), 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        gene = line[0]
        mutations = line[1][1:-1]
        depths = line[2][1:-1]
        if mutations.split(',') != ['']:
            for i in range(len(mutations.split(','))):
                current_mutation = mutations.split(',')[i].strip().replace('\'', '')
                current_depth = depths.split(',')[i].strip().replace('\'', '')
                if gene in p5_1028222:
                    p5_1028222[gene][0].append(current_mutation)
                    p5_1028222[gene][1].append(current_depth)
                else:
                    p5_1028222[gene] = [[current_mutation], [current_depth]]
        else:
            p5_1028222[gene] = [[], []]

path_1028247 = '/home/people/malhal/papers/rmlst/benchmarking/sample_benchmarks/1028247'

p1_1028247 = dict()
p3_1028247 = dict()
p5_1028247 = dict()

with open('{}/1p_1028247.txt'.format(path_1028247), 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        gene = line[0]
        mutations = line[1][1:-1]
        depths = line[2][1:-1]
        if mutations.split(',') != ['']:
            for i in range(len(mutations.split(','))):
                current_mutation = mutations.split(',')[i].strip().replace('\'', '')
                current_depth = depths.split(',')[i].strip().replace('\'', '')
                if gene in p1_1028247:
                    p1_1028247[gene][0].append(current_mutation)
                    p1_1028247[gene][1].append(current_depth)
                else:
                    p1_1028247[gene] = [[current_mutation], [current_depth]]
        else:
            p1_1028247[gene] = [[], []]

with open('{}/3p_1028247.txt'.format(path_1028247), 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        gene = line[0]
        mutations = line[1][1:-1]
        depths = line[2][1:-1]
        if mutations.split(',') != ['']:
            for i in range(len(mutations.split(','))):
                current_mutation = mutations.split(',')[i].strip().replace('\'', '')
                current_depth = depths.split(',')[i].strip().replace('\'', '')
                if gene in p3_1028247:
                    p3_1028247[gene][0].append(current_mutation)
                    p3_1028247[gene][1].append(current_depth)
                else:
                    p3_1028247[gene] = [[current_mutation], [current_depth]]
        else:
            p3_1028247[gene] = [[], []]

with open('{}/5p_1028247.txt'.format(path_1028247), 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        gene = line[0]
        mutations = line[1][1:-1]
        depths = line[2][1:-1]
        if mutations.split(',') != ['']:
            for i in range(len(mutations.split(','))):
                current_mutation = mutations.split(',')[i].strip().replace('\'', '')
                current_depth = depths.split(',')[i].strip().replace('\'', '')
                if gene in p5_1028247:
                    p5_1028247[gene][0].append(current_mutation)
                    p5_1028247[gene][1].append(current_depth)
                else:
                    p5_1028247[gene] = [[current_mutation], [current_depth]]
        else:
            p5_1028247[gene] = [[], []]

final_mutation_file = '/home/people/malhal/papers/rmlst/benchmarking/sample_benchmarks/final_mutation_lists/biological_rmlst_snvs.txt'

biological_rmlst_snvs = dict()

with open(final_mutation_file, 'r') as f:
    for line in f:
        if line.startswith('B'):
            line = line.strip().split('>')
            gene = line[0].strip().split('_')[0]
            biological_rmlst_snvs[gene] = [line[0].strip(), line[1].strip(), []] #line[2] is 240_G_A format. Not 240_A
        else:
            line = line.strip()
            new_gene_format = line.split('_')[0] + '_' + line.split('_')[2]
            biological_rmlst_snvs[gene][2].append(new_gene_format)


def derive_biological_rmlst_snvs_results(biological_rmlst_snvs, confindr_dict):
    all_results_depths = np.array([])
    for gene in confindr_dict:
        depths = confindr_dict[gene][1]
        for depth in depths:
            all_results_depths = np.append(all_results_depths, float(depth))
    threshold = np.percentile(all_results_depths, 50) / 5 #Works fine for fast filtering

    total_snvs = 0
    total_bio_snvs_found_correctly = 0
    total_bio_snvs_found_incorrectly = 0
    found_bio_snvs = list()
    total_primary_snvs_found = 0
    remaining_snvs = dict()
    remaining_bio_snvs = dict()

    #Also consider non biological rMLST SNVs for primary and secondary samples
    #Consider removing low depth snps to give a better representation of resutls.
    for gene in confindr_dict:
        mutations = confindr_dict[gene][0]
        depths = confindr_dict[gene][1]
        remaining_snvs[gene] = [[], []]
        #print (gene, mutations)
        for mutation in mutations:
            if float(depths[mutations.index(mutation)]) > threshold:
                if gene == 'BACT000040' and mutation == '240_G': #In the test data BACT000040 is not complete in the primary samples
                    total_snvs += 1
                    total_bio_snvs_found_correctly += 1
                else:
                    if mutation in biological_rmlst_snvs[gene][2]:
                        total_snvs += 1
                        total_bio_snvs_found_correctly += 1
                        found_bio_snvs.append(mutation)
                    else:
                        if mutation in p5_1028222[gene][0]:
                            total_snvs += 1
                            total_primary_snvs_found += 1
                        else:
                            total_snvs += 1
                            total_bio_snvs_found_incorrectly += 1
                            remaining_snvs[gene][0].append(mutation)
                            remaining_snvs[gene][1].append(depths[mutations.index(mutation)])

    for item in biological_rmlst_snvs:
        for mutation in biological_rmlst_snvs[item][2]:
            if mutation not in found_bio_snvs:
                if mutation != '240_A':
                    remaining_bio_snvs[item] = mutation


    return total_snvs, total_bio_snvs_found_correctly, total_bio_snvs_found_incorrectly, remaining_snvs, remaining_bio_snvs, total_primary_snvs_found


rates = [1, 2, 3, 4, 5]
file_depths = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

print ('experiment,totalsnvs,totalbiologicalsnvscorrect,totalbiologicalsnvsincorrect,totalprimarysnvsfound')
for rate in rates:
    for file_depth in file_depths:
        confindr_file = '/home/people/malhal/papers/rmlst/benchmarking/confindr/{}/{}/intra_contamination.csv'.format(rate, file_depth)
        confindr_dict = dict()
        with open(confindr_file, 'r') as f:
            for line in f:
                if line.startswith('Gene'):
                    headers = line.strip().split(',')
                else:
                    line = line.strip().split(',')
                    allele = line[0]
                    gene = allele.split('_')[0]
                    mutation = line[1] + '_' + line[4].split(':')[0]
                    if gene not in confindr_dict:
                        confindr_dict[gene] = [[mutation], [line[7]]]
                    else:
                        confindr_dict[gene][0].append(mutation)
                        confindr_dict[gene][1].append(line[7])

        total_snvs, total_bio_snvs_found_correctly, total_bio_snvs_found_incorrectly, remaining_snvs, remaining_bio_snvs, total_primary_snvs_found = derive_biological_rmlst_snvs_results(biological_rmlst_snvs, confindr_dict)
        print (str(rate) + '/' + str(file_depth), total_snvs, total_bio_snvs_found_correctly, total_bio_snvs_found_incorrectly, total_primary_snvs_found, sep=',')
        #for item in remaining_snvs:
#    print (item, remaining_snvs[item])

#What about non biological mutations but true positives.
#print ('Remaining Biological SNVs')
#for item in remaining_bio_snvs:
#    print (item, remaining_bio_snvs[item])


#p1_1028247
#p3_1028247
#p5_1028247

#biological_rmlst_snvs