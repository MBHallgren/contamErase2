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
        allele = line[0]
        gene = line[0].split('_')[0]
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
        allele = line[0]
        gene = line[0].split('_')[0]
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
        allele = line[0]
        gene = line[0].split('_')[0]
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
        allele = line[0]
        gene = line[0].split('_')[0]
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
        allele = line[0]
        gene = line[0].split('_')[0]
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
        allele = line[0]
        gene = line[0].split('_')[0]
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
            allele = line[0]
            gene = line[1].split('_')[0]
            biological_rmlst_snvs[gene] = [line[0].strip(), line[1].strip(), []] #line[2] is 240_G_A format. Not 240_A
        else:
            line = line.strip()
            new_allele_format = line.split('_')[0] + '_' + line.split('_')[2]
            biological_rmlst_snvs[gene][2].append(new_allele_format)


def derive_biological_rmlst_snvs_results(biological_rmlst_snvs, confindr_dict):
    all_results_depths = np.array([])
    for gene in confindr_dict:
        depths = confindr_dict[gene][1]
        for depth in depths:
            all_results_depths = np.append(all_results_depths, float(depth))
    threshold = np.percentile(all_results_depths, 50) / 5 #Works fine for fast filtering

    total_snvs = 0
    total_bio_snvs_found_correctly = list()
    missing_bio_snvs = list()
    found_bio_snvs = list()
    total_primary_snvs_found = list()
    remaining_snvs = dict()
    remaining_bio_snvs = dict()
    noise = list()

    #Also consider non biological rMLST SNVs for primary and secondary samples
    #Consider removing low depth snps to give a better representation of resutls.
    for gene in confindr_dict:
        mutations = confindr_dict[gene][0]
        depths = confindr_dict[gene][1]
        remaining_snvs[gene] = [[], []]
        #print (allele, mutations)
        for mutation in mutations:
            if float(depths[mutations.index(mutation)]) > threshold:
                if gene == 'BACT000040' and mutation == '240_G': #In the test data BACT000040 is not complete in the primary samples
                    total_snvs += 1
                    total_bio_snvs_found_correctly.append(gene + '_' + mutation)
                else:
                    if mutation in biological_rmlst_snvs[gene][2]:
                        total_snvs += 1
                        total_bio_snvs_found_correctly.append(gene + '_' + mutation)
                        found_bio_snvs.append(mutation)
                    else:
                        if mutation in p5_1028222[gene][0]:
                            total_snvs += 1
                            total_primary_snvs_found.append(gene + '_' + mutation)
                        else:
                            total_snvs += 1
                            noise.append(gene + '_' + mutation)
                            remaining_snvs[gene][0].append(mutation)
                            remaining_snvs[gene][1].append(depths[mutations.index(mutation)])

    for item in biological_rmlst_snvs:
        for mutation in biological_rmlst_snvs[item][2]:
            if mutation not in found_bio_snvs:
                if mutation != '240_A':
                    remaining_bio_snvs[item] = mutation

    return total_snvs, total_bio_snvs_found_correctly, total_primary_snvs_found, noise, remaining_snvs, remaining_bio_snvs

rates = [1, 2, 3, 4, 5]
file_depths = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

print ('experiment,totalsnvs,totalbiologicalsnvscorrect,missingbiologicalsnvscorrect,totalprimarysnvsfound, noise')
for rate in rates:
    for file_depth in file_depths:
        confindr_file = '/home/people/malhal/papers/rmlst/benchmarking/confindr/{}/{}/intra_contamination.csv'.format(rate, file_depth)
        confindr_dict = dict()
        with open(confindr_file, 'r') as f:
            for line in f:
                if line.startswith('allele'):
                    headers = line.strip().split(',')
                else:
                    line = line.strip().split(',')
                    allele = line[0]
                    gene = allele.split('_')[0]
                    mutation = line[1] + '_' + line[4].split(':')[0]
                    if gene.startswith('BACT'):
                        if gene not in confindr_dict:
                            confindr_dict[gene] = [[mutation], [line[7]]]
                        else:
                            confindr_dict[gene][0].append(mutation)
                            confindr_dict[gene][1].append(line[7])


        total_snvs, total_bio_snvs_found_correctly, total_primary_snvs_found, noise, remaining_snvs, remaining_bio_snvs = derive_biological_rmlst_snvs_results(biological_rmlst_snvs, confindr_dict)

        all_biological_snvs = list()
        missing_biological_snvs = list()
        for gene in biological_rmlst_snvs:
            for mutation in biological_rmlst_snvs[gene][2]:
                all_biological_snvs.append(gene + '_' + mutation)

        for item in all_biological_snvs:
            if item not in total_bio_snvs_found_correctly:
                missing_biological_snvs.append(item)
        print (missing_biological_snvs)
        #BACT000038 is a problem. 4 not IDd. Are these called correctly? check again. TBD
        #BACT000030 2 mutations. ID'ed with co-occurence!
        #BACT000040_240_A #This is known and needs to be fixed.' #Not perfect rmlst alignment
        #['BACT000030_489_T', 'BACT000030_492_T', 'BACT000038_222_T', 'BACT000038_225_C', 'BACT000038_369_A', 'BACT000038_370_T', 'BACT000040_240_A']

        #30 doesn't find any of the 2. Check multiple alignment again of consensus sequences.
        #38 doesn't find the 4. Check multiple alignment again of consensus sequences.
        #40, 240 A is found. Look into why this is an error. Are we expecting something else?

        print (str(rate) + '/' + str(file_depth), total_snvs, len(total_bio_snvs_found_correctly), len(missing_biological_snvs), len(total_primary_snvs_found), len(noise), sep=',')

        #TBD current we use 5p mutations for primary SNVs, should also bench are 1-4?
        #TBD Which missing biological SNVs do we see and why?
