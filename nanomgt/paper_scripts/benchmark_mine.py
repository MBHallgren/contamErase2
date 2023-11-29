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

#final_mutation_file = '/home/people/malhal/papers/rmlst/benchmarking/sample_benchmarks/final_mutation_lists/biological_rmlst_snvs.txt'

final_mutation_file = '/home/people/malhal/papers/rmlst/benchmarking/sample_benchmarks/mix/mix_mutations.txt'

biological_rmlst_snvs = dict()

with open(final_mutation_file, 'r') as f:
    for line in f:
        if line.startswith('B'):
            gene = line.strip()
            biological_rmlst_snvs[gene] = [] #line[2] is 240_G_A format. Not 240_A
        else:
            line = line.strip()
            mutation = line.split('_')[0] + '_' + line.split('_')[2]
            biological_rmlst_snvs[gene].append(mutation)


def derive_biological_rmlst_snvs_results(biological_rmlst_snvs, results):

    #sys.exit()
    #Rework how we derive biological_rmlst_snvs

    #all_results_depths = np.array([])
    #for gene in results:
    ##    depths = results[gene][1]
    #    for depth in depths:
    #        all_results_depths = np.append(all_results_depths, float(depth))
    #threshold = np.percentile(all_results_depths, 50) / 5 #Works fine for fast filtering

    total_snvs = 0
    found_bio_snvs = list()
    remaining_snvs = dict()
    total_bio_snvs_found_correctly = list()
    remaining_bio_snvs = dict()
    noise_novel = list()
    noise_known = list()

    for gene in results:
        if gene in biological_rmlst_snvs:
            mutations = results[gene][0]
            remaining_snvs[gene] = [[], []]
            for mutation in mutations:
                #if gene == 'BACT000040' and mutation == '240_G': #In the test data BACT000040 is not complete in the primary samples
                    #total_snvs += 1
                    #total_bio_snvs_found_correctly.append(gene + '_' + mutation)
                #    pass
                #else:
                    #print (gene, mutation)
                if mutation in biological_rmlst_snvs[gene]:
                    total_snvs += 1
                    total_bio_snvs_found_correctly.append(gene + '_' + mutation)
                    found_bio_snvs.append(mutation)
                else:
                    if 'never seen' in results[gene][2][mutations.index(mutation)]: #This is a correction of a mistake, remove later
                        total_snvs += 1
                        noise_novel.append(gene + '_' + mutation)
                        remaining_snvs[gene][0].append(mutation)
                    else:
                        total_snvs += 1
                        noise_known.append(gene + '_' + mutation)
                        remaining_snvs[gene][0].append(mutation)

    #GO BACK AND CHECK THIS
    for item in biological_rmlst_snvs:
        for mutation in biological_rmlst_snvs[item]:
            if mutation not in found_bio_snvs:
                if mutation != '240_A':
                    remaining_bio_snvs[item] = mutation


    return total_snvs, total_bio_snvs_found_correctly, remaining_snvs, remaining_bio_snvs, noise_novel, noise_known

rates = [1, 2, 3, 4, 5]
file_depths = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

print ('experiment,totalsnvs,totalbiologicalsnvscorrect,missingbiologicalsnvscorrect,noise_novel,noise_known')
#print ('MajorityAlelle,Position,MajorityBase,MutationBase,MutationDepth,TotalDepth')
for rate in rates:
    for file_depth in file_depths:
        results_file = '/home/people/malhal/papers/rmlst/benchmarking/mix/{}/{}_results.csv'.format(rate, file_depth)
        #confindr_file = '/home/people/malhal/papers/rmlst/benchmarking/confindr/{}/{}/intra_contamination.csv'.format(rate, file_depth)
        results = dict()
        with open(results_file, 'r') as f:
            for line in f:
                if line.startswith('Gene'):
                    headers = line.strip().split(',')
                    headers = line[1:]
                else:
                    line=line.replace('. ', ',')
                    line = line.strip().split(',')
                    allele = line[0]
                    gene = allele.split('_')[0]
                    mutation = line[1] + '_' + line[3]
                    if gene.startswith('BACT'):
                        if gene not in results:
                            results[gene] = [[mutation], [line[4]], [line[-1]]]
                        else:
                            results[gene][0].append(mutation)
                            results[gene][1].append(line[4])
                            results[gene][2].append(line[-1])

        total_snvs, total_bio_snvs_found_correctly, remaining_snvs, remaining_bio_snvs, noise_novel, noise_known = derive_biological_rmlst_snvs_results(biological_rmlst_snvs, results)

        all_biological_snvs = list()
        missing_biological_snvs = list()
        for gene in biological_rmlst_snvs:
            for mutation in biological_rmlst_snvs[gene]:
                all_biological_snvs.append(gene + '_' + mutation)

        for item in all_biological_snvs:
            if item not in total_bio_snvs_found_correctly:
                missing_biological_snvs.append(item)
                #print (item)
        #print (missing_biological_snvs)
        #BACT000038 is a problem. 4 not IDd. Are these called correctly? check again. TBD
        #BACT000030 2 mutations. ID'ed with co-occurence!
        #BACT000040_240_A #This is known and needs to be fixed.'    240_A' is not in list
        print (str(rate) + '/' + str(file_depth), total_snvs, len(total_bio_snvs_found_correctly), len(missing_biological_snvs), len(noise_novel), len(noise_known), sep=',')
        #TBD current we use 5p mutations for primary SNVs, should also bench are 1-4?
        #TBD Which missing biological SNVs do we see and why?