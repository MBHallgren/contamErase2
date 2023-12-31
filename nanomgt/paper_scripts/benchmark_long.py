import os
import sys
import numpy as np

#final_mutation_file = '/home/people/malhal/papers/rmlst/benchmarking/sample_benchmarks/final_mutation_lists/biological_rmlst_snvs.txt'

final_mutation_file = '/home/people/malhal/papers/rmlst/benchmarking/sample_benchmarks/mix/mix_mutations.txt'

biological_rmlst_snvs = dict()

with open(final_mutation_file, 'r') as f:
    for line in f:
        if line.startswith('B'):
            if '>' in line:
                gene = line.strip().split('>')[0].split('_')[0]
            else:
                gene = line.strip()
            biological_rmlst_snvs[gene] = [] #line[2] is 240_G_A format. Not 240_A
        else:
            line = line.strip()
            mutation = line.split('_')[0] + '_' + line.split('_')[2]
            biological_rmlst_snvs[gene].append(mutation)


def derive_biological_rmlst_snvs_results(biological_rmlst_snvs, results):

    total_snvs = 0
    found_bio_snvs = list()
    remaining_snvs = dict()
    total_bio_snvs_found_correctly = list()
    remaining_bio_snvs = dict()
    noise_novel = list()
    noise_known = list()

    for gene in results:
        if gene in biological_rmlst_snvs:
            mutations = results[gene]
            remaining_snvs[gene] = []
            for mutation in mutations:
                total_snvs += 1
                if mutation in biological_rmlst_snvs[gene]:
                    total_bio_snvs_found_correctly.append(gene + '_' + mutation)
                    found_bio_snvs.append(mutation)
                else:
                    noise_known.append(gene + '_' + mutation)
                    remaining_snvs[gene].append(mutation)


    return total_snvs, total_bio_snvs_found_correctly, remaining_snvs, remaining_bio_snvs, noise_novel, noise_known

rates = [1, 2, 3, 4, 5]
file_depths = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

#print (biological_rmlst_snvs)

print ('experiment,totalsnvs,totalbiologicalsnvscorrect,missingbiologicalsnvscorrect,noise_novel,noise_known')
#print ('MajorityAlelle,Position,MajorityBase,MutationBase,MutationDepth,TotalDepth')
for rate in rates:
    for file_depth in file_depths:
        #results_file = '/home/people/malhal/papers/rmlst/benchmarking/mix/{}/{}_results.csv'.format(rate, file_depth)
        #results_file = '/home/people/malhal/papers/rmlst/benchmarking/confindr/{}/{}/intra_contamination.csv'.format(rate, file_depth)
        #confindr_file = '/home/people/malhal/papers/rmlst/benchmarking/confindr/{}/{}/intra_contamination.csv'.format(rate, file_depth)
        results_file = '/home/people/malhal/papers/rmlst/benchmarking/ecoli_long/{}/{}/output.vcf'.format(rate, file_depth)
        #print (results_file)
        results = dict()
        with open(results_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    gene = line[0].split('_')[0]
                    if gene not in results:
                        results[gene] = []
                    mutation = line[1] + '_' + line[4]
                    results[gene].append(mutation)

        total_snvs, total_bio_snvs_found_correctly, remaining_snvs, remaining_bio_snvs, noise_novel, noise_known = derive_biological_rmlst_snvs_results(biological_rmlst_snvs, results)

        all_biological_snvs = list()
        missing_biological_snvs = list()
        for gene in biological_rmlst_snvs:
            for mutation in biological_rmlst_snvs[gene]:
                all_biological_snvs.append(gene + '_' + mutation)

        for item in all_biological_snvs:
            if item not in total_bio_snvs_found_correctly:
                missing_biological_snvs.append(item)

        print (str(rate) + '/' + str(file_depth), total_snvs, len(total_bio_snvs_found_correctly), len(missing_biological_snvs), len(noise_novel), len(noise_known), sep=',')
        #TBD current we use 5p mutations for primary SNVs, should also bench are 1-4?
        #TBD Which missing biological SNVs do we see and why?
