import os
import sys
import numpy as np


final_mutation_file = '/home/people/malhal/papers/rmlst/benchmarking/sample_benchmarks/final_mutation_lists/biological_rmlst_snvs.txt'

#final_mutation_file = '/home/people/malhal/papers/rmlst/benchmarking/sample_benchmarks/mix/mix_mutations.txt'

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
    total_bio_snvs_found_correctly = list()
    noise_novel = list()
    noise_known = list()

    #print (biological_rmlst_snvs)

    for gene in results:
        for mutation in results[gene][0]:
            total_snvs += 1
        if gene in biological_rmlst_snvs:
            mutations = results[gene][0]
            #print (mutations)
            for mutation in mutations:
                #if gene == 'BACT000040' and mutation == '240_G': #In the test data BACT000040 is not complete in the primary samples
                if mutation in biological_rmlst_snvs[gene]:
                    total_bio_snvs_found_correctly.append(gene + '_' + mutation)
                    found_bio_snvs.append(mutation)
                else:
                    if 'Novel' in results[gene][1][mutations.index(mutation)]: #This is a correction of a mistake, remove later
                        noise_novel.append(gene + '_' + mutation)
                    else:
                        noise_known.append(gene + '_' + mutation)

    #GO BACK AND CHECK THIS
    #for item in biological_rmlst_snvs:
    #    for mutation in biological_rmlst_snvs[item]:
    #        if mutation not in found_bio_snvs:
    #            if mutation != '240_A':
    #                remaining_bio_snvs[item] = mutation


    return total_snvs, total_bio_snvs_found_correctly, noise_novel, noise_known

rates = [1, 2, 3, 4, 5]
file_depths = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

print ('experiment,totalsnvs,totalbiologicalsnvscorrect,missingbiologicalsnvscorrect,noise_novel,noise_known')
#print ('MajorityAlelle,Position,MajorityBase,MutationBase,MutationDepth,TotalDepth')
for rate in rates:
    for file_depth in file_depths:
        results_file = '/home/people/malhal/papers/rmlst/benchmarking/saures_mine/{}/{}_results.csv'.format(rate, file_depth)
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
                            results[gene] = [[mutation], [line[-2]]] #mutation, depth, confindr
                        else:
                            results[gene][0].append(mutation)
                            results[gene][1].append(line[-2])

        total_snvs, total_bio_snvs_found_correctly, noise_novel, noise_known = derive_biological_rmlst_snvs_results(biological_rmlst_snvs, results)

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
