import os
import sys

def build_consensus_dict(res_file, mat_file):
    gene_score_dict, top_alleles, allele_pair_dict = extract_max_scored_alleles(res_file)
    non_alignment_matches = {}
    consensus_dict = {}
    odd_size_alleles = set()
    correct_size_alleles = set()
    with open(res_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                allele = line[0]
                gene = allele.split('_')[0]
                if int(line[3]) != gene_score_dict[gene][-1][0]:
                    odd_size_alleles.add(line[0])
                else:
                    correct_size_alleles.add(allele)
                    if gene not in consensus_dict:
                        consensus_dict[gene] = []
                        for i in range(int(line[3])):
                            consensus_dict[gene].append([0, 0, 0, 0, 0, 0]) #[A, C, G, T, N, -]

    with open(mat_file, 'r') as f:
        correct_size_flag = False
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                allele = line[1:]
                gene = allele.split('_')[0]
                if allele in correct_size_alleles:
                    correct_size_flag = True
                    index = 0
                else:
                    correct_size_flag = False
            else:
                if line != '':
                    if correct_size_flag:
                        line = line.split('\t')
                        if line[0] != '-': #excludes read gaps. Reconsider?
                            line = line[1:]
                            for i in range(len(line)):
                                consensus_dict[gene][index][i] += int(line[i])
                            index += 1
    return odd_size_alleles, non_alignment_matches, consensus_dict, top_alleles, allele_pair_dict, gene_score_dict

def extract_max_scored_alleles(res_file):
    gene_score_dict = {}
    max_scored_alleles = set()
    allele_pair_dict = {}

    with open(res_file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                line_split = line.split("\t")
                allele = line_split[0]
                gene = allele.split("_")[0]
                allele_pair_dict[gene] = ''
                mapscore = int(line_split[1])
                value = [int(line_split[3]), int(line_split[1]), line_split[0]]

                # If gene is not present in dict or has lower score, update/insert
                if gene not in gene_score_dict or mapscore > gene_score_dict[gene][1]:
                    gene_score_dict[gene] = (allele, mapscore, value)

    # Create dictionary with allele as key and value as specified
    for gene, (allele, _, value) in gene_score_dict.items():
        max_scored_alleles.add(allele)
        allele_pair_dict[gene] = allele



    return gene_score_dict, max_scored_alleles, allele_pair_dict


def derive_mutation_positions2(consensus_dict, min_n, urd, lrd):
    upper_confirmed_mutation_dict = {}
    lower_confirmed_mutation_dict = {}
    all_confirmed_mutation_dict = {}
    for item in consensus_dict:
        upper_confirmed_mutation_dict[item] = [[], []]
        lower_confirmed_mutation_dict[item] = [[], []]
        all_confirmed_mutation_dict[item] = [[], []]
        for i in range(len(consensus_dict[item])):
            positions = consensus_dict[item][i][:4]
            max_number = max(positions)
            index_of_max = positions.index(max_number)
            nucleotide_index = ['A', 'C', 'G', 'T']
            for t in range(len(positions)):
                if t != index_of_max: #Use this to parse majority calls in the future if needed
                    if positions[t] >= min_n:
                        total_depth = sum(positions)
                        relative_depth = positions[t] / total_depth
                        if relative_depth >= urd:
                            upper_confirmed_mutation_dict[item][0].append(
                                '{}_{}'.format(i + 1, nucleotide_index[t]))
                            upper_confirmed_mutation_dict[item][1].append(positions[t])
                            all_confirmed_mutation_dict[item][0].append(
                                '{}_{}'.format(i + 1, nucleotide_index[t]))
                            all_confirmed_mutation_dict[item][1].append(positions[t])

                        elif relative_depth >= lrd:
                            lower_confirmed_mutation_dict[item][0].append(
                                '{}_{}'.format(i + 1, nucleotide_index[t]))
                            lower_confirmed_mutation_dict[item][1].append(positions[t])
                            all_confirmed_mutation_dict[item][0].append(
                                '{}_{}'.format(i + 1, nucleotide_index[t]))
                            all_confirmed_mutation_dict[item][1].append(positions[t])
    return upper_confirmed_mutation_dict, lower_confirmed_mutation_dict, all_confirmed_mutation_dict


odd_size_alleles, non_alignment_matches, consensus_dict, top_alleles, allele_pair_dict, gene_score_dict = build_consensus_dict(
'1028247.res',
'1028247.mat')

upper_confirmed_mutation_dict, lower_confirmed_mutation_dict, all_confirmed_mutation_dict = derive_mutation_positions2(
    consensus_dict, 2, 0.05, 0.01)

for item in all_confirmed_mutation_dict:
    print(item, all_confirmed_mutation_dict[item][0], all_confirmed_mutation_dict[item][1], sep='\t')