import os
import sys
import logging
import subprocess
from nanodecon import kmergenetyper
from nanodecon import util
from nanodecon import kma
import yaml

def quick_determine_multiple_rmlst_alleles_from_res_files(res_file, cutoff):
    rmlst_dict = {}
    with open(res_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                gene = line.split()[0].split('_')[0]
                if gene not in rmlst_dict:
                    depth = float(line.split('\t')[-3])
                    if depth >= cutoff:
                        rmlst_dict[gene] = {}
                        allele = line.split()[0]
                        rmlst_dict[gene][allele] = depth
                else:
                    depth = float(line.split('\t')[-3])
                    if depth >= cutoff:
                        allele = line.split()[0]
                        rmlst_dict[gene][allele] = depth
    return rmlst_dict


def determine_intra_species_contamination_nanopore(arguments, rmlst_candidates, primary, candidate_dict, total_bases, ugt):
    candidate_rmlst_dict = {}
    primary_species = primary.split()[1].lower() + '_' + primary.split()[2].lower()
    os.system('mkdir ' + arguments.output + '/rmlst_' + primary_species)
    os.system('mkdir ' + arguments.output + '/rmlst_' + primary_species + '/fastas')
    for item in rmlst_candidates:
        candidate_rmlst_dict = util.derive_genes_ids_from_rmlst_candidate(item, arguments, primary_species, candidate_rmlst_dict)
    intra_species_reads = '{}/{}_read_pool.fastq'.format(arguments.output, primary_species)
    base_count = util.number_of_bases_in_file(intra_species_reads)
    relative_minimum_depth = (base_count / candidate_dict[primary][1]) * 0.01
    print('relative_minimum_depth', relative_minimum_depth)
    if relative_minimum_depth < 1.5:
        print('Relative minimum depth is less than 1.5, setting to 1.5')
        relative_minimum_depth = 1.5

    kma.KMARunner('{}'.format(intra_species_reads),
                  arguments.output + '/rmlst_' + primary_species + '/rmlst_' + primary_species,
                  arguments.db_dir + '/rmlst_db',
                  "-t {} -ID 10 -ont -md {} -matrix".format(arguments.threads, relative_minimum_depth))\
        .run()

    os.system('gunzip ' + arguments.output + '/rmlst_' + primary_species + '/rmlst_' + primary_species + '.mat.gz')
    candidate_rmlst_dict_results = build_rmlst_consensus_seuqences(arguments,
                                    arguments.output + '/rmlst_' + primary_species + '/rmlst_' + primary_species + '.res',
                                    arguments.output + '/rmlst_' + primary_species + '/rmlst_' + primary_species + '.mat',
                                    candidate_rmlst_dict,
                                    primary)

    return candidate_rmlst_dict_results

def build_rmlst_consensus_seuqences(arguments, rmlst_res_file, rmlst_mat_file, candidate_rmlst_dict, primary):
    length_dict = {}
    consensus_dict = {}
    #TBD CONSIDER GAPS AND NON LENGTH HITS
    with open(rmlst_res_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                gene = line.split()[0].split('_')[0]
                length = line.split()[3]
                comb = gene + '_' + length
                if comb not in consensus_dict:
                    consensus_dict[comb] = list()
                    for i in range(int(length)):
                        consensus_dict[comb].append([0, 0, 0, 0, 0, 0]) #A, C, G, T, N, -
                length_dict[line.split()[0]] = length

    #rework
    with open(rmlst_mat_file, 'r') as f:
        for line in f:
            if not line.strip() == '' and not line.startswith('-'):
                if line.startswith('#'):
                    index = 0
                    header = line.strip()[1:]
                    current_gene_length = length_dict[header]
                    comp = header.split('_')[0] + '_' + current_gene_length
                else:
                    ad6 = line.strip().split('\t')[1:]
                    for i in range(len(ad6)):
                        consensus_dict[comp][index][i] += int(ad6[i])
                    index += 1

    adjusted_error_rate_per_base = calculate_random_error_rate(consensus_dict)
    print ('adjusted_error_rate_per_base', adjusted_error_rate_per_base)
    primary_species = primary.split()[1].lower() + '_' + primary.split()[2].lower()

    candidate_rmlst_dict_results = dict()
    with open(arguments.output + '/rmlst_' + primary_species + '/rmlst_results.txt', 'w') as f:
        for item in candidate_rmlst_dict:
            candidate_rmlst_dict_results[item] = dict()
            if item != primary:
                print ('#' + item, file=f)
                validated_count = 0
                validated_unique_gene_count = 0
                unvalidated_count = 0
                validated_genes = set()
                #validation step doesnt work. Prioritze template_cov > depth
                #Consider if genes is not 100 complete in template_cov
                    #Is it the same as primary, if so is it incomplete?
                for gene in candidate_rmlst_dict[item]:
                    sequence = derive_sequence_from_template(gene, arguments.db_dir + '/rmlst_db')
                    result, score, zeros = validate_sequence(arguments, sequence, consensus_dict, adjusted_error_rate_per_base, gene, primary)
                    if gene.split('_')[0] not in validated_genes: #Don't varify multiple allels of the same gene in case of multiple IDs from database error
                        if result:
                            validated_genes.add(gene.split('_')[0])
                            validated_count += 1
                            if gene not in candidate_rmlst_dict[primary]:
                                validated_unique_gene_count += 1
                            print (gene, score, file=f)
                        else:
                            print (gene, score, file=f)
                            unvalidated_count += 1

                if arguments.conditions is None:
                    conditions = {
                        'scenario_1': {
                            'total_varified_genes': 53,
                            'total_unique_genes': 1
                        },
                        'scenario_2': {
                            'total_varified_genes': 52,
                            'total_unique_genes': 3
                        },
                        'scenario_3': {
                            'total_varified_genes': 50,
                            'total_unique_genes': 5
                        },
                        'scenario_4': {
                            'total_varified_genes': 40,
                            'total_unique_genes': 15
                        }

                    }
                else:
                    with open(arguments.conditions, 'r') as yaml_file:
                        conditions = yaml.safe_load(yaml_file)

                for condition in conditions:
                    if validated_count > conditions[condition]['total_varified_genes'] and validated_unique_gene_count > conditions[condition]['total_unique_genes']:
                        candidate_rmlst_dict_results[item]['rmlst_id'] = True
                        candidate_rmlst_dict_results[item]['validated_count'] = validated_count
                        candidate_rmlst_dict_results[item]['validated_unique_gene_count'] = validated_unique_gene_count
                        candidate_rmlst_dict_results[item]['unvalidated_count'] = unvalidated_count
                        break

                if 'rmlst_id' not in candidate_rmlst_dict_results[item]:
                    candidate_rmlst_dict_results[item]['rmlst_id'] = False
                    candidate_rmlst_dict_results[item]['validated_count'] = validated_count
                    candidate_rmlst_dict_results[item]['validated_unique_gene_count'] = validated_unique_gene_count
                    candidate_rmlst_dict_results[item]['unvalidated_count'] = unvalidated_count
            else:
                candidate_rmlst_dict_results[item]['rmlst_id'] = True
            if item != primary:
                print ('#Validated genes: {}'.format(candidate_rmlst_dict_results[item]['validated_count']), file=f)
                print ('#Validated unique genes: {}'.format(candidate_rmlst_dict_results[item]['validated_unique_gene_count']), file=f)
                print ('#Unvalidated genes: {}'.format(candidate_rmlst_dict_results[item]['unvalidated_count']), file=f)
        return candidate_rmlst_dict_results

def validate_sequence(arguments, sequence, consensus_dict, adjusted_error_rate_per_base, gene, primary):
    primary_species = primary.split()[1].lower() + '_' + primary.split()[2].lower()
    comb = gene.split('_')[0] + '_' + str(len(sequence))
    threshold = 1
    hits = 0
    misses = 0
    zeros = 0
    depth_list = []
    hit_list = []
    if comb in consensus_dict:
        for index in range(len(sequence)):
            new_count = [max(val - adjusted_error_rate_per_base, 0) for val in consensus_dict[comb][index][0:4]] #consider removing adjusted_error_rate_per_base
            total_count = sum(new_count)
            if total_count == 0:
                zeros += 1
                depth_list.append(consensus_dict[comb][index])
                hit_list.append('Zero')
            else:
                base_call = sequence[index]
                if base_call == 'A':
                    base_index = 0
                elif base_call == 'C':
                    base_index = 1
                elif base_call == 'G':
                    base_index = 2
                elif base_call == 'T':
                    base_index = 3

                depth_score = (new_count[base_index] / total_count) #Modify this to scale with depth
                if depth_score > arguments.isd:
                    hits += 1
                    depth_list.append(consensus_dict[comb][index][0:4])
                    hit_list.append('|')
                else:
                    misses += 1
                    depth_list.append(consensus_dict[comb][index][0:4])
                    hit_list.append('X')
        score = hits / (hits + misses + zeros)
    else:
        score = 0
        zeros = 0

    with open(arguments.output + '/rmlst_' + primary_species + '/rmlst_genes_depths.txt', 'a') as f:
        print (gene, depth_list, hit_list, file=f)

    if score >= threshold:
        return True, score, zeros
    else:
        return False, score, zeros




def derive_sequence_from_template(template, database):
    number = util.findTemplateNumber(template, database)
    cmd = 'kma seq2fasta -seqs {} -t_db {}'.format(number, database)
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode == 0:
        return result.stdout.strip().split('\n')[1]
    else:
        print("Error:", result.stderr)
        return None

def derive_gene_dict(arguments, rmlst_res_file, rmlst_mat_file):
    pass


def get_two_lowest(numbers):
    sorted_numbers = sorted(numbers)
    return sorted_numbers[:2]


def calculate_random_error_rate(consensus_dict):
    base_count = 0
    error_count = 0
    for item in consensus_dict:
        for index in range(len(consensus_dict[item])):
            base_count += sum(consensus_dict[item][index][0:4])
            error_count += sum(get_two_lowest(consensus_dict[item][index][0:4]))
    return (error_count/2) / base_count

def validate_rmlst_genes_illumina():
    #Only include alelle of same length as primary's alelle
    #Do direct validation from .res file, if depth > threshold
    pass

def determine_intra_species_contamination_illumina(arguments, rmlst_candidates, primary, candidate_dict, total_bases, ugt):
    candidate_rmlst_dict = {}
    primary_species = primary.split()[1].lower() + '_' + primary.split()[2].lower()
    os.system('mkdir ' + arguments.output + '/rmlst_' + primary_species)
    os.system('mkdir ' + arguments.output + '/rmlst_' + primary_species + '/fastas')
    for item in rmlst_candidates:
        candidate_rmlst_dict = util.derive_genes_ids_from_rmlst_candidate(item, arguments, primary_species, candidate_rmlst_dict)
    intra_species_reads_1 = '{}/{}_read_pool_R1.fastq'.format(arguments.output, primary_species)
    intra_species_reads_2 = '{}/{}_read_pool_R2.fastq'.format(arguments.output, primary_species)
    base_count = util.number_of_bases_in_file(intra_species_reads_1) + util.number_of_bases_in_file(intra_species_reads_2)
    relative_minimum_depth = (base_count / candidate_dict[primary][1]) * 0.01
    print('relative_minimum_depth', relative_minimum_depth)
    if relative_minimum_depth < 1.5:
        print('Relative minimum depth is less than 1.5, setting to 1.5')
        relative_minimum_depth = 1.5

    kma.KMARunner('{} {}'.format(intra_species_reads_1, intra_species_reads_2),
                  arguments.output + '/rmlst_' + primary_species + '/rmlst_' + primary_species,
                  arguments.db_dir + '/rmlst_db',
                  "-1t1 -t {} -ID 10 -md {} -matrix -vcf -oa".format(arguments.threads, relative_minimum_depth))\
        .run()

    os.system('gunzip ' + arguments.output + '/rmlst_' + primary_species + '/rmlst_' + primary_species + '.mat.gz')
    candidate_rmlst_dict_results = build_rmlst_consensus_seuqences(arguments,
                                    arguments.output + '/rmlst_' + primary_species + '/rmlst_' + primary_species + '.res',
                                    arguments.output + '/rmlst_' + primary_species + '/rmlst_' + primary_species + '.mat',
                                    candidate_rmlst_dict,
                                    primary)

    return candidate_rmlst_dict_results