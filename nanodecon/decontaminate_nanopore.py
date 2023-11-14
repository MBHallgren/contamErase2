import os
import sys
import gzip
from Bio import SeqIO
from nanodecon import kma
from nanodecon import util
from itertools import combinations


from nanodecon.intra_species_detection import determine_intra_species_contamination_nanopore
from nanodecon.nanopore_mutations import parse_sam_and_find_mutations

def nanopore_decontamination(arguments):
    os.system('mkdir ' + arguments.output)
    kma.KMARunner(arguments.nanopore,
                  arguments.output + "/bacteria_alignment",
                  arguments.db_dir + "/bac_db",
                  "-mem_mode -1t1 -t {} -ID 10 -ont -eq 14 -mct 0.5".format(arguments.threads)).run()

    total_bacteria_aligning_bases = util.number_of_bases_in_file(arguments.output + "/bacteria_alignment.fsa")
    primary, candidate_dict = drive_bacteria_results(arguments, total_bacteria_aligning_bases)
    primary_species = primary.split()[1] + ' ' + primary.split()[2]

    #TBD make this non-verbose
    headers = produce_species_specific_kma_db(primary_species,
                                    '/home/people/malhal/contamErase_db/rmlst.fsa',
                                    '/home/people/malhal/contamErase_db/rmlst_scheme.txt',
                                    arguments.output)

    kma.KMARunner(arguments.nanopore,
                  arguments.output + "/initial_rmlst_alignment",
                  arguments.output + '/specie_db',
                  "-t {} -ID 10 -ont -md 1.5 -matrix -eq 14 -mct 0.5".format(
                      arguments.threads, arguments.output)).run()

    #kma.KMARunner(arguments.nanopore,
    #              arguments.output + "/initial_rmlst_alignment",
    #              arguments.output + '/specie_db',
    #              "-mem_mode -t 16 -ID 90 -ont".format(
    #                  arguments.threads, arguments.output)).run()

    os.system('gunzip ' + arguments.output + '/initial_rmlst_alignment.frag.gz')

    extract_mapped_rmlst_read(arguments.output, arguments.nanopore)

    os.system('cat {}/rmlst_reads.fastq | NanoFilt -q 14 -l 500 > {}/trimmed_rmlst_reads.fastq'.format(arguments.output, arguments.output))

    index_top_hits_db(arguments.output)

    arguments.nanopore = arguments.output + '/trimmed_rmlst_reads.fastq'
    #

    kma.KMARunner(arguments.nanopore,
                  arguments.output + "/rmlst_alignment",
                  arguments.output + '/top_hits_db',
                  "-t {} -ID 10 -ont -md 1.5 -matrix -eq 14 -mct 0.5 -sam 2096> {}/rmlst_alignment.sam".format(arguments.threads, arguments.output)).run()
    #TMP removed -oa above
    os.system('gunzip ' + arguments.output + '/rmlst_alignment.mat.gz')

    consensus_dict = build_consensus_dict(arguments.output + '/rmlst_alignment.res',
                                          arguments.output + '/rmlst_alignment.mat')
    for item in consensus_dict:
        print (item, consensus_dict[item])
    sys.exit()

    consensus_dict, read_positions_blacklisted_dict = adjust_consensus_dict_for_individual_qscores(consensus_dict, arguments.output + '/rmlst_alignment.sam', arguments.output + '/trimmed_rmlst_reads.fastq')

    sys.exit()
    confirmed_mutation_dict = derive_mutation_positions2(consensus_dict, arguments)
    number = 0
    for item in confirmed_mutation_dict:
        number += len(confirmed_mutation_dict[item][0])
    print('Number of mutations: ' + str(number))

    #Consider this, can we exclude all novel mutations?
        #Can we do something to include novel mutations if the signal is strong enough?
        #Add significant but non biological mutations.
    confirmed_mutation_dict = validate_mutations(arguments, confirmed_mutation_dict, gene_score_dict, arguments.output + '/specie.fsa')

    #TBD: make this non-verbose
    number = 0
    for item in confirmed_mutation_dict:
        number += len(confirmed_mutation_dict[item][0])
    print ('Number of mutations: ' + str(number))

    confirmed_mutation_dict = upper_co_occuring_mutations_in_reads(arguments, confirmed_mutation_dict, gene_score_dict, arguments.output + '/specie.fsa', allele_pair_dict, consensus_dict)
    number = 0
    for item in confirmed_mutation_dict:
        number += len(confirmed_mutation_dict[item][0])
    print ('Number of mutations: ' + str(number))

    format_output(confirmed_mutation_dict, top_alleles, allele_pair_dict, consensus_dict)

    sys.exit()

def index_top_hits_db(output):
    infile = output + '/initial_rmlst_alignment.res'

    top_hits = {}
    with open(infile, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                allele = line[0]
                gene = allele.split('_')[0]
                if gene not in top_hits:
                    top_hits[gene] = [allele, line[1], 0]
                else:
                    if line[1] > top_hits[gene][1]:
                        top_hits[gene] = [allele, line[1], 0]

    name_file = output + '/specie_db.name'
    t = 1
    with open(name_file, 'r') as f:
        for line in f:
            line = line.strip()
            allele = line
            gene = allele.split('_')[0]
            if gene in top_hits:
                if top_hits[gene][0] == allele:
                    top_hits[gene][-1] = t
            t += 1

    seqs = ''
    for item in top_hits:
        seqs += str(top_hits[item][-1]) + ','
    seqs = seqs[:-1]
    os.system('kma seq2fasta -seqs {} -t_db {} > {}/top_hits.fsa'.format(seqs, output + '/specie_db', output))
    os.system('kma index -i {} -o {}' .format(output + '/top_hits.fsa', output + '/top_hits_db'))

def adjust_consensus_dict_for_individual_qscores(consensus_dict, sam_file, fastq_file):
    black_listed_positions = blacklist_positions(fastq_file, 14)

    adjusted_consensus_dict = {}

    for gene in consensus_dict:
        adjusted_consensus_dict[gene] = [[], '']
        derive_aligned_reads = derive_aligned_reads_for_gene(sam_file, gene)

def derive_aligned_reads_for_gene(sam_file, gene):
    aligned_reads = []
    for line in open(sam_file):
        if line[0] == '@':
            continue
        cols = line.strip().split('\t')
        qname, flag, rname, pos, mapq, cigar_str, rnext, pnext, tlen, seq = cols[:10]
        read_id = qname.split(' ')[0]
        gene_name = rname.split('_')[0]
        #if gene_name == gene:

    return aligned_reads

def blacklist_positions(fastq_file, quality_threshold):
    blacklist_dict = {}
    for record in SeqIO.parse(fastq_file, "fastq"):
        # Initialize the blacklist for this read
        blacklist = []

        # Check each quality score
        for pos, score in enumerate(record.letter_annotations["phred_quality"]):
            if score < quality_threshold:
                # Add the position to the blacklist (positions are 0-indexed)
                blacklist.append(pos)

        # Add the blacklist to the dictionary with the read ID as the key
        blacklist_dict[record.id] = blacklist

    return blacklist_dict

def format_output(confirmed_mutation_dict, top_alleles, allele_pair_dict, consensus_dict):
    header = 'Gene,MajorityAlelle,Position,MajorityBase,MutationBase,MutationDepth,TotalDepth'
    print (header)
    for gene in confirmed_mutation_dict:
        for mutation in zip(confirmed_mutation_dict[gene][0], confirmed_mutation_dict[gene][1]):
            position = mutation[0].split('_')[0]
            mutation_base = mutation[0].split('_')[1]
            mutation_depth = mutation[1]
            majority_base = consensus_dict[gene][1][int(position) - 1]
            total_depth = sum(consensus_dict[gene][0][int(position) - 1])
            print ('{},{},{},{},{},{},{}'.format(gene, allele_pair_dict[gene], position, majority_base, mutation_base, mutation_depth, total_depth))

def extract_mapped_rmlst_read(output, nanopore):
    read_set = set()
    with open(output + '/initial_rmlst_alignment.frag', 'r') as frag:
        for line in frag:
            line = line.rstrip()
            line = line.split('\t')
            read_set.add(line[-1])
    with open(output + '/rmlst_reads.txt', 'w') as f:
        for item in read_set:
            f.write(item + '\n')
    os.system('seqtk subseq {} {} > {}'.format(nanopore, output + '/rmlst_reads.txt',
                                               output + '/rmlst_reads.fastq'))


def realign_rmlst_to_hits(res_file, name_file):
    rmlst_alleles = set()
    rmlst_genes = dict()
    with open(res_file, 'r') as res:
        for line in res:
            line = line.rstrip()
            if line[0] != '#':
                allele = line.split()[0]
                gene = allele.split('_')[0]
                if gene not in rmlst_genes:
                    rmlst_genes[gene] = [line[1], allele]
                else:
                    if rmlst_genes[gene][0] < line[1]:
                        rmlst_genes[gene] = [line[1], allele]

    for item in rmlst_genes:
        rmlst_alleles.add(rmlst_genes[item][1])
    
    ids = []
    with open(name_file, 'r') as names:
        t = 1
        for line in names:
            line = line.rstrip()
            if line in rmlst_alleles:
                ids.append(t)
            t += 1
    
    seqs = ",".join(str(x) for x in ids)
    
    return seqs
    
            
    




def derive_mutation_positions2(consensus_dict, arguments):
    all_confirmed_mutation_dict = {}
    for item in consensus_dict:
        all_confirmed_mutation_dict[item] = [[], []]
        for i in range(len(consensus_dict[item][0])):
            positions = consensus_dict[item][0][i][:4]
            max_number = max(positions)
            index_of_max = positions.index(max_number)
            nucleotide_index = ['A', 'C', 'G', 'T']
            for t in range(len(positions)):
                if t != index_of_max: #Use this to parse majority calls in the future if needed
                    if positions[t] >= arguments.min_n:
                        total_depth = sum(positions)
                        relative_depth = positions[t] / total_depth
                        if relative_depth >= arguments.coc * arguments.mrd and positions[t] >= 2: #Only consider mutations with min depth >=2
                            all_confirmed_mutation_dict[item][0].append(
                                '{}_{}'.format(i + 1, nucleotide_index[t]))
                            all_confirmed_mutation_dict[item][1].append(positions[t])
    return all_confirmed_mutation_dict





def upper_co_occuring_mutations_in_reads(arguments, confirmed_mutation_dict, gene_score_dict, fsa_file, allele_pair_dict, consensus_dict):
    #TBD why not just get the mutation list from the confirmed_mutation_dict?
    #HERE
    reads_mutation_dict = parse_sam_and_find_mutations(arguments.output + '/rmlst_alignment.sam',
                                                       arguments.output + '/specie.fsa',
                                                       allele_pair_dict,
                                                       confirmed_mutation_dict,
                                                       arguments.output + '/rmlst_alignment.fsa',
                                                       consensus_dict)

    co_occurence_matrix_dict = {}
    for gene in confirmed_mutation_dict:
        mutation_list = confirmed_mutation_dict[gene][0]
        depth_list = confirmed_mutation_dict[gene][1]
        num_mutations = len(mutation_list)
        if num_mutations > 1:
            co_occurrence_matrix = [[0] * num_mutations for _ in range(num_mutations)]
            for read in reads_mutation_dict:
                read_gene = read.split(' ')[1].split('_')[0]
                if read_gene == gene:
                    read_mutations = reads_mutation_dict[read]
                    valid_mutations = [mutation for mutation in read_mutations if mutation in mutation_list]
                    if len(valid_mutations) > 1:
                        # Increment the matrix for pairs of mutations in the read
                        for i in range(len(valid_mutations)):
                            for j in range(i + 1, len(valid_mutations)):
                                mutation1 = mutation_list.index(valid_mutations[i])
                                mutation2 = mutation_list.index(valid_mutations[j])
                                co_occurrence_matrix[mutation1][mutation2] += 1
                                co_occurrence_matrix[mutation2][mutation1] += 1

            # Print the co-occurrence matrix with mutation names
            print ("gene:", gene)
            print("Mutation names:", mutation_list)
            #print ("Depth:", depth_list)
            average_depth = sum(confirmed_mutation_dict[gene][1]) / len(confirmed_mutation_dict[gene][1])
            print ("Average depth:", average_depth)
            print ("Threshold:", average_depth * 0.5) #Here, TBD look at threshold. Is more 0.5 really fine? Or should we do something similar to the benchmarking script
            for i, row in enumerate(co_occurrence_matrix):
                mutation_name = mutation_list[i]
                print(f"{mutation_name}: {row}")

            co_occurence_matrix_dict[gene] = co_occurrence_matrix

    adjusted_mutation_dict = {}
    for gene in co_occurence_matrix_dict:
        adjusted_mutation_dict[gene] = [[], []]
        average_depth = sum(confirmed_mutation_dict[gene][1]) / len(confirmed_mutation_dict[gene][1])
        threshold = average_depth * 0.5 # TBD reconsider
        if threshold < 3: # Can we justifiably set this to 3?
            threshold = 3
        for i, row in enumerate(co_occurence_matrix_dict[gene]):
            for number_of_co_occurences in row:
                total_depth = sum(consensus_dict[gene][0][i])
                relative_depth = confirmed_mutation_dict[gene][1][i] / total_depth
                if float(number_of_co_occurences) >= float(threshold):
                    adjusted_mutation_dict[gene][0].append(confirmed_mutation_dict[gene][0][i])
                    adjusted_mutation_dict[gene][1].append(confirmed_mutation_dict[gene][1][i])
                    break
                elif (relative_depth >= arguments.mrd):
                    adjusted_mutation_dict[gene][0].append(confirmed_mutation_dict[gene][0][i])
                    adjusted_mutation_dict[gene][1].append(confirmed_mutation_dict[gene][1][i])
                    break

    return adjusted_mutation_dict


    #pass

def determine_mutation_sets(reads_mutation_dict, mutation_position_dict):
    mutation_count_dict = {}
    for read in reads_mutation_dict:
        reference = reads_mutation_dict[read][1]
        read_mutation_set = set(reads_mutation_dict[read][0])
        gene_mutation_set = set(mutation_position_dict[reference][1])
        if read_mutation_set != set():
            common_elements = read_mutation_set & gene_mutation_set
            if common_elements != set():
                if reference + '_' + str(common_elements) not in mutation_count_dict:
                    mutation_count_dict[reference + '_' + str(common_elements)] = 1
                else:
                    mutation_count_dict[reference + '_' + str(common_elements)] += 1

    sorted_items = sorted(mutation_count_dict.items(), key=lambda item: item[1], reverse=True)

def validate_mutations(arguments, mutation_position_dict, gene_score_dict, fsa_file):
    correct_length_dict = derive_correct_length_headers(arguments, gene_score_dict, fsa_file)
    mutations_found_in_rmlst_genes = dict()
    for allele in mutation_position_dict:
        gene = allele.split('_')[0]
        mutations_found_in_rmlst_genes[gene] = set()
        for template_gene in correct_length_dict:
            if template_gene == gene:
                for template_allele in correct_length_dict[template_gene]:
                    #print (allele, template_allele)
                    confirmed_mutations = validate_mutation_positions(mutation_position_dict[allele][0], correct_length_dict[template_gene][template_allele], allele)
                    if confirmed_mutations != set():
                        if gene not in mutations_found_in_rmlst_genes:
                            mutations_found_in_rmlst_genes[gene] = confirmed_mutations
                        else:
                            mutations_found_in_rmlst_genes[gene] = mutations_found_in_rmlst_genes[gene] | confirmed_mutations
    return_dict = dict()
    for item in mutation_position_dict:
        return_dict[item] = [[], [], 0]
        for i in range(len(mutation_position_dict[item][0])):
            if mutation_position_dict[item][0][i] in mutations_found_in_rmlst_genes[item.split('_')[0]]:
                return_dict[item][0].append(mutation_position_dict[item][0][i])
                return_dict[item][1].append(mutation_position_dict[item][1][i])
        if return_dict[item][1] != []:
            return_dict[item][2] = min(return_dict[item][1])
    return return_dict

def validate_mutation_positions(mutations, sequence, allele):
    gene = allele.split('_')[0]
    confirmed_mutations = set()
    for mutation in mutations:
        position, mutant = mutation.split('_')
        if sequence[int(position) - 1] == mutant:
            confirmed_mutations.add(mutation)
    return confirmed_mutations
def derive_correct_length_headers(arguments, gene_score_dict, fsa_file):
    correct_length_dict = {}
    sequence = ''
    gene = None
    with open(fsa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if gene != None:
                    if gene in gene_score_dict:
                        if sequence != '':
                            if len(sequence) == gene_score_dict[gene][-1][0]:
                                if gene not in correct_length_dict:
                                    correct_length_dict[gene] = {}
                                correct_length_dict[gene][allele] = sequence
                header = line.strip()[1:]
                allele = header
                gene = allele.split('_')[0]
                sequence = ''
            else:
                sequence += line.strip()
    if gene != None:
        if gene in gene_score_dict:
            if sequence != '':
                if len(sequence) == gene_score_dict[gene][-1][0]:
                    if gene not in correct_length_dict:
                        correct_length_dict[gene] = {}
                    correct_length_dict[gene][allele] = sequence

    return correct_length_dict

def produce_species_specific_kma_db(species, fsa_file, scheme_file, output):
    gene_set = set()
    t = 0
    with open(scheme_file, 'r') as f:
        for line in f:
            if line.startswith('rST'):
                headers = line.strip().split('\t')[1:54]
            else:
                if line.strip().split('\t') != ['']:
                    if line.strip().split('\t')[55] == species:
                        t += 1
                        for i in range(len(headers)):
                            allele = headers[i] + '_' + line.strip().split('\t')[i+1]
                            gene_set.add(allele)
    produce_species_fsa_file(fsa_file, gene_set, output)
    os.system('kma index -i {}/specie.fsa -o {}/specie_db'.format(output, output))
    return headers

def produce_species_fsa_file(fsa_file, gene_set, output):
    with open(output + '/specie.fsa', 'w') as outfile:
        with open(fsa_file, 'r') as f:
            write_sequence = False  # A flag to indicate whether the sequence should be written to the output
            for line in f:
                if line.startswith('>'):
                    # Check if the gene_id (without '>') is in the gene_set
                    write_sequence = line.strip().split()[0][1:] in gene_set
                # Write the line (header or sequence) if write_sequence is True
                if write_sequence:
                    outfile.write(line)

def derive_mutation_positions(consensus_dict, fsa_file, headers, arguments, top_alleles):
    upper_confirmed_mutation_dict = {}
    lower_confirmed_mutation_dict = {}
    with open(fsa_file, 'r') as f:
        #TBD do we change the method of looking for mutations?
        #Main_call = max(list)
        #Mutations = any positions > t and < Main_call
        sequence = ''
        min_depth = 100000
        upper_mutation_list = []
        upper_mutation_depth = []
        lower_mutation_list = []
        lower_mutation_depth = []
        for line in f:
            if line.startswith('>'):
                if sequence != '':
                    if allele in top_alleles:
                        upper_mutation_list = []
                        upper_mutation_depth = []
                        lower_mutation_list = []
                        lower_mutation_depth = []
                        for i in range(len(sequence)):
                            if sequence[i] == 'A':
                                index = 0
                            elif sequence[i] == 'C':
                                index = 1
                            elif sequence[i] == 'G':
                                index = 2
                            elif sequence[i] == 'T':
                                index = 3
                            else:
                                index = 4
                                sys.exit('Check here')
                            #relative_depth = consensus_dict[allele][i][index] / total_base_count
                            nucleotide_index = ['A', 'C', 'G', 'T']
                            a_depth = consensus_dict[gene][i][0]
                            c_depth = consensus_dict[gene][i][1]
                            g_depth = consensus_dict[gene][i][2]
                            t_depth = consensus_dict[gene][i][3]
                            depths = [a_depth, c_depth, g_depth, t_depth]
                            for t in range(len(depths)):
                                if t != index:
                                    if depths[t] >= arguments.min_n:
                                        total_depth = sum(consensus_dict[gene][i][:4])
                                        relative_depth = depths[t] / total_depth
                                        if relative_depth >= arguments.urd:
                                            upper_mutation_list.append(
                                                '{}_{}_{}'.format(i + 1, sequence[i], nucleotide_index[t]))
                                            upper_mutation_depth.append(depths[t])
                                        elif relative_depth >= arguments.lrd:
                                            lower_mutation_list.append(
                                                '{}_{}_{}'.format(i + 1, sequence[i], nucleotide_index[t]))
                                            lower_mutation_depth.append(depths[t])


                        if upper_mutation_list != []:
                            upper_confirmed_mutation_dict[allele] = [min(upper_mutation_depth), upper_mutation_list, upper_mutation_depth, total_depth]
                        else:
                            upper_confirmed_mutation_dict[allele] = [0, [], [], total_depth]
                        if lower_mutation_list != []:
                            lower_confirmed_mutation_dict[allele] = [min(lower_mutation_depth), lower_mutation_list, lower_mutation_depth, total_depth]
                        else:
                            lower_confirmed_mutation_dict[allele] = [0, [], [], total_depth]

                allele = line.strip()[1:]
                gene = allele.split('_')[0]
                sequence = ''
                min_depth = 100000
            else:
                sequence += line.strip()
    if allele in top_alleles:
        upper_mutation_list = []
        upper_mutation_depth = []
        lower_mutation_list = []
        lower_mutation_depth = []
        for i in range(len(sequence)):
            if sequence[i] == 'A':
                index = 0
            elif sequence[i] == 'C':
                index = 1
            elif sequence[i] == 'G':
                index = 2
            elif sequence[i] == 'T':
                index = 3
            else:
                index = 4
                sys.exit('Check here')
            # relative_depth = consensus_dict[allele][i][index] / total_base_count
            nucleotide_index = ['A', 'C', 'G', 'T']
            a_depth = consensus_dict[gene][i][0]
            c_depth = consensus_dict[gene][i][1]
            g_depth = consensus_dict[gene][i][2]
            t_depth = consensus_dict[gene][i][3]
            depths = [a_depth, c_depth, g_depth, t_depth]
            for t in range(len(depths)):
                if t != index:
                    if depths[t] >= arguments.min_n:
                        total_depth = sum(consensus_dict[gene][i][:4])
                        relative_depth = depths[t] / total_depth
                        if relative_depth >= arguments.urd:
                            upper_mutation_list.append(
                                '{}_{}_{}'.format(i + 1, sequence[i], nucleotide_index[t]))
                            upper_mutation_depth.append(depths[t])
                        elif relative_depth >= arguments.lrd:
                            lower_mutation_list.append(
                                '{}_{}_{}'.format(i + 1, sequence[i], nucleotide_index[t]))
                            lower_mutation_depth.append(depths[t])

        if upper_mutation_list != []:
            upper_confirmed_mutation_dict[allele] = [min(upper_mutation_depth), upper_mutation_list, upper_mutation_depth, total_depth]
        else:
            upper_confirmed_mutation_dict[allele] = [0, [], [], total_depth]
        if lower_mutation_list != []:
            lower_confirmed_mutation_dict[allele] = [min(lower_mutation_depth), lower_mutation_list, lower_mutation_depth, total_depth]
        else:
            lower_confirmed_mutation_dict[allele] = [0, [], [], total_depth]

    return upper_confirmed_mutation_dict, lower_confirmed_mutation_dict
def check_all_species_alleles_against_consensus_dict(consensus_dict, fsa_file, headers, arguments):
    confirmed_alleles = {}
    if arguments.min_n != None:
        threshold = arguments.min_n
    else:
        threshold = None
    with open(fsa_file, 'r') as f:
        #TBD Doesn't work for #2 gene, template 48 not in resfile
        sequence = ''
        min_depth = 100000
        for line in f:
            if line.startswith('>'):
                if sequence != '':
                    if len(sequence) == len(consensus_dict[allele]):
                        mutation_list = []
                        mutation_depth = []
                        for i in range(len(sequence)):
                            if sequence[i] == 'A':
                                index = 0
                            elif sequence[i] == 'C':
                                index = 1
                            elif sequence[i] == 'G':
                                index = 2
                            elif sequence[i] == 'T':
                                index = 3
                            else:
                                index = 4
                                sys.exit('Check here')
                            #relative_depth = consensus_dict[allele][i][index] / total_base_count
                            nucleotide_index = ['A', 'C', 'G', 'T']
                            a_depth = consensus_dict[allele][i][0]
                            c_depth = consensus_dict[allele][i][1]
                            g_depth = consensus_dict[allele][i][2]
                            t_depth = consensus_dict[allele][i][3]
                            depths = [a_depth, c_depth, g_depth, t_depth]
                            for t in range(len(depths)):
                                if t != index:
                                    if threshold == None:
                                        total_depth = sum(consensus_dict[allele][i][:4])
                                        relative_depth = depths[t] / total_depth
                                        if relative_depth > arguments.min_rd:
                                            mutation_list.append(
                                                '{}_{}_{}'.format(i + 1, sequence[i], nucleotide_index[t]))
                                            mutation_depth.append(depths[t])
                                    else:
                                        if depths[t] > threshold:
                                            total_depth = sum(consensus_dict[allele][i][:4])
                                            relative_depth = depths[t] / total_depth
                                            if relative_depth > arguments.min_rd:
                                                mutation_list.append('{}_{}_{}'.format(i + 1, sequence[i], nucleotide_index[t]))
                                                mutation_depth.append(depths[t])

                        if mutation_list != []:
                            confirmed_alleles[gene] = [min(mutation_depth), mutation_list, mutation_depth]

                    sequence = ''
                    min_depth = 100000
                gene = line.strip()[1:]
                allele = gene.split('_')[0]
            else:
                sequence += line.strip()
    if mutation_list != []:
        confirmed_alleles[gene] = [min(mutation_depth), mutation_list, mutation_depth]

    final_allleles = {}
    for gene in headers:
        top_score = 0
        top_allele = ''
        for allele in confirmed_alleles:
            if allele.startswith(gene):
                if confirmed_alleles[allele][0] > top_score:
                    top_score = confirmed_alleles[allele][0]
                    top_allele = allele
        if top_allele != '':
            final_allleles[top_allele] = confirmed_alleles[top_allele]

    for gene in headers:
        contenders = {}
        for allele in confirmed_alleles:
            if allele.startswith(gene):
                if allele not in final_allleles:
                    contenders[allele] = confirmed_alleles[allele]
        most_mutated_allele = []
        most_mutated_allele_score = 0
        for allele in contenders:
            if len(contenders[allele][1]) > most_mutated_allele_score:
                most_mutated_allele_score = len(contenders[allele][1])
                most_mutated_allele = [allele]
            elif len(contenders[allele][1]) == most_mutated_allele_score:
                most_mutated_allele.append(allele)

        if len(most_mutated_allele) == 1:
            # TBD Check depths of mutation positions
            final_allleles[most_mutated_allele[0]] = confirmed_alleles[most_mutated_allele[0]]
            for allele in contenders:
                if allele != most_mutated_allele[0]:
                    if not is_subset(contenders[allele][1], confirmed_alleles[most_mutated_allele[0]][1]):
                        final_allleles[allele] = confirmed_alleles[allele]

        elif len(most_mutated_allele) > 1:
            #TBD Check depths of mutation positions
            for allele in most_mutated_allele:
                final_allleles[allele] = confirmed_alleles[allele]
            for allele in contenders:
                if allele not in most_mutated_allele:
                    not_subset = True
                    for most_mutated in most_mutated_allele:
                        if is_subset(contenders[allele][1], confirmed_alleles[most_mutated][1]):
                            not_subset = False
                    if not_subset:
                        final_allleles[allele] = confirmed_alleles[allele]
    return final_allleles

def is_subset(list_a, list_b):
    for element in list_a:
        if element not in list_b:
            return False
    return True

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


def build_consensus_dict(res_file, mat_file):
    consensus_dict = {}
    with open(res_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                allele = line[0]
                consensus_dict[allele] = [[], '']
                consensus_dict[allele][0].append([0, 0, 0, 0, 0, 0]) #[A, C, G, T, N, -]

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
                                consensus_dict[gene][0][index][i] += int(line[i])
                            index += 1
    for gene in consensus_dict:
        for position in consensus_dict[gene][0]:
            consensus_dict[gene][1] += 'ACGTN-'[position.index(max(position))]
    return consensus_dict


def produce_final_output_nanopore(arguments, frag_file, primary, candidate_rmlst_dict_results, black_list_plasmid, black_list_viral, black_list_human):
    rmlst_hits = output_primary_reads_nanopore(arguments, frag_file, primary, candidate_rmlst_dict_results, black_list_plasmid)
    if len(rmlst_hits) >= 1:
        output_rmlst_hits_nanopore(arguments, frag_file, rmlst_hits)
    if len(black_list_plasmid) > 0:
        output_black_list_plasmid_nanopore(arguments, black_list_plasmid)
    if len(black_list_viral) > 0:
        output_black_list_viral_nanopore(arguments, black_list_viral)
    if len(black_list_human) > 0:
        output_black_list_human_nanopore(arguments, black_list_human)
    clean_up(arguments)

def output_black_list_human_nanopore(arguments, black_list_human):
    with open(arguments.output + 'human_remove.txt', 'w') as write_file:
        for read in black_list_human:
            write_file.write(read + '\n')
    os.system('seqtk subseq {} {} > {}'.format(arguments.nanopore, arguments.output + 'human_remove.txt', arguments.output + '/output/contaminations/human_reads.fastq'))
def output_black_list_viral_nanopore(arguments, black_list_viral):
    with open(arguments.output + 'viral_remove.txt', 'w') as write_file:
        for read in black_list_viral:
            write_file.write(read + '\n')
    os.system('seqtk subseq {} {} > {}'.format(arguments.nanopore, arguments.output + 'viral_remove.txt', arguments.output + '/output/contaminations/viral_reads.fastq'))

def output_black_list_plasmid_nanopore(arguments, black_list_plasmid):
    with open(arguments.output + 'plasmid_remove.txt', 'w') as write_file:
        for read in black_list_plasmid:
            write_file.write(read + '\n')
    os.system('seqtk subseq {} {} > {}'.format(arguments.nanopore, arguments.output + 'plasmid_remove.txt', arguments.output + '/output/contaminations/plasmid_reads.fastq'))

def clean_up(arguments):
    os.system('rm -r ' + arguments.output + '/*.fastq')
    os.system('mv ' + arguments.output + '/*.frag ' + arguments.output + '/alignments/.')
    os.system('mv ' + arguments.output + '/*.fsa ' + arguments.output + '/alignments/.')
    os.system('mv ' + arguments.output + '/*.res ' + arguments.output + '/alignments/.')
    os.system('mv ' + arguments.output + '/*.aln ' + arguments.output + '/alignments/.')
    os.system('tar -czvf ' + arguments.output + '/alignments.tar.gz ' + arguments.output + '/alignments')
    os.system('rm -r ' + arguments.output + '/alignments')

def output_rmlst_hits_nanopore(arguments, frag_file, rmlst_hits):
    intra_species_reads = set()
    with open(frag_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[-2] in rmlst_hits:
                intra_species_reads.add(line[-1])
    print (len(intra_species_reads), 'intra species reads')
    with open('{}/{}.txt'.format(arguments.output, 'intra_species_reads'), 'w') as write_file:
        for read in intra_species_reads:
            write_file.write(read + '\n')
    os.system('seqtk subseq {} {} > {}'.format(arguments.nanopore, arguments.output + '/{}.txt'.format('intra_species_reads'), arguments.output + '/output/contaminations/intra_species_reads.fastq'))

def output_primary_reads_nanopore(arguments, frag_file, primary, candidate_rmlst_dict_results, black_list_plasmids):
    black_list_templates = set()
    if candidate_rmlst_dict_results is not None:
        for item in candidate_rmlst_dict_results:
            if item != primary and candidate_rmlst_dict_results[item]['rmlst_id'] == True:
                black_list_templates.add(item)

        read_set = set()
        with open(frag_file, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if line[-2] not in black_list_templates:
                    read_set.add(line[-1])
        for item in black_list_plasmids:
            read_set.add(item)
    else:
        read_set = set()
        with open(frag_file, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                read_set.add(line[-1])
        for item in black_list_plasmids:
            read_set.add(item)

    with open('{}/{}.txt'.format(arguments.output, 'primary_reads_list'), 'w') as write_file:
        for read in read_set:
            write_file.write(read + '\n')
    os.system('seqtk subseq {} {} > {}'.format(arguments.nanopore, arguments.output + '/{}.txt'.format('primary_reads_list'), arguments.output + '/output/primary_reads.fastq'))

    return black_list_templates


def derive_all_typed_rmlst_genes_from_input_reads(filename):
    all_rmlst_genes = dict()
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                all_rmlst_genes[line.split('\t')[0]] = line.split('\t')[8]
    return all_rmlst_genes

def number_of_bases_in_file(filename):
    gzipped, type = determine_file_type(filename)
    #print (gzipped, type)
    if type == 'fasta':
        sum = 0
        with open(filename, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sum += len(line.strip())
        return sum

    elif type == 'fastq':
        if gzipped:
            line_count = 1
            sum = 0
            with gzip.open(filename, 'r') as f:
                for line in f:
                    if line_count == 2:
                        sum += len(line.strip())
                    line_count += 1
                    if line_count == 5:
                        line_count = 1
            return sum
        else:
            line_count = 1
            sum = 0
            with open(filename, 'r') as f:
                for line in f:
                    if line_count == 2:
                        sum += len(line.strip())
                    line_count += 1
                    if line_count == 5:
                        line_count = 1
            return sum

def derive_non_bacterial_black_list(output):
    bac_list = set()
    plasmid_list = set()
    viral_list = set()
    human_list = set()
    os.system('gunzip {}/bacteria_alignment.frag.gz'.format(output))
    with open('{}/bacteria_alignment.frag'.format(output), 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            bac_list.add(line[-1])
    os.system('gunzip {}/plasmid_alignment.frag.gz'.format(output))
    with open ('{}/plasmid_alignment.frag'.format(output), 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            plasmid_list.add(line[-1])
    os.system('gunzip {}/viral_alignment.frag.gz'.format(output))
    with open ('{}/viral_alignment.frag'.format(output), 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            viral_list.add(line[-1])
    os.system('gunzip {}/human_alignment.frag.gz'.format(output))
    with open ('{}/human_alignment.frag'.format(output), 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            human_list.add(line[-1])
    black_list_plasmid = plasmid_list - bac_list #Remove bac from plasmid
    black_list_viral = viral_list - bac_list #Remove bac from viral
    black_list_human = human_list - bac_list #Remove bac from human
    return black_list_plasmid, black_list_viral, black_list_human

def eval_bacteria_results(results, total_bacteria_aligning_bases):
    primary = results[0]['#Template']
    candidate_dict = dict()
    for item in results:
        bases_hit = int(item['Template_length']) * float(item['Depth'].strip())
        template_id = float(item['Template_Identity'])
        relative_template_depth = bases_hit / total_bacteria_aligning_bases
        if relative_template_depth > 0.01 or template_id > 20.00:
            candidate_dict[item['#Template']] = [relative_template_depth, int(item['Template_length'])]
    return primary, candidate_dict

def drive_bacteria_results(arguments, total_bacteria_aligning_bases):
    results = sort_lines_by_score(arguments.output + "/bacteria_alignment.res")
    primary, candidate_dict = eval_bacteria_results(results, total_bacteria_aligning_bases)
    return primary, candidate_dict

def derive_rmlst_candidates(primary, candidate_dict):
    rmlst_candidate_set = set()
    primary_species = primary.split()[1] + ' ' + primary.split()[2]
    for item in candidate_dict:
        species = item.split()[1] + ' ' + item.split()[2]
        if species == primary_species:
            rmlst_candidate_set.add(item)
    return rmlst_candidate_set


def derive_read_pools_nanopore(candidate_dict, arguments, primary, black_list_human):
    primary_species = primary.split()[1] + ' ' + primary.split()[2]
    pool_dict = dict()
    for item in candidate_dict:
        species = item.split()[1] + ' ' + item.split()[2]
        if species not in pool_dict:
            pool_dict[species] = set()

    with open('{}/bacteria_alignment.frag'.format(arguments.output), 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            line_species = line[-2].split()[1] + ' ' + line[-2].split()[2]
            if line_species in pool_dict:
                pool_dict[line_species].add(line[-1])

    for item in pool_dict:
        if item != primary_species:
            with open('{}/{}.txt'.format(arguments.output, item.replace(' ', '_').lower()), 'w') as write_file:
                for read in pool_dict[item]:
                    write_file.write('{}\n'.format(read))
            os.system('seqtk subseq {} {}/{}.txt > {}/{}.fastq'.format(arguments.nanopore, arguments.output, item.replace(' ', '_').lower(), arguments.output, item.replace(' ', '_').lower()))
            for read in pool_dict[item]:
                black_list_human.add(read)

    primary_species_reads = set()
    with open('{}'.format(arguments.nanopore), 'r') as f:
        for line in f:
            if line[0] == '@':
                read = line.strip()[1:]
                if read not in black_list_human:
                    primary_species_reads.add(read)

    with open('{}/remove.txt'.format(arguments.output), 'w') as write_file:
        with open('{}'.format(arguments.nanopore), 'r') as f:
            for line in f:
                if line[0] == '@':
                    if line.strip()[1:] not in black_list_human:
                        write_file.write(line.strip()[1:] + '\n')
    os.system('seqtk subseq {} {}/remove.txt > {}/{}_read_pool.fastq'.format(arguments.nanopore, arguments.output, arguments.output, primary_species.replace(' ', '_').lower()))

def determine_file_type(file):
    gzipped = False
    type = None
    if file.endswith('.gz'):
        gzipped = True
        file = file[:-3]
    if file.endswith('.fastq') or file.endswith('.fq'):
        type = 'fastq'
    elif file.endswith('.fasta') or file.endswith('.fa') or file.endswith('.fna') or file.endswith('.fsa'):
        type = 'fasta'
    return gzipped, type

def sort_lines_by_score(filename):
    data = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        header = lines[0].strip().split('\t')
        for line in lines[1:]:
            values = line.strip().split('\t')
            data.append(dict(zip(header, values)))
    sorted_data = sorted(data, key=lambda x: int(x['Score']), reverse=True)
    return sorted_data
