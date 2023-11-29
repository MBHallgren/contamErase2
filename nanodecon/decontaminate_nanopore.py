import os
import sys
import gzip
from Bio import SeqIO
from nanodecon import kma
from nanodecon import util
from itertools import combinations


from nanodecon.intra_species_detection import determine_intra_species_contamination_nanopore
from nanodecon.nanopore_mutations import parse_sam_and_find_mutations
from nanodecon.nanopore_mutations import extract_alignment
from nanodecon.nanopore_mutations import create_mutation_vector
from nanodecon.nanopore_mutations import identify_mutations

def nanopore_decontamination(arguments):
    check_arguments(arguments)
    os.system('mkdir ' + arguments.output)
    if arguments.nanopore.endswith('.gz'):
        os.system('zcat {} | NanoFilt -q {} > {}/trimmed_reads.fastq'.format(arguments.nanopore, 14, arguments.output))
    else:
        os.system('cat {} | NanoFilt -q {} > {}/trimmed_reads.fastq'.format(arguments.nanopore, 14, arguments.output))
    arguments.nanopore = arguments.output + '/trimmed_reads.fastq'


    kma.KMARunner(arguments.nanopore,
                  arguments.output + "/bacteria_alignment",
                  arguments.db_dir + "/bac_db",
                  "-mem_mode -1t1 -t {} -ID 10 -ont -eq 14 -mct 0.5".format(arguments.threads)).run()

    total_bacteria_aligning_bases = util.number_of_bases_in_file(arguments.output + "/bacteria_alignment.fsa")
    primary, candidate_dict = drive_bacteria_results(arguments, total_bacteria_aligning_bases)
    primary_species = primary.split()[1] + ' ' + primary.split()[2]

    produce_species_specific_kma_db(primary_species,
                                    '/home/people/malhal/contamErase_db/rmlst.fsa',
                                    '/home/people/malhal/contamErase_db/rmlst_scheme.txt',
                                    arguments.output)

    kma.KMARunner(arguments.nanopore,
                  arguments.output + "/initial_rmlst_alignment",
                  arguments.output + '/specie_db',
                  "-t {} -ID 10 -ont -md 1.5 -matrix -eq 14 -mct 0.5".format(arguments.threads, arguments.output)).run()

    os.system('gunzip ' + arguments.output + '/initial_rmlst_alignment.frag.gz')

    extract_mapped_rmlst_read(arguments.output, arguments.nanopore)

    index_top_hits_db(arguments.output)

    arguments.nanopore = arguments.output + '/trimmed_rmlst_reads.fastq'

    kma.KMARunner(arguments.nanopore,
                  arguments.output + "/rmlst_alignment",
                  arguments.output + '/top_hits_db',
                  "-t {} -ID 10 -ont -md 1.5 -matrix -eq 14 -mct 0.5 -sam 2096> {}/rmlst_alignment.sam".format(arguments.threads, arguments.output)).run()
    os.system('gunzip ' + arguments.output + '/rmlst_alignment.mat.gz')

    consensus_dict = build_consensus_dict(arguments.output + '/rmlst_alignment.res',
                                          arguments.output + '/rmlst_alignment.mat')

    consensus_dict, read_positions_blacklisted_dict = adjust_consensus_dict_for_individual_qscores(consensus_dict, arguments.output + '/rmlst_alignment.sam', arguments.nanopore)

    confirmed_mutation_dict = derive_mutation_positions(consensus_dict, arguments)

    bio_validation_dict = bio_validation_mutations(consensus_dict, arguments.output + '/specie.fsa', confirmed_mutation_dict)

    confirmed_mutation_dict = co_occurence_until_convergence(arguments, confirmed_mutation_dict, consensus_dict, read_positions_blacklisted_dict, bio_validation_dict)
    format_output(confirmed_mutation_dict, consensus_dict, bio_validation_dict)

    sys.exit()

def co_occurence_until_convergence(arguments, confirmed_mutation_dict, consensus_dict, read_positions_blacklisted_dict, bio_validation_dict):
    """
    Iteratively call upper_co_occuring_mutations_in_reads until no new mutations are found.
    """

    current_count = count_mutations_in_mutations_dict(confirmed_mutation_dict)
    iteration_count = 0
    while True:
        confirmed_mutation_dict =\
            upper_co_occuring_mutations_in_reads(arguments, confirmed_mutation_dict, consensus_dict,read_positions_blacklisted_dict, bio_validation_dict)
        new_count = count_mutations_in_mutations_dict(confirmed_mutation_dict)
        iteration_count += 1
        if new_count == current_count:
            break
        current_count = new_count
    return confirmed_mutation_dict

def check_arguments(arguments):
    if arguments.cor > 1:
        print ('cor must be between 0 and 1, otherwise it is not a reward.')
        sys.exit()
    if arguments.bp < 1:
        print ('bp must be greater than 1, otherwise it is not a penalty.')
        sys.exit()
    if arguments.pp < 1:
        print ('pp must be greater than 1, otherwise it is not a penalty.')
        sys.exit()


def count_mutations_in_mutations_dict(mutation_dict):
    count = 0
    for item in mutation_dict:
        count += len(mutation_dict[item][0])
    return count

def bio_validation_mutations(consensus_dict, fsa_file, confirmed_mutation_dict):
    correct_length_dict = derive_correct_length_headers(consensus_dict, fsa_file)
    mutation_dict = dict()
    for gene in correct_length_dict:
        mutation_dict[gene] = set()
        for sequence in correct_length_dict[gene][1]:
            for i in range(len(sequence)):
                mutation_dict[gene].add(str(i+1) + '_' + sequence[i])
    return mutation_dict

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
    os.system('kma index -i {} -o {} 2>/dev/null' .format(output + '/top_hits.fsa', output + '/top_hits_db'))

def adjust_consensus_dict_for_individual_qscores(consensus_dict, sam_file, fastq_file):
    #Both of the below are in the consensus dict
    black_listed_positions = blacklist_positions(fastq_file, 2) #TBD fix this. Bases in the start of genes are removed with a higher value
    #TBD also check the last position in the consensus dict. all of them return 6 0's right now.
    adjusted_consensus_dict = {}

    for allele in consensus_dict:
        adjusted_consensus_dict[allele] = [[], consensus_dict[allele][1]]
        for position in consensus_dict[allele][0]:
            adjusted_consensus_dict[allele][0].append([0, 0, 0, 0, 0, 0])

    total_black_list_count = 0
    with open(sam_file, 'r') as sam_file:
        for alignment in sam_file:
            if alignment[0] == '@':
                continue
            cols = alignment.strip().split('\t')
            qname, flag, rname, pos, mapq, cigar_str, rnext, pnext, tlen, seq = cols[:10]
            read_id = qname.split(' ')[0]
            template_seq = adjusted_consensus_dict[rname][1]
            pos = int(pos)
            tlen = int(tlen)

            if pos == 1 and len(seq) >= tlen: #We will only consider read that span the entire gene.
                # Obtaining the alignment using your function
                aligned_ref, aligned_query = extract_alignment(template_seq[pos - 1:pos - 1 + tlen], seq, cigar_str)
                mutation_vector = create_mutation_vector(aligned_ref, aligned_query)
                for i in range(len(mutation_vector)):
                    nucleotide_list = ['A', 'C', 'G', 'T', 'N', '-']
                    if i not in black_listed_positions[read_id]:
                        adjusted_consensus_dict[rname][0][i][nucleotide_list.index(mutation_vector[i])] += 1
                    else:
                        total_black_list_count += 1
    #print ('Total blacklisted positions: ' + str(total_black_list_count))
    return adjusted_consensus_dict, black_listed_positions

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

def format_output(confirmed_mutation_dict, consensus_dict, bio_validation_dict):
    header = 'Gene,MajorityAlelle,Position,MajorityBase,MutationBase,MutationDepth,TotalDepth,MutationComment'
    print (header)
    for allele in confirmed_mutation_dict:
        for mutation in zip(confirmed_mutation_dict[allele][0], confirmed_mutation_dict[allele][1]):
            position = mutation[0].split('_')[0]
            mutation_base = mutation[0].split('_')[1]
            mutation_depth = mutation[1]
            majority_base = consensus_dict[allele][1][int(position) - 1]
            total_depth = sum(consensus_dict[allele][0][int(position) - 1])
            biological_existance = check_single_mutation_exisistance(bio_validation_dict, allele, mutation[0])
            if biological_existance:
                print ('{},{},{},{},{},{}. {}'.format(allele, position, majority_base, mutation_base, mutation_depth, total_depth, 'Previously seen mutation in the database.'))
            else:
                print ('{},{},{},{},{},{}. {}'.format(allele, position, majority_base, mutation_base, mutation_depth, total_depth, 'Novel mutation, never seen before in the database.'))

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
                                               output + '/trimmed_rmlst_reads.fastq'))


def derive_mutation_positions(consensus_dict, arguments):
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
                        if relative_depth >= arguments.cor * arguments.mrd: #Only consider mutations with min depth >=2
                            all_confirmed_mutation_dict[item][0].append(
                                '{}_{}'.format(i + 1, nucleotide_index[t]))
                            all_confirmed_mutation_dict[item][1].append(positions[t])
    return all_confirmed_mutation_dict





def upper_co_occuring_mutations_in_reads(arguments, confirmed_mutation_dict, consensus_dict, read_positions_blacklisted_dict, bio_validation_dict):
    #TBD why not just get the mutation list from the confirmed_mutation_dict?
    #HERE

    reads_mutation_dict = parse_sam_and_find_mutations(arguments.output + '/rmlst_alignment.sam',
                                                       confirmed_mutation_dict,
                                                       consensus_dict,
                                                       read_positions_blacklisted_dict)

    co_occurence_matrix_dict = {}
    for allele in confirmed_mutation_dict:
        mutation_list = confirmed_mutation_dict[allele][0]
        depth_list = confirmed_mutation_dict[allele][1]
        num_mutations = len(mutation_list)
        if num_mutations > 1:
            co_occurrence_matrix = [[0] * num_mutations for _ in range(num_mutations)]
            for read in reads_mutation_dict:
                read_allele = read.split(' ')[1]
                if read_allele == allele:
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
            #print ("allele:", allele)
            #print("Mutation names:", mutation_list)
            #print ("Depth:", depth_list)
            average_depth = sum(confirmed_mutation_dict[allele][1]) / len(confirmed_mutation_dict[allele][1])
            #positional_depth = sum(consensus_dict[allele][0][0]) / len(consensus_dict[allele][0][0])
            #total_gene_depth = 0
            #for i in range(len(consensus_dict[allele][0])):
            #    total_gene_depth += sum(consensus_dict[allele][0][i])
            #average_depth = total_gene_depth / len(consensus_dict[allele][0])
            #threshold = average_depth * arguments.mrd * arguments.cor
            #print ("Threshold:", average_depth * 0.5) #Here, TBD look at threshold. Is more 0.5 really fine? Or should we do something similar to the benchmarking script
            #for i, row in enumerate(co_occurrence_matrix):
            #    mutation_name = mutation_list[i]
            #    print(f"{mutation_name}: {row} {check_single_mutation_exisistance(bio_validation_dict, allele, mutation_name)}")

            co_occurence_matrix_dict[allele] = [co_occurrence_matrix, mutation_list]


    #TBD continue here.
    adjusted_mutation_dict = {}
    for allele in confirmed_mutation_dict:
        if allele in co_occurence_matrix_dict:
            adjusted_mutation_dict[allele] = [[], []]
            matrix = co_occurence_matrix_dict[allele][0]
            mutation_list = co_occurence_matrix_dict[allele][1]
            depth_list = confirmed_mutation_dict[allele][1]
            for i in range(len(matrix)):
                row = matrix[i]
                mutation = mutation_list[i]
                position = int(mutation.split('_')[0])
                position_depth = sum(consensus_dict[allele][0][position - 1])
                mutation_depth = depth_list[i]
                proxi_mutations = find_mutations_proximity_specific_mutation(mutation_list, mutation, arguments.proxi
                biological_existance = check_single_mutation_exisistance(bio_validation_dict, allele, mutation)


                mutation_threshold = position_depth * arguments.mrd
                co_occurence_list = check_mutation_co_occurrence(row, mutation_list, mutation,
                                                                 position_depth, arguments.cor, arguments.pp, arguments.mrd, proxi_mutations)

                if co_occurence_list != []:
                    mutation_threshold = mutation_threshold * arguments.cor
                if not biological_existance:
                    mutation_threshold = mutation_threshold + (arguments.bp-1) * position_depth * arguments.mrd
                if proxi_mutations != []:
                    mutation_threshold = mutation_threshold + (arguments.pp-1) * position_depth * arguments.mrd
                if mutation_depth >= mutation_threshold:
                    adjusted_mutation_dict[allele][0].append(confirmed_mutation_dict[allele][0][i])
                    adjusted_mutation_dict[allele][1].append(confirmed_mutation_dict[allele][1][i])


        else:
            adjusted_mutation_dict[allele] = [[], []]
            if confirmed_mutation_dict[allele][0] != []:
                mutation = confirmed_mutation_dict[allele][0][0]
                position = int(mutation.split('_')[0])
                position_depth = sum(consensus_dict[allele][0][position - 1])
                mutation_threshold = position_depth * arguments.mrd
                depth = confirmed_mutation_dict[allele][1][0]
                biological_existance = check_single_mutation_exisistance(bio_validation_dict, allele, mutation)
                if not biological_existance:
                    mutation_threshold = mutation_threshold + (arguments.bp-1) * position_depth * arguments.mrd

                #Check for mutation_threshold
                if depth >= mutation_threshold:
                    #print ('Added single mutation_2: ', allele, mutation)
                    adjusted_mutation_dict[allele][0].append(confirmed_mutation_dict[allele][0][0])
                    adjusted_mutation_dict[allele][1].append(confirmed_mutation_dict[allele][1][0])
    return adjusted_mutation_dict

def check_mutation_co_occurrence(list_of_mutation_co_occurrence, mutation_list, mutation,
                                 position_depth, cor, pp, mrd, proxi_mutations):
    """
    Checks if a given mutation co-occurs with another mutation above a specified threshold.

    :param list_of_mutation_co_occurrence: List of co-occurrence counts for each mutation.
    :param threshold: The threshold value for determining significant co-occurrence.
    :param mutation_list: List of mutations corresponding to the co-occurrence counts.
    :param mutation: The mutation to check for co-occurrence.
    :return: True if the mutation co-occurs with any other mutation above the threshold, False otherwise.
    """
    if mutation not in mutation_list:
        #Should never happen
        return [] #no co-occurence and not in proximity

    co_threshold = position_depth * mrd * cor #default co_threshold
    if co_threshold < 3:
        co_threshold = 3

    # Find the index of the mutation in the mutation list
    mutation_index = mutation_list.index(mutation)

    co_occurence_list = []
    # Check if the co-occurrence count of the mutation with any other mutation is above the threshold
    for i, count in enumerate(list_of_mutation_co_occurrence):
        if mutation_list[i] in proxi_mutations:
            # Add penalty for proximity to make it harder to get the co-occurence reward
            # for mutations within the proximity
            co_threshold = co_threshold * pp
        if i != mutation_index and count >= co_threshold:
            co_occurence_list.append(mutation_list[i])

    return co_occurence_list

def check_single_mutation_exisistance(bio_validation_dict, allele, specific_mutation):
    gene = allele.split('_')[0]
    if specific_mutation in bio_validation_dict[gene]:
        return True
    return False

def find_mutations_proximity_specific_mutation(mutations, specific_mutation, proxi):
    """
    Find mutations that are within a specified proximity of a specific mutation.

    :param mutations: A list of mutations.
    :param specific_mutation: The specific mutation to compare others against.
    :param proxi: The proximity within which another mutation can't be.
    :return: A list of mutations which are within the specified proximity of the specific mutation.
    """
    specific_mutation_pos = int(specific_mutation.split('_')[0])
    proximity_mutations = []

    # Split mutations into position and base, and convert positions to integers
    split_mutations = [(int(mutation.split('_')[0]), mutation) for mutation in mutations]

    for pos, mutation in split_mutations:
        # Check if the mutation is within 'proxi' positions of the specific mutation
        if abs(pos - specific_mutation_pos) <= proxi:
            if mutation != specific_mutation:
                proximity_mutations.append(mutation)

    return proximity_mutations

def derive_correct_length_headers(consensus_dict, fsa_file):
    correct_length_dict = {}
    for allele in consensus_dict:
        gene = allele.split('_')[0]
        correct_length_dict[gene] = [len(consensus_dict[allele][0]), []]
    sequence = ''
    with open(fsa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if gene != None:
                    if sequence != '':
                        if len(sequence) == correct_length_dict[gene][0]:
                            correct_length_dict[gene][1].append(sequence)
                header = line.strip()[1:]
                allele = header
                gene = allele.split('_')[0]
                sequence = ''
            else:
                sequence += line.strip()
    if gene != None:
        if sequence != '':
            if len(sequence) == correct_length_dict[gene][0]:
                correct_length_dict[gene][1].append(sequence)

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

def build_consensus_dict(res_file, mat_file):
    consensus_dict = {}
    with open(res_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                allele = line[0]
                length = int(line[3])
                consensus_dict[allele] = [[], '']
                for i in range(length):
                    consensus_dict[allele][0].append([0, 0, 0, 0, 0, 0]) #[A, C, G, T, N, -]

    with open(mat_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                allele = line[1:].strip()
                index = 0
            elif line != '':
                line = line.split('\t')
                if line[0] != '-': #excludes read gaps. Reconsider?
                    line = line[1:]
                    for i in range(len(line)):
                        consensus_dict[allele][0][index][i] += int(line[i])
                    index += 1

    for allele in consensus_dict:
        for position in consensus_dict[allele][0]:
            consensus_dict[allele][1] += 'ACGT'[position[:4].index(max(position[:4]))]
    return consensus_dict


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
