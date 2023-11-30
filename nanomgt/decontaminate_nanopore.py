import os
import sys
import gzip
from Bio import SeqIO
from nanomgt import kma
from nanomgt import util
from itertools import combinations

from nanomgt.nanopore_mutations import parse_sam_and_find_mutations
from nanomgt.nanopore_mutations import extract_alignment
from nanomgt.nanopore_mutations import create_mutation_vector
from nanomgt.nanopore_mutations import identify_mutations

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
    primary_specie = primary.split()[1] + ' ' + primary.split()[2]

    produce_specie_specific_kma_db(primary_specie,
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

    confirmed_mutation_dict, co_occurrence_tmp_dict = co_occurrence_until_convergence(arguments, confirmed_mutation_dict, consensus_dict, read_positions_blacklisted_dict, bio_validation_dict)
    format_output(confirmed_mutation_dict, consensus_dict, bio_validation_dict, co_occurrence_tmp_dict)

    sys.exit()

def co_occurrence_until_convergence(arguments, confirmed_mutation_dict, consensus_dict, read_positions_blacklisted_dict, bio_validation_dict):
    """
    Iteratively identifies co-occurring mutations in reads until convergence is reached.

    Args:
    arguments (object): An object containing various parameters for the function.
    confirmed_mutation_dict (dict): A dictionary of confirmed mutations.
    consensus_dict (dict): A dictionary containing consensus data.
    read_positions_blacklisted_dict (dict): A dictionary with read positions that are blacklisted.
    bio_validation_dict (dict): A dictionary used for biological validation.

    Returns:
    dict: The updated dictionary of confirmed mutations after convergence.
    """

    current_count = count_mutations_in_mutations_dict(confirmed_mutation_dict)
    iteration_count = 0

    # Iterate until no new mutations are found
    while True:
        confirmed_mutation_dict, co_occurrence_tmp_dict = upper_co_occuring_mutations_in_reads(
            arguments,
            confirmed_mutation_dict,
            consensus_dict,
            read_positions_blacklisted_dict,
            bio_validation_dict
        )

        new_count = count_mutations_in_mutations_dict(confirmed_mutation_dict)
        iteration_count += 1

        # Check for convergence
        if new_count == current_count:
            break
        current_count = new_count

    return confirmed_mutation_dict, co_occurrence_tmp_dict


def check_arguments(arguments):
    """
    Validates the given arguments to ensure they are within expected ranges.

    Args:
    arguments (object): An object containing various parameters to validate.

    This function exits the program if any argument is outside its valid range.
    """

    # Validate 'cor' argument
    if arguments.cor > 1:
        print('cor must be between 0 and 1, otherwise it is not a reward.')
        sys.exit()

    # Validate 'bp' argument
    if arguments.bp < 1:
        print('bp must be greater than 1, otherwise it is not a penalty.')
        sys.exit()

    # Validate 'pp' argument
    if arguments.pp < 1:
        print('pp must be greater than 1, otherwise it is not a penalty.')
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


def index_top_hits_db(output_directory):
    """
    Indexes the top hits from an initial RMLST alignment result and creates a FASTA file and database for these hits.

    Args:
    output_directory (str): The directory where the input files are located and where the output files will be saved.

    The function reads two files:
    - 'initial_rmlst_alignment.res' to identify the top hits.
    - 'specie_db.name' to match and index these hits.

    It then creates a FASTA file of the top hits and indexes it for use with the KMA tool.
    """
    # Read and process initial alignment results to identify top hits
    initial_alignment_file = os.path.join(output_directory, 'initial_rmlst_alignment.res')
    top_hits = {}
    with open(initial_alignment_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                elements = line.strip().split('\t')
                allele, score = elements[:2]
                gene = allele.split('_')[0]

                # Update the top hits for each gene
                if gene not in top_hits or score > top_hits[gene][1]:
                    top_hits[gene] = [allele, score, 0]

    # Match the top hits with their respective sequence number
    species_name_file = os.path.join(output_directory, 'specie_db.name')
    with open(species_name_file, 'r') as file:
        for sequence_number, line in enumerate(file, start=1):
            allele = line.strip()
            gene = allele.split('_')[0]

            if gene in top_hits and top_hits[gene][0] == allele:
                top_hits[gene][-1] = sequence_number

    # Prepare the sequences for FASTA conversion and indexing
    sequence_numbers = ','.join(str(top_hits[gene][-1]) for gene in top_hits)

    # Generate FASTA file from the sequences
    fasta_file = os.path.join(output_directory, 'top_hits.fsa')
    specie_db = os.path.join(output_directory, 'specie_db')
    os.system(f'kma seq2fasta -seqs {sequence_numbers} -t_db {specie_db} > {fasta_file}')

    # Index the generated FASTA file
    top_hits_db = os.path.join(output_directory, 'top_hits_db')
    os.system(f'kma index -i {fasta_file} -o {top_hits_db} 2>/dev/null')

def adjust_consensus_dict_for_individual_qscores(consensus_dict, sam_file, fastq_file):
    """
    Adjusts the consensus dictionary with individual quality scores from SAM and FASTQ files.

    Args:
    consensus_dict (dict): Dictionary containing consensus data.
    sam_file (str): Path to the SAM file.
    fastq_file (str): Path to the FASTQ file.

    Returns:
    tuple: A tuple containing the adjusted consensus dictionary and blacklisted positions.
    """
    # Blacklist certain positions based on the FASTQ file
    black_listed_positions = blacklist_positions(fastq_file, 2)

    # Initialize the adjusted consensus dictionary
    adjusted_consensus_dict = {}

    # Set up the consensus dictionary structure
    for allele, data in consensus_dict.items():
        allele_positions, allele_seq = data
        adjusted_consensus_dict[allele] = [[[0, 0, 0, 0, 0, 0] for _ in range(len(allele_positions))], allele_seq]

    total_black_list_count = 0

    # Process alignments from the SAM file
    with open(sam_file, 'r') as file:
        for alignment in file:
            # Skip header lines
            if alignment.startswith('@'):
                continue

            # Parse the alignment data
            cols = alignment.strip().split('\t')
            qname, _, rname, pos, _, cigar_str, _, _, _, seq = cols[:10]
            read_id = qname.split(' ')[0]
            template_seq = adjusted_consensus_dict[rname][1]
            pos = int(pos)
            tlen = len(template_seq)

            # Process only full-length reads
            if pos == 1 and len(seq) >= tlen:
                aligned_ref, aligned_query = extract_alignment(template_seq[pos - 1:pos - 1 + tlen], seq, cigar_str)
                mutation_vector = create_mutation_vector(aligned_ref, aligned_query)

                # Update the consensus dictionary
                for i, mutation in enumerate(mutation_vector):
                    nucleotide_list = ['A', 'C', 'G', 'T', 'N', '-']
                    if i not in black_listed_positions[read_id]:
                        adjusted_consensus_dict[rname][0][i][nucleotide_list.index(mutation)] += 1
                    else:
                        total_black_list_count += 1

    return adjusted_consensus_dict, black_listed_positions
def blacklist_positions(fastq_file, quality_threshold):
    """
    Create a blacklist of low-quality positions in a FASTQ file for each read.

    Args:
        fastq_file (str): The path to the FASTQ file.
        quality_threshold (int): The quality score threshold for blacklisting positions.

    Returns:
        dict: A dictionary where keys are read IDs and values are lists of blacklisted positions.
    """
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
def format_output(confirmed_mutation_dict, consensus_dict, bio_validation_dict, co_occurrence_tmp_dict):
    """
    Format and print the output of confirmed mutations with additional information.

    Args:
        confirmed_mutation_dict (dict): A dictionary containing confirmed mutations for alleles.
        consensus_dict (dict): A dictionary containing consensus information for alleles.
        bio_validation_dict (dict): A dictionary containing biological validation data for genes.

    Returns:
        None
    """
    header = 'Gene,Position,MajorityBase,MutationBase,MutationDepth,TotalDepth,GeneLength,MutationComment,CoOccurrence'
    print(header)

    for allele in confirmed_mutation_dict:
        for mutation in zip(confirmed_mutation_dict[allele][0], confirmed_mutation_dict[allele][1]):
            position = mutation[0].split('_')[0]
            mutation_base = mutation[0].split('_')[1]
            mutation_depth = mutation[1]
            majority_base = consensus_dict[allele][1][int(position) - 1]
            total_depth = sum(consensus_dict[allele][0][int(position) - 1])
            biological_existence = check_single_mutation_existence(bio_validation_dict, allele, mutation[0])
            gene_length = len(consensus_dict[allele][1])
            if mutation[0] in co_occurrence_tmp_dict[allele]:
                co_occurrence = 'Yes'
            else:
                co_occurrence = 'No'

            if biological_existence:
                print('{},{},{},{},{},{},{},{},{}'.format(allele, position, majority_base, mutation_base, mutation_depth, total_depth, gene_length, 'Mutation seen in database', co_occurrence))
            else:
                print('{},{},{},{},{},{},{},{},{}'.format(allele, position, majority_base, mutation_base, mutation_depth, total_depth, gene_length, 'Novel mutation', co_occurrence))

def extract_mapped_rmlst_read(output_directory, nanopore_fastq):
    """
    Extract and trim mapped rMLST reads from an initial alignment file.

    Args:
        output_directory (str): The directory where output files will be saved.
        nanopore_fastq (str): The path to the nanopore FASTQ file.

    Returns:
        None
    """
    read_set = set()

    # Extract read IDs from the initial rMLST alignment file
    with open(output_directory + '/initial_rmlst_alignment.frag', 'r') as frag:
        for line in frag:
            line = line.rstrip()
            line = line.split('\t')
            read_set.add(line[-1])

    # Write the extracted read IDs to a text file
    with open(output_directory + '/rmlst_reads.txt', 'w') as f:
        for item in read_set:
            f.write(item + '\n')

    # Use seqtk to extract the corresponding reads from the nanopore FASTQ file
    os.system('seqtk subseq {} {} > {}'.format(nanopore_fastq, output_directory + '/rmlst_reads.txt',
                                               output_directory + '/trimmed_rmlst_reads.fastq'))



def derive_mutation_positions(consensus_dict, arguments):
    """
    Derive mutation positions and their depths from a consensus dictionary.

    Args:
        consensus_dict (dict): A dictionary containing consensus information for alleles.
        arguments: Arguments containing parameters for mutation derivation.

    Returns:
        dict: A dictionary containing derived mutation positions and depths for each allele.
    """
    all_confirmed_mutation_dict = {}

    for allele, allele_data in consensus_dict.items():
        all_confirmed_mutation_dict[allele] = [[], []]

        for i in range(len(allele_data[0])):
            positions = allele_data[0][i][:4]
            max_number = max(positions)
            index_of_max = positions.index(max_number)
            nucleotide_index = ['A', 'C', 'G', 'T']

            for t in range(len(positions)):
                if t != index_of_max:
                    if positions[t] >= arguments.min_n:
                        total_depth = sum(positions)
                        relative_depth = positions[t] / total_depth

                        if relative_depth >= arguments.cor * arguments.mrd:
                            # Only consider mutations with minimum depth >= 2
                            all_confirmed_mutation_dict[allele][0].append(
                                '{}_{}'.format(i + 1, nucleotide_index[t]))
                            all_confirmed_mutation_dict[allele][1].append(positions[t])

    return all_confirmed_mutation_dict






def upper_co_occuring_mutations_in_reads(arguments, confirmed_mutation_dict, consensus_dict,
                                         read_positions_blacklisted_dict, bio_validation_dict):
    """
    Filter and adjust confirmed mutations based on co-occurrence, depth, and biological validation.

    Args:
        arguments: Arguments containing parameters for filtering.
        confirmed_mutation_dict (dict): A dictionary containing confirmed mutations for alleles.
        consensus_dict (dict): A dictionary containing consensus information for alleles.
        read_positions_blacklisted_dict (dict): A dictionary of blacklisted positions.
        bio_validation_dict (dict): A dictionary containing biological validation data for genes.

    Returns:
        dict: A filtered and adjusted mutation dictionary for alleles.
    """
    reads_mutation_dict = parse_sam_and_find_mutations(arguments.output + '/rmlst_alignment.sam',
                                                       confirmed_mutation_dict,
                                                       consensus_dict,
                                                       read_positions_blacklisted_dict)

    co_occurrence_tmp_dict = {}
    co_occurrence_matrix_dict = {}
    for allele in confirmed_mutation_dict:
        mutation_list = confirmed_mutation_dict[allele][0]
        num_mutations = len(mutation_list)
        if num_mutations > 1:
            co_occurrence_matrix = [[0] * num_mutations for _ in range(num_mutations)]
            for read in reads_mutation_dict:
                read_allele = read.split(' ')[1]
                if read_allele == allele:
                    read_mutations = reads_mutation_dict[read]
                    valid_mutations = [mutation for mutation in read_mutations if mutation in mutation_list]
                    if len(valid_mutations) > 1:
                        for i in range(len(valid_mutations)):
                            for j in range(i + 1, len(valid_mutations)):
                                mutation1 = mutation_list.index(valid_mutations[i])
                                mutation2 = mutation_list.index(valid_mutations[j])
                                co_occurrence_matrix[mutation1][mutation2] += 1
                                co_occurrence_matrix[mutation2][mutation1] += 1
            co_occurrence_matrix_dict[allele] = [co_occurrence_matrix, mutation_list]

    adjusted_mutation_dict = {}
    for allele in confirmed_mutation_dict:
        co_occurrence_tmp_dict[allele] = []
        if allele in co_occurrence_matrix_dict:
            adjusted_mutation_dict[allele] = [[], []]
            matrix = co_occurrence_matrix_dict[allele][0]
            mutation_list = co_occurrence_matrix_dict[allele][1]
            depth_list = confirmed_mutation_dict[allele][1]
            for i in range(len(matrix)):
                row = matrix[i]
                mutation = mutation_list[i]
                position = int(mutation.split('_')[0])
                position_depth = sum(consensus_dict[allele][0][position - 1])
                mutation_depth = depth_list[i]
                proxi_mutations = find_mutations_proximity_specific_mutation(mutation_list, mutation, arguments.proxi)
                biological_existence = check_single_mutation_existence(bio_validation_dict, allele, mutation)

                mutation_threshold = position_depth * arguments.mrd
                co_occurrence_list = check_mutation_co_occurrence(row, mutation_list, mutation,
                                                                 position_depth, arguments.cor, arguments.pp, arguments.mrd, proxi_mutations)

                if co_occurrence_list != []:
                    co_occurrence_tmp_dict[allele].append(co_occurrence_list)
                    mutation_threshold = mutation_threshold * arguments.cor
                if not biological_existence:
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
                biological_existence = check_single_mutation_existence(bio_validation_dict, allele, mutation)
                if not biological_existence:
                    mutation_threshold = mutation_threshold + (arguments.bp-1) * position_depth * arguments.mrd

                if depth >= mutation_threshold:
                    adjusted_mutation_dict[allele][0].append(confirmed_mutation_dict[allele][0][0])
                    adjusted_mutation_dict[allele][1].append(confirmed_mutation_dict[allele][1][0])
    return adjusted_mutation_dict, co_occurrence_tmp_dict



def check_mutation_co_occurrence(list_of_mutation_co_occurrence, mutation_list, mutation,
                                 position_depth, correlation_coefficient, proximity_penalty, mrd, proximity_mutations):
    """
    Check for co-occurrence of a mutation with other mutations in a list.

    Args:
        list_of_mutation_co_occurrence (list): A list of co-occurrence counts for each mutation.
        mutation_list (list): A list of mutations.
        mutation (str): The mutation for which co-occurrence is being checked.
        position_depth (int): The depth at which the mutation occurs.
        correlation_coefficient (float): The correlation coefficient used for threshold calculation.
        proximity_penalty (float): The penalty factor for mutations within proximity.
        mrd (float): The mutation rate difference.
        proximity_mutations (list): A list of mutations within proximity.

    Returns:
        list: A list of mutations that co-occur with the given mutation.
    """
    if mutation not in mutation_list:
        # Should never happen
        return []  # No co-occurrence and not in proximity

    co_threshold = position_depth * mrd * correlation_coefficient  # Default co_threshold
    if co_threshold < 3:
        co_threshold = 3

    # Find the index of the mutation in the mutation list
    mutation_index = mutation_list.index(mutation)

    co_occurrence_list = []
    # Check if the co-occurrence count of the mutation with any other mutation is above the threshold
    for i, count in enumerate(list_of_mutation_co_occurrence):
        if mutation_list[i] in proximity_mutations:
            # Add penalty for proximity to make it harder to get the co-occurrence reward
            # for mutations within the proximity
            co_threshold = co_threshold * proximity_penalty
        if i != mutation_index and count >= co_threshold:
            co_occurrence_list.append(mutation_list[i])

    return co_occurrence_list

def check_single_mutation_existence(bio_validation_dict, allele, specific_mutation):
    """
    Check if a specific mutation exists for a given allele and gene in a biological validation dictionary.

    Args:
        bio_validation_dict (dict): A dictionary containing biological validation data for genes.
        allele (str): The allele for which the existence of a specific mutation is checked.
        specific_mutation (str): The specific mutation to check for.

    Returns:
        bool: True if the specific mutation exists for the allele, False otherwise.
    """
    gene = allele.split('_')[0]

    if specific_mutation in bio_validation_dict.get(gene, []):
        return True

    return False

def find_mutations_proximity_specific_mutation(mutations, specific_mutation, proximity):
    """
    Find mutations that are in proximity (within a certain number of positions) to a specific mutation.

    Args:
        mutations (list): A list of mutation strings in the format "position_base".
        specific_mutation (str): The specific mutation to which proximity is determined.
        proximity (int): The maximum number of positions for mutations to be considered in proximity.

    Returns:
        list: A list of mutations that are in proximity to the specific mutation.
    """
    specific_mutation_pos = int(specific_mutation.split('_')[0])
    proximity_mutations = []

    # Split mutations into position and base, and convert positions to integers
    split_mutations = [(int(mutation.split('_')[0]), mutation) for mutation in mutations]

    for pos, mutation in split_mutations:
        # Check if the mutation is within 'proximity' positions of the specific mutation
        if abs(pos - specific_mutation_pos) <= proximity:
            if mutation != specific_mutation:
                proximity_mutations.append(mutation)

    return proximity_mutations

def derive_correct_length_headers(consensus_dict, fsa_file):
    """
    Derive correct length headers and sequences from a FASTA file based on a consensus dictionary.

    Args:
        consensus_dict (dict): A dictionary containing consensus information for alleles.
        fsa_file (str): The path to the input FASTA file.

    Returns:
        dict: A dictionary mapping gene names to correct length sequences.
    """
    correct_length_dict = {}

    for allele in consensus_dict:
        gene = allele.split('_')[0]
        correct_length_dict[gene] = [len(consensus_dict[allele][0]), []]

    sequence = ''
    gene = None

    with open(fsa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if gene is not None:
                    if sequence != '':
                        if len(sequence) == correct_length_dict[gene][0]:
                            correct_length_dict[gene][1].append(sequence)
                header = line.strip()[1:]
                allele = header
                gene = allele.split('_')[0]
                sequence = ''
            else:
                sequence += line.strip()

    if gene is not None:
        if sequence != '':
            if len(sequence) == correct_length_dict[gene][0]:
                correct_length_dict[gene][1].append(sequence)

    return correct_length_dict

def produce_specie_specific_kma_db(specie, fsa_file, scheme_file, output_directory):
    """
    Produce a specie-specific KMA database using the provided FASTA file and scheme file.

    Args:
        specie (str): The target specie for which the KMA database is being created.
        fsa_file (str): The path to the input FASTA file.
        scheme_file (str): The path to the scheme file containing gene information.
        output_directory (str): The directory where the output KMA database will be created.

    Returns:
        None
    """
    gene_set = set()
    t = 0

    with open(scheme_file, 'r') as f:
        for line in f:
            if line.startswith('rST'):
                headers = line.strip().split('\t')[1:54]
            else:
                if line.strip().split('\t') != ['']:
                    if line.strip().split('\t')[55] == specie:
                        t += 1
                        for i in range(len(headers)):
                            allele = headers[i] + '_' + line.strip().split('\t')[i + 1]
                            gene_set.add(allele)

    # Create a specie-specific FASTA file with the selected genes
    produce_specie_fsa_file(fsa_file, gene_set, output_directory)

    # Create a KMA database from the specie-specific FASTA file
    os.system('kma index -i {}/specie.fsa -o {}/specie_db 2>/dev/null'.format(output_directory, output_directory))


def produce_specie_fsa_file(fsa_file, gene_set, output_directory):
    """
    Produce a specie-specific FASTA file containing sequences for genes in the given gene set.

    Args:
        fsa_file (str): The path to the input FASTA file.
        gene_set (set): A set containing gene IDs to include in the output.
        output_directory (str): The directory where the output file will be saved.

    Returns:
        None
    """
    output_file = output_directory + '/specie.fsa'

    with open(output_file, 'w') as outfile:
        with open(fsa_file, 'r') as f:
            write_sequence = False  # A flag to indicate whether the sequence should be written to the output
            for line in f:
                if line.startswith('>'):
                    # Check if the gene_id (without '>') is in the gene_set
                    gene_id = line.strip().split()[0][1:]
                    write_sequence = gene_id in gene_set
                # Write the line (header or sequence) if write_sequence is True
                if write_sequence:
                    outfile.write(line)

#Here
def build_consensus_dict(res_file, mat_file):
    """
    Build a consensus dictionary from result and matrix files.

    Args:
        res_file (str): The name of the result file.
        mat_file (str): The name of the matrix file.

    Returns:
        dict: A dictionary containing consensus information for alleles.
    """
    consensus_dict = {}

    # Read and process the result file
    with open(res_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                allele = line[0]
                length = int(line[3])
                consensus_dict[allele] = [[], '']
                for i in range(length):
                    consensus_dict[allele][0].append([0, 0, 0, 0, 0, 0])  # [A, C, G, T, N, -]

    # Read and process the matrix file
    with open(mat_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                allele = line[1:].strip()
                index = 0
            elif line != '':
                line = line.split('\t')
                if line[0] != '-':  # Excludes read gaps. Reconsider?
                    line = line[1:]
                    for i in range(len(line)):
                        consensus_dict[allele][0][index][i] += int(line[i])
                    index += 1

    # Generate consensus sequences for alleles
    for allele in consensus_dict:
        for position in consensus_dict[allele][0]:
            consensus_dict[allele][1] += 'ACGT'[position[:4].index(max(position[:4]))]

    return consensus_dict



def number_of_bases_in_file(filename):
    """
    Calculate the total number of bases in a FASTA or FASTQ file.

    Args:
        filename (str): The name of the input file.

    Returns:
        int: The total number of bases in the file.
    """
    gzipped, file_type = determine_file_type(filename)

    if file_type == 'fasta':
        total_bases = 0
        with open(filename, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    total_bases += len(line.strip())
        return total_bases

    elif file_type == 'fastq':
        total_bases = 0
        line_count = 1
        with (gzip.open(filename, 'r') if gzipped else open(filename, 'r')) as f:
            for line in f:
                if line_count == 2:
                    total_bases += len(line.strip())
                line_count += 1
                if line_count == 5:
                    line_count = 1
        return total_bases

def eval_bacteria_results(results, total_bacteria_aligning_bases):
    """
    Evaluate bacterial alignment results to determine primary and candidate results.

    Args:
        results (list of dict): List of dictionaries containing bacterial alignment results.
        total_bacteria_aligning_bases (int): Total number of bases aligning to bacteria.

    Returns:
        tuple: A tuple containing the primary result and a dictionary of candidate results.
    """
    primary = results[0]['#Template']
    candidate_dict = dict()

    for item in results:
        # Calculate the number of bases hit based on depth
        bases_hit = int(item['Template_length']) * float(item['Depth'].strip())

        # Calculate the template identity
        template_id = float(item['Template_Identity'])

        # Calculate relative template depth
        relative_template_depth = bases_hit / total_bacteria_aligning_bases

        # Check if the result qualifies as a candidate
        if relative_template_depth > 0.01 or template_id > 20.00:
            candidate_dict[item['#Template']] = [relative_template_depth, int(item['Template_length'])]

    return primary, candidate_dict


def drive_bacteria_results(arguments, total_bacteria_aligning_bases):
    """
    Drive the analysis of bacterial alignment results.

    Args:
        arguments: Parsed command-line arguments.
        total_bacteria_aligning_bases (int): Total number of bases aligning to bacteria.

    Returns:
        tuple: A tuple containing primary results and a dictionary of candidate results.
    """
    # Sort bacterial alignment results by score
    results = sort_lines_by_score(arguments.output + "/bacteria_alignment.res")

    # Evaluate the bacterial results to determine primary and candidates
    primary, candidate_dict = eval_bacteria_results(results, total_bacteria_aligning_bases)

    return primary, candidate_dict


def determine_file_type(file):
    """
    Determine the file type and whether it is gzipped.

    Args:
        file (str): The name of the input file.

    Returns:
        tuple: A tuple containing a boolean indicating if the file is gzipped and a string indicating the file type.
    """
    gzipped = False
    file_type = None

    # Check if the file has a '.gz' extension, indicating it is gzipped
    if file.endswith('.gz'):
        gzipped = True
        file = file[:-3]  # Remove the '.gz' extension

    # Check the file type based on its extension
    if file.endswith('.fastq') or file.endswith('.fq'):
        file_type = 'fastq'
    elif file.endswith('.fasta') or file.endswith('.fa') or file.endswith('.fna') or file.endswith('.fsa'):
        file_type = 'fasta'

    return gzipped, file_type


def sort_lines_by_score(filename):
    """
    Sort lines in a tab-delimited file by the 'Score' column in descending order.

    Args:
        filename (str): The name of the input file.

    Returns:
        list: A list of dictionaries representing the sorted data.
    """
    data = []

    # Read data from the input file
    with open(filename, 'r') as file:
        lines = file.readlines()
        header = lines[0].strip().split('\t')

        # Parse and store data as dictionaries
        for line in lines[1:]:
            values = line.strip().split('\t')
            data.append(dict(zip(header, values)))

    # Sort the data by the 'Score' column in descending order
    sorted_data = sorted(data, key=lambda x: int(x['Score']), reverse=True)

    return sorted_data
