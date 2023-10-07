import re


def extract_alignment(reference, query, cigar_string):
    """
    Extracts the alignment between reference and query based on the provided CIGAR string.

    Parameters:
    - reference (str): The reference sequence.
    - query (str): The query sequence.
    - cigar_string (str): The CIGAR string.

    Returns:
    - str, str: The portions of the reference and query sequences that align, respectively.
    """
    # Parse the CIGAR string
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)

    aligned_ref = ""
    aligned_query = ""

    ref_index = 0
    query_index = 0

    # Loop through CIGAR tuples and build the aligned sequences
    for length, operation in cigar_tuples:
        length = int(length)

        if operation == "M":  # Match or mismatch
            aligned_ref += reference[ref_index:ref_index + length]
            aligned_query += query[query_index:query_index + length]
            ref_index += length
            query_index += length
        elif operation == "I":  # Insertion
            aligned_ref += "-" * length
            aligned_query += query[query_index:query_index + length]
            query_index += length
        elif operation == "D":  # Deletion
            aligned_ref += reference[ref_index:ref_index + length]
            aligned_query += "-" * length
            ref_index += length
        elif operation == "S":  # Soft clipping
            query_index += length
        elif operation == "=":  # Match
            aligned_ref += reference[ref_index:ref_index + length]
            aligned_query += query[query_index:query_index + length]
            ref_index += length
            query_index += length
        elif operation == "X":  # Mismatch
            aligned_ref += reference[ref_index:ref_index + length]
            aligned_query += query[query_index:query_index + length]
            ref_index += length
            query_index += length
        # Optionally handle other CIGAR operations (N, H, P) as per requirements

    return aligned_ref, aligned_query

def create_mutation_vector(aligned_ref, aligned_query):
    """
    Creates a mutation vector based on the aligned reference and query sequences.

    Parameters:
    - aligned_ref (str): Aligned portion of the reference sequence.
    - aligned_query (str): Aligned portion of the query sequence.

    Returns:
    - list[str]: A mutation vector representing the alignment with respect to the original reference sequence.
    """
    mutation_vector = []

    # Loop through the aligned sequences
    for ref_nt, query_nt in zip(aligned_ref, aligned_query):
        if ref_nt == "-":
            # Skip insertions in the reference
            continue
        elif query_nt == "-":
            # Represent deletions in the query with "-"
            mutation_vector.append("-")
        else:
            # Otherwise, represent the nucleotide from the query
            mutation_vector.append(query_nt)

    return mutation_vector


def identify_mutations(mutation_vector, reference_sequence):
    """
    Identify all mutation positions from a mutation vector compared to the reference.

    Parameters:
    - mutation_vector (list[str]): The mutation vector.
    - reference_sequence (str): The original reference sequence.

    Returns:
    - list[str]: A list where each mutation is described as a string in the format "POSITION_NUCLEOTIDE".
    """
    mutations = []

    # Ensure that mutation vector and reference have equal lengths
    if len(mutation_vector) != len(reference_sequence):
        raise ValueError("The mutation vector and reference sequence must have the same length.")

    # Loop through the mutation vector and reference sequence
    for i, (mv_nt, ref_nt) in enumerate(zip(mutation_vector, reference_sequence)):
        # If the nucleotide in the mutation vector is different from the reference
        # and is not a "-" (representing a deletion)
        if mv_nt != ref_nt and mv_nt != "-":
            # Add to mutations in the format "POSITION_NUCLEOTIDE"
            mutations.append(f"{i + 1}_{mv_nt}")

    return mutations

def parse_sam_and_find_mutations(sam_filepath):
    """
    Parses a SAM file, extracts necessary information and finds mutations in each read.

    Parameters:
    - sam_filepath (str): The path to the SAM file.

    Returns:
    - dict: A dictionary where keys are read names and values are lists of mutation strings.
    """
    mutations_dict = {}

    with open(sam_filepath, 'r') as sam_file:
        for line in sam_file:
            # Skip header lines
            if line.startswith('@'):
                continue

            # Extract relevant columns: [QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ]
            cols = line.strip().split('\t')
            qname, flag, rname, pos, mapq, cigar_str, rnext, pnext, tlen, seq = cols[:10]

            # Convert string columns to appropriate types
            pos = int(pos)
            tlen = int(tlen)

            # Assuming reference is provided or it is known in another way
            reference = "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG"

            # Obtaining the alignment using your function
            aligned_ref, aligned_query = align_sequences(reference, seq, cigar_str)

            # Creating a mutation vector using your function
            mutation_vector = get_mutation_vector(aligned_ref, aligned_query)

            # Identifying mutations using your function
            mutations = identify_mutations(mutation_vector, reference[pos - 1:pos - 1 + tlen])

            # Storing mutations in the dictionary
            mutations_dict[qname] = mutations

    return mutations_dict

def parse_sam_get_references(sam_file_path):
    """
    Parse the SAM file and extract all unique reference names.

    :param sam_file_path: Path to the SAM file.
    :return: Set of reference names.
    """
    references = set()
    with open(sam_file_path, 'r') as sam_file:
        for line in sam_file:
            # Skip header lines
            if line.startswith("@"):
                continue

            fields = line.split("\t")
            if len(fields) > 2:
                references.add(fields[2])  # The reference name is in the 3rd field
    return references


def load_references_from_fasta(fasta_file_path, references):
    """
    Load reference sequences from a FASTA file that are in the provided reference names set.

    :param fasta_file_path: Path to the FASTA file.
    :param references: Set of reference names to load.
    :return: Dictionary with reference names as keys and sequences as values.
    """
    ref_sequences = {}
    with open(fasta_file_path, 'r') as fasta_file:
        current_ref_name = None
        current_sequence = []
        for line in fasta_file:
            line = line.strip()
            # Check if the line is a header line
            if line.startswith(">"):
                # If we were building a sequence, add it to the dictionary
                if current_ref_name and (current_ref_name in references):
                    ref_sequences[current_ref_name] = "".join(current_sequence)

                # Get the new reference name
                current_ref_name = re.sub(r"^>\s*", "", line).split()[0]  # Remove ">" and extract the first word

                # Reset the sequence
                current_sequence = []
            else:
                # Add the line to the sequence
                current_sequence.append(line)

        # Add the last sequence
        if current_ref_name and (current_ref_name in references):
            ref_sequences[current_ref_name] = "".join(current_sequence)
    return ref_sequences


# Example usage:
sam_file_path = "test/test.sam"
fasta_file_path = "specie.fsa"

references = parse_sam_get_references(sam_file_path)
ref_sequences = load_references_from_fasta(fasta_file_path, references)

for ref_name, ref_sequence in ref_sequences.items():
    print(ref_name, ref_sequence)


# Example usage:
# sam_filepath = "path_to_your_file.sam"
# mutations_dict = parse_sam_and_find_mutations(sam_filepath)
# print(mutations_dict)


# Example usage
"""
reference = "ATGAAAGTAAGACCATCAGTAAAACCTATTTGCGAAAAATGTAAAATCATTAAACGTAAAGGTAAAGTAATGGTAATTTGTGAAAATCCAAAACACAAACAAAGACAAGGTTAA"
query = "ACTTTTGGTGGCGGAGCCTTGAAACGATTTTCTATTTTTGCAATGGTATTGTACCCTACATCACTGCATCAATCGTAATGCAATTATTACAAATGGATATTGTCCCTAAATTCTCGAGAATGGGCAAAACAAGGTGAAGTAGGTAGAAGAAAGTTAAATAACGTTACTCGTTATTTAGCAATTTCTTTAGCATTTATCCAATCTATAGGTATGGCATTCCAATTTAATAATTATCTCAAAGGTGCTGATTATCAATCAGTCAATTATGAGTTATTTATTAATAGCACTAGTTTTCACAGCAGGAACTGCTTTCTTAATATGGCTTGGTGATCAAATCACTCAGTTCGGTGTTGGTGGTAATGGTATTTCTATTATCATATTTGCGGGTATTTTATCAACGTTGCCAGCCAGTCTAATTGAATTGGACAAACTGCGTTCGTTGGGCAAAGAAGATACATCATTAGCATGGTTGAAAGTATTAGGACTTTTAGTGTCACTGATACTATTAACAGTAGGTGCTATTTATGTACTTGAAGCCGTTCGTAAAATTCCGATTCAGATGAAAGAACAAACTGCTCAAAGATTAGGTTCACAGGCAACTTATCTACCTTTAAAAGTGAACTCAGCCGGTGTTATTCCTGTAATCTTTGCAATGGCGTTCTTCTTACTGCCTAGAACATTAACGTTGTTCTATCCAGATAAAGAATGGGCACAAAACATTGCGAATGCTGCAGATCCTTCAAGTAATGTTGGTATGGTTGCAGACATCGTATTAATTATACTATTCACATACTTCTATGCATTCGTAGAGCTTAATCCAAAAAATGGCTGATAACCTTAAGAAACAAGGTAGCTATGTTCCAGGTATTAGACCTGGTGAACAAACCAAAATATATTACTAAAGTTTTATATCAATCAAGTCTTCGCAGGTTTATCTCCTTAGCCCGTTATTTCAATTCTACCGACTAGCAATGGAAAAATTCTTATGTGGGATTACCACAATCAATTCAGATTGGTGGTACAAGTTTACTTATCGTTATTGGTGTAGCGATTGAAACAATGAAATCTTTAGAAGCGCAAGATTCTCAAAAAAGAATATAAAGGCTTTGTGGTAGATAATTTTAGGAGGGCACTTATGAATATCATTTTGATGGGTTTACCTGGCGCAGGTAAAGGAACTCAAGCAAGTGAAATTGTCAAGAAATTCCCATATGCCATTTCAACTGGTGACATGTTCAGAAAAGCTATAAAAGAAGAAACTGAATTAGGTAAAGAAGCTAAGTCATGATGGACCGTGGCGAATTAGTTCCTGATGAAGTGACTGTAGGTATCGTTAAGGAAAGAATTTCTGAAGACGATGCAAAAAAGGCTTTTATTAGATGGCTTCCCAAGAAACAATCGAGCAAGCTGAGGCATTAAATAATATTATGTCTGAGCTTGACAGAAACATTGATGCTGTCATCAGATCGAAGTTCCGGAAGAAGAATTAATGAACCGTCTTACAGGTCGTCGAATCTGTGAGTCATGTGGTACAACGTATCATCTTGTATTTAATCCTCCGAAGTCGAAGGTATTTGTGATATCGATGGTGGTAAATTGTATCAACGAGAAACACGATATACCAAAACGGTAGCTAATCGTTTAGTGTTACATTAAACAATCTAAACCTATTTTAGATTTCTATGATCAAAAAGGTGTATTCTAAAATATTGATGGTTCAAAAGATATTAGTGATGTTACCAAAGATGTCATTGATATTTTAGATCATTTGTAATAGATCAACTTTTAGCGGTCTATTTTAAGTGAGGCATGATTTAATCATATTAGATAGATGAGAGTGGCAATTAACGATAATTAATTCCCGTAATTTTTTGTTTAACTTATAGAAAAATCTAATCACTGATAGATTCGTAACAATAAAGTTGTCTATTGACGATTACCTTACTATTGCAAAAAGGGGAAGTTAATCAATGGCTAAACAAGATGTAATTGAATTAGAAGGTACTGTATTAGATACTTTACCGAACGCAATGTTTAAAGTAGAATTAGAAAATGGTCATGAGATTTTAGCTCACGTAAGTGGTAAAATCAGAATGAATTACATTCGTATTCTACCTGGCGACAAAGTAACTGTTGAAGATGTCTCCGTACGATTTAACACGCGGAAGAATTACTTATCGTTATAAATAATCGTCACTCCATAATATAGGGAGGTATAAAAATGAAAGTAAGACCATCAGTAAAACCATATTTGCTGAAAAATGTAAAATCATTAAACGTAAAAAGGCAAAGAATGGTAATTTGTGAAAATCCAAAACACAAACAAAGACAAGGTTAATAAAAGAGAGGTGTAAATTAATATGGCACGTATTGCAGGAGTAGATATTCCACGTGAAAAACGCGTAGTTATCTCATTAACTTATATATACGGTATCGGTACGTCAACTGCTCAAAAAGTCTCAAAAGAAAGCTAACGTATCAAGTGACAACAGTCGTGTGAAGATTTAACTGATGACGAATCAAGGTCGTATCCCTGTGAAGTTGTAGACGGTTATAAAGTCGAAGGTGACTTACGTCGTGAAACTAACTTAAATATCAAACGTTTAATGGAAAGTTCATCATACCGTGGTATCCGTCACCGTCGTGGTTTACCAGTTCGTGGTCAAAAAACGAAAAACAACGCGCGTACTCGTAAAATACCAGTTAAAACGGTAGCTAACAAGAAAAAATAATAGGTAAAGGAGGCAAATTTTAAATGGCACGTAAAACAAGAATCGCAAAAAGGAGCGTAAAAAAGAATATTGAAAATGGTGTAGCACACATCCGTTCAACATTCAACAACACTATTGTAACTATCACTGATGAGTTCGGTAATGCTTTATCATGGTCATCAGCTGGTGCATTAGGATTCAAAGAAGGATCAAAAATCAACACCATTTGCAGCACAAATGGCTTCTGAAACTGCATCTAAATCAGCTATGGAGCATGGTTTAAAAACAGTTGAAGTAACAGTTAAATAATCTGTCAGGCTGTGGACAATCAGCTATTCGTGCATTACAATCTGCAGGTTTAGAAGTAACTGCGATCAGAGACGTTACTCCAGTACCTCATAACGGTTGTCGTCCACCAAAACGTCGTCGTGTATAATTTATGATGGTATTGTTACAGGTCACTGAGCAAACATAGCAATACGTAAC"
cigar_string = "2221S26=1I7=1I27=2I2=1X4=1D46=859S"

aligned_ref, aligned_query = extract_alignment(reference, query, cigar_string)
print("Aligned Reference: ", aligned_ref)
print("Aligned Query:     ", aligned_query)
mutation_vector = create_mutation_vector(aligned_ref, aligned_query)
print( "reference:         ", reference)
print("Mutation Vector:   ", "".join(mutation_vector))
mutations = identify_mutations(mutation_vector, reference)
print("Mutations: ", mutations)"""