import os
import sys
import gzip
from nanodecon import kma
from nanodecon import util

from nanodecon.intra_species_detection import determine_intra_species_contamination_nanopore

def nanopore_decontamination(arguments):
    os.system('mkdir ' + arguments.output)
    #os.system('mkdir ' + arguments.output + '/output/contaminations')
    #os.system('mkdir ' + arguments.output + '/alignments')
    #kma.KMARunner(arguments.nanopore,
    #              arguments.output + "/bacteria_alignment",
    #              arguments.db_dir + "/bac_db",
    #              "-mem_mode -1t1 -t {} -ID 10 -ont".format(arguments.threads)).run()

    total_bacteria_aligning_bases = util.number_of_bases_in_file(arguments.output + "/bacteria_alignment.fsa")
    primary, candidate_dict = drive_bacteria_results(arguments, total_bacteria_aligning_bases)
    primary_species = primary.split()[1] + ' ' + primary.split()[2]
    print(primary_species)
    print(primary_species)

    #produce_species_specific_kma_db(primary_species,
    #                                '/home/people/malhal/contamErase_db/rmlst.fsa',
    #                                '/home/people/malhal/contamErase_db/rmlst_scheme.txt',
    #                                arguments.output)
    #kma.KMARunner(arguments.nanopore,
    #              arguments.output + "/rmlst_alignment",
    #              arguments.output + '/specie_db',
    #              "-t {} -ID 10 -ont -md 1.5 -matrix -mp 14".format(arguments.threads)).run()

    #os.system('gunzip ' + arguments.output + '/rmlst_alignment.mat.gz')

    odd_size_alleles, non_alignment_matches, consensus_dict = build_consensus_dict(arguments,
                                                                                   arguments.output + '/rmlst_alignment.res',
                                                                                   arguments.output + '/rmlst_alignment.mat')

    confirmed_alleles = check_all_species_alleles_against_consensus_dict(consensus_dict,
                                                                         arguments.output + '/specie.fsa')
    for item in confirmed_alleles:
        if confirmed_alleles[item][0] >= 10:
            print(item, confirmed_alleles[item])
    sys.exit()
    calculate_rmlst_scheme_matches(confirmed_alleles, arguments.db_dir + '/rmlst_scheme.txt')

    sys.exit()
    sys.exit()
    produce_final_output_nanopore(arguments, arguments.output + '/bacteria_alignment.frag', primary, candidate_rmlst_dict_results, black_list_plasmid, black_list_viral, black_list_human)
    #produce_contamination_report #TBD

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

def check_all_species_alleles_against_consensus_dict(consensus_dict, fsa_file):
    confirmed_alleles = {}
    relative_threshold = 0.01
    with open(fsa_file, 'r') as f:
        sequence = ''
        min_depth = 100000
        for line in f:
            if line.startswith('>'):
                if sequence != '':
                    if len(sequence) == len(consensus_dict[allele]):
                        mutation_list = []
                        for i in range(len(sequence)):
                            #total_base_count = sum(consensus_dict[allele][i][:4])
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
                            if consensus_dict[allele][i][index] < min_depth:
                                min_depth = consensus_dict[allele][i][index]
                            major_nucleotide = max(consensus_dict[allele][i][:4])
                            if consensus_dict[allele][i][index] < major_nucleotide:
                                mutation_list.append('{}_{}'.format(consensus_dict[allele][i][index], sequence[i]))
                        if min_depth > 0 and min_depth != 100000:
                            confirmed_alleles[gene] = [min_depth, '', mutation_list]
                    sequence = ''
                    min_depth = 100000
                gene = line.strip()[1:]
                allele = gene.split('_')[0]
            else:
                sequence += line.strip()
        if len(sequence) == len(consensus_dict[allele]):
            if min_depth > 0 and min_depth != 100000:
                confirmed_alleles[gene] = [min_depth, '', mutation_list]
    return confirmed_alleles


def build_consensus_dict(arguments, res_file, mat_file):
    top_allele_dict = {}
    non_alignment_matches = {}
    consensus_dict = {}
    with open(res_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                allele = line[0].split('_')[0]
                if allele not in top_allele_dict:
                    top_allele_dict[allele] = [int(line[3]), int(line[1]), line[0]]
                else:
                    if int(line[1]) > top_allele_dict[allele][1]:
                        top_allele_dict[allele] = [int(line[3]), int(line[1]), line[0]]
    odd_size_alleles = set()
    correct_size_alleles = set()
    with open(res_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                allele = line[0].split('_')[0]
                if int(line[3]) != top_allele_dict[allele][0]:
                    odd_size_alleles.add(line[0])
                else:
                    correct_size_alleles.add(line[0])
                    if allele not in consensus_dict:
                        consensus_dict[allele] = []
                        for i in range(int(line[3])):
                            consensus_dict[allele].append([0, 0, 0, 0, 0, 0]) #[A, C, G, T, N, -]

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
    #Consider how we handle gaps in reads and template
    #missing odd size alleles
    #missing non alignment matches with gaps
    return odd_size_alleles, non_alignment_matches, consensus_dict


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
    print (gzipped, type)
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
