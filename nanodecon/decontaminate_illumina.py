import os
import gzip
import sys
from nanodecon import kma
from nanodecon import util

from nanodecon.intra_species_detection import determine_intra_species_contamination_illumina

def illumina_decontamination(arguments):
    input_string = " ".join(arguments.illumina)

    os.system('mkdir ' + arguments.output + '/output')
    #os.system('mkdir ' + arguments.output + '/output/contaminations')
    #os.system('mkdir ' + arguments.output + '/alignments')

    #kma.KMARunner(input_string,
    #              arguments.output + "/bacteria_alignment",
    #              arguments.db_dir + "/bac_db",
    #              "-mem_mode -1t1 -t {} -ID 10".format(arguments.threads)).run()

    total_bacteria_aligning_bases = util.number_of_bases_in_file(arguments.output + "/bacteria_alignment.fsa")
    primary, candidate_dict = drive_bacteria_results(arguments, total_bacteria_aligning_bases)

    primary_species = primary.split()[1] + ' ' + primary.split()[2]
    print (primary_species)
    print (primary_species)
    #produce_species_specific_kma_db(primary_species,
    #                                '/home/people/malhal/contamErase_db/rmlst.fsa',
    #                                '/home/people/malhal/contamErase_db/rmlst_scheme.txt',
    #                                arguments.output)
    #kma.KMARunner(input_string,
    #             arguments.output + "/rmlst_alignment",
    #             arguments.output + '/specie_db',
    #             "-1t1 -t {} -ID 10 -md 1.5 -matrix -vcf -oa".format(arguments.threads)).run()
    #os.system('gunzip ' + arguments.output + '/rmlst_alignment.mat.gz')

    odd_size_alleles, non_alignment_matches, consensus_dict = build_consensus_dict(arguments, arguments.output + '/rmlst_alignment.res', arguments.output + '/rmlst_alignment.mat')
    sys.exit()
    allele_lengths = check_allele_lengths(arguments.output)
    for item in allele_lengths:
        print (item, allele_lengths[item])
    sys.exit()

    #black_list_plasmid, black_list_viral, black_list_human = derive_non_bacterial_black_list(arguments.output)
    rmlst_candidates = derive_rmlst_candidates(primary, candidate_dict)


    derive_read_pools_illumina(candidate_dict, arguments, primary, black_list_human)


    if len(rmlst_candidates) >= 2:
        candidate_rmlst_dict_results = determine_intra_species_contamination_illumina(arguments, rmlst_candidates, primary, candidate_dict, total_bacteria_aligning_bases, arguments.ugt)
    else:
        candidate_rmlst_dict_results = None
    sys.exit()
    produce_final_output_illumina(arguments, arguments.output + '/bacteria_alignment.frag', primary, candidate_rmlst_dict_results, black_list_plasmid, black_list_viral, black_list_human)
    #produce_contamination_report #TBD
    sys.exit()

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
                print(line)
                allele = line[1:]
                gene = allele.split('_')[0]
                if allele in correct_size_alleles:
                    print ('found correct size allele')
                    correct_size_flag = True
                    index = 0
                else:
                    correct_size_flag = False
            else:
                if line != '':
                    if correct_size_flag:
                        line = line.split('\t')
                        line = line[1:]
                        print (index)
                        print (allele)
                        for i in range(len(line)):
                            consensus_dict[gene][index][i] += int(line[i])
                            if i == 5 and int(line[i]) > 0:
                                print (line)
                        index += 1
    for item in consensus_dict:
        print (item, consensus_dict[item])






    return odd_size_alleles, non_alignment_matches, consensus_dict

def check_allele_lengths(output):
    allele_lengths = {}
    with open(output + '/rmlst_alignment.res', 'r') as f:
        for line in f:
            if not line.startswith('#'):
                allele = line.strip().split('\t')[0].split('_')[0]
                if allele not in allele_lengths:
                    allele_lengths[allele] = set()
                    allele_lengths[allele].add(line.strip().split('\t')[3])
                else:
                    allele_lengths[allele].add(line.strip().split('\t')[3])
    return allele_lengths


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

def produce_final_output_illumina(arguments, frag_file, primary, candidate_rmlst_dict_results, black_list_plasmid, black_list_viral, black_list_human):
    primary_species = primary.split()[1].lower() + '_' + primary.split()[2].lower()
    rmlst_hits = output_primary_reads(arguments, frag_file, primary, candidate_rmlst_dict_results, black_list_plasmid, primary_species)
    if len(rmlst_hits) >= 1:
        output_rmlst_hits(arguments, frag_file, rmlst_hits)
    if len(black_list_plasmid) > 0:
        output_black_list_plasmid(arguments, black_list_plasmid)
    if len(black_list_viral) > 0:
        output_black_list_viral(arguments, black_list_viral)
    if len(black_list_human) > 0:
        output_black_list_human(arguments, black_list_human)
    clean_up(arguments)

def output_black_list_human(arguments, black_list_human):
    human_reads_1 = set()
    human_reads_2 = set()
    for item in black_list_human:
        if item.endswith('1'):
            human_reads_1.add(item)
        else:
            human_reads_2.add(item)

    with open('{}/{}_R1.txt'.format(arguments.output, 'human_reads'), 'w') as write_file:
        for item in human_reads_1:
            write_file.write(item + '\n')
    with open('{}/{}_R2.txt'.format(arguments.output, 'human_reads'), 'w') as write_file:
        for item in human_reads_2:
            write_file.write(item + '\n')
    os.system('seqtk subseq {} {} > {}/output/contaminations/human_reads_1.fastq'.format(arguments.illumina[0], arguments.output + '/human_reads_R1.txt', arguments.output))
    os.system('seqtk subseq {} {} > {}/output/contaminations/human_reads_2.fastq'.format(arguments.illumina[1], arguments.output + '/human_reads_R2.txt', arguments.output))

def output_black_list_viral(arguments, black_list_viral):
    viral_reads_1 = set()
    viral_reads_2 = set()
    for item in black_list_viral:
        if item.endswith('1'):
            viral_reads_1.add(item)
        else:
            viral_reads_2.add(item)

    with open('{}/{}_R1.txt'.format(arguments.output, 'viral_reads'), 'w') as write_file:
        for item in viral_reads_1:
            write_file.write(item + '\n')
    with open('{}/{}_R2.txt'.format(arguments.output, 'viral_reads'), 'w') as write_file:
        for item in viral_reads_2:
            write_file.write(item + '\n')
    os.system('seqtk subseq {} {} > {}/output/contaminations/viral_reads_1.fastq'.format(arguments.illumina[0], arguments.output + '/viral_reads_R1.txt', arguments.output))
    os.system('seqtk subseq {} {} > {}/output/contaminations/viral_reads_2.fastq'.format(arguments.illumina[1], arguments.output + '/viral_reads_R2.txt', arguments.output))

def output_black_list_plasmid(arguments, black_list_plasmid):
    plasmid_reads_1 = set()
    plasmid_reads_2 = set()
    for item in black_list_plasmid:
        if item.endswith('1'):
            plasmid_reads_1.add(item)
        else:
            plasmid_reads_2.add(item)

    with open('{}/{}_R1.txt'.format(arguments.output, 'plasmid_reads'), 'w') as write_file:
        for item in plasmid_reads_1:
            write_file.write(item + '\n')
    with open('{}/{}_R2.txt'.format(arguments.output, 'plasmid_reads'), 'w') as write_file:
        for item in plasmid_reads_2:
            write_file.write(item + '\n')
    os.system('seqtk subseq {} {} > {}/output/contaminations/plasmid_reads_1.fastq'.format(arguments.illumina[0], arguments.output + '/plasmid_reads_R1.txt', arguments.output))
    os.system('seqtk subseq {} {} > {}/output/contaminations/plasmid_reads_2.fastq'.format(arguments.illumina[1], arguments.output + '/plasmid_reads_R2.txt', arguments.output))

def clean_up(arguments):
    os.system('rm -r ' + arguments.output + '/*.fastq')
    os.system('mv ' + arguments.output + '/*.frag ' + arguments.output + '/alignments/.')
    os.system('mv ' + arguments.output + '/*.fsa ' + arguments.output + '/alignments/.')
    os.system('mv ' + arguments.output + '/*.res ' + arguments.output + '/alignments/.')
    os.system('mv ' + arguments.output + '/*.aln ' + arguments.output + '/alignments/.')
    os.system('tar -czvf ' + arguments.output + '/alignments.tar.gz ' + arguments.output + '/alignments')
    os.system('rm -r ' + arguments.output + '/alignments')

def output_rmlst_hits(arguments, frag_file, rmlst_hits):
    intra_species_reads_1 = set()
    intra_species_reads_2 = set()
    with open(frag_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[-2] in rmlst_hits:
                if line[-1].endswith('1'):
                    intra_species_reads_1.add(line[-1])
                elif line[-1].endswith('2'):
                    intra_species_reads_2.add(line[-1])
    with open('{}/{}_R1.txt'.format(arguments.output, 'intra_species_contamination'), 'w') as write_file:
        for item in intra_species_reads_1:
            write_file.write(item + '\n')
    with open('{}/{}_R2.txt'.format(arguments.output, 'intra_species_contamination'), 'w') as write_file:
        for item in intra_species_reads_2:
            write_file.write(item + '\n')
    os.system('seqtk subseq {} {} > {}/output/contaminations/intra_species_contamination_R1.fastq'.format(arguments.illumina[0], arguments.output + '/intra_species_contamination_R1.txt', arguments.output))
    os.system('seqtk subseq {} {} > {}/output/contaminations/intra_species_contamination_R2.fastq'.format(arguments.illumina[1], arguments.output + '/intra_species_contamination_R2.txt', arguments.output))

def output_primary_reads(arguments, frag_file, primary, candidate_rmlst_dict_results, black_list_plasmid, primary_species):
    black_list_templates = set()
    if candidate_rmlst_dict_results is not None:
        for item in candidate_rmlst_dict_results:
            if item != primary and candidate_rmlst_dict_results[item]['rmlst_id'] == True:
                black_list_templates.add(item)

        blast_list_reads = set()
        with open(frag_file, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if line[-2] not in black_list_templates:
                    blast_list_reads.add(line[-1])

        read_set_1 = set()
        read_set_2 = set()
        with open('{}/{}_read_pool_R1.fastq'.format(arguments.output, primary_species), 'r') as f:
            for line in f:
                if line[0] == '@':
                    read = line.strip()[1:]
                    if read not in blast_list_reads:
                        read_set_1.add(read)
        with open('{}/{}_read_pool_R2.fastq'.format(arguments.output, primary_species), 'r') as f:
            for line in f:
                if line[0] == '@':
                    read = line.strip()[1:]
                    if read not in blast_list_reads:
                        read_set_2.add(read)

        with open('{}/{}.txt'.format(arguments.output, primary_species), 'w') as write_file:
            for item in read_set_1:
                write_file.write(item + '\n')
        with open('{}/{}.txt'.format(arguments.output, primary_species), 'w') as write_file:
            for item in read_set_2:
                write_file.write(item + '\n')
        os.system('seqtk subseq {} {} > {}/output/contaminations/{}_R1.fastq'.format(arguments.illumina[0], arguments.output + '/' + primary_species + '.txt', arguments.output, primary_species))
        os.system('seqtk subseq {} {} > {}/output/contaminations/{}_R2.fastq'.format(arguments.illumina[1], arguments.output + '/' + primary_species + '.txt', arguments.output, primary_species))
    else:
        os.system('mv {}/{}_read_pool_R1.fastq {}/output/{}_R1.fastq'.format(arguments.output, primary_species, arguments.output, primary_species))
        os.system('mv {}/{}_read_pool_R2.fastq {}/output/{}_R2.fastq'.format(arguments.output, primary_species, arguments.output, primary_species))

    return black_list_templates


def derive_all_typed_rmlst_genes_from_input_reads(filename):
    all_rmlst_genes = dict()
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                all_rmlst_genes[line.split('\t')[0]] = line.split('\t')[8]
    return all_rmlst_genes
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


def derive_read_pools_illumina(candidate_dict, arguments, primary, black_list_human):
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
            with open('{}/{}_R1.txt'.format(arguments.output, item.replace(' ', '_').lower()), 'w') as write_file:
                for read in pool_dict[item]:
                    if read.endswith('/1'):
                        write_file.write(read + '\n')

            with open('{}/{}_R2.txt'.format(arguments.output, item.replace(' ', '_').lower()), 'w') as write_file:
                for read in pool_dict[item]:
                    if read.endswith('/2'):
                        write_file.write(read + '\n')
            os.system('seqtk subseq {} {} > {}/{}_R1.fastq'.format(arguments.illumina[0], '{}/{}_R1.txt'.format(arguments.output, item.replace(' ', '_').lower()), arguments.output, item.replace(' ', '_').lower()))
            os.system('seqtk subseq {} {} > {}/{}_R2.fastq'.format(arguments.illumina[1], '{}/{}_R2.txt'.format(arguments.output, item.replace(' ', '_').lower()), arguments.output, item.replace(' ', '_').lower()))
            for read in pool_dict[item]:
                black_list_human.add(read)

    primary_species_reads_1 = set()
    primary_species_reads_2 = set()
    with open('{}'.format(arguments.illumina[0]), 'r') as f:
        for line in f:
            if line[0] == '@':
                read = line.strip()[1:]
                if read not in black_list_human:
                    primary_species_reads_1.add(read)
    with open('{}'.format(arguments.illumina[1]), 'r') as f:
        for line in f:
            if line[0] == '@':
                read = line.strip()[1:]
                if read not in black_list_human:
                    primary_species_reads_2.add(read)
    with open('{}/{}_R1.txt'.format(arguments.output, primary_species.replace(' ', '_').lower()), 'w') as write_file:
        for read in primary_species_reads_1:
            write_file.write(read + '\n')
    with open('{}/{}_R2.txt'.format(arguments.output, primary_species.replace(' ', '_').lower()), 'w') as write_file:
        for read in primary_species_reads_2:
            write_file.write(read + '\n')
    os.system('seqtk subseq {} {} > {}/{}_read_pool_R1.fastq'.format(arguments.illumina[0], '{}/{}_R1.txt'.format(arguments.output, primary_species.replace(' ', '_').lower()), arguments.output, primary_species.replace(' ', '_').lower()))
    os.system('seqtk subseq {} {} > {}/{}_read_pool_R2.fastq'.format(arguments.illumina[1], '{}/{}_R2.txt'.format(arguments.output, primary_species.replace(' ', '_').lower()), arguments.output, primary_species.replace(' ', '_').lower()))


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
