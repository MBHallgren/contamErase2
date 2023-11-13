import sys
import os

from nanodecon import kmergenetyper
def number_of_bases_in_file(filename):
    gzipped, type = determine_file_type(filename)
    #print (gzipped, type)
    #determine type#
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


def derive_genes_ids_from_rmlst_candidate(item, arguments, primary_species, candidate_rmlst_dict):
    number = findTemplateNumber(item, arguments.db_dir + "/bac_db")
    identifier = item.split()[0].split('.')[0]
    cmd = 'kma seq2fasta -seqs {} -t_db {} > {}/{}.fasta' \
        .format(number, arguments.db_dir + '/bac_db', arguments.output + '/rmlst_' + primary_species + '/fastas',
                identifier)
    os.system(cmd)

    kmergenetyper.KmergenetyperRunner(None,
                                      None,
                                      '{}/{}.fasta'.format(
                                          arguments.output + '/rmlst_' + primary_species + '/fastas', identifier),
                                      '-t_db {} -o {} -md {} -t {}' \
                                      .format(arguments.db_dir + '/rmlst_db',
                                              arguments.output + '/rmlst_' + primary_species + '/fastas/' + identifier,
                                              0.1,
                                              arguments.threads)
                                      ).run()

    mlst_genes = list()
    gene_lengths = list()
    with open(
            arguments.output + '/rmlst_' + primary_species + '/fastas/' + identifier + '/' + identifier + '.res',
            'r') as f:
        for line in f:
            if not line.startswith('#'):
                mlst_genes.append(line.split('\t')[0])
                gene_lengths.append(int(line.split('\t')[3]))
    candidate_rmlst_dict[item] = {}
    for i in range(len(mlst_genes)):
        candidate_rmlst_dict[item][mlst_genes[i]] = [False, gene_lengths[i]]  #initialize as unconfirmed
    return candidate_rmlst_dict


def findTemplateNumber(name, database):
    if os.path.exists(database + ".name"):
        with open(database + ".name") as f:
            t = 1
            for line in f:
                if line.rstrip().startswith(name):
                    return t
                else:
                    t += 1