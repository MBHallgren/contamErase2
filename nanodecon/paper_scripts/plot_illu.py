import os
import sys

path = '/home/people/malhal/contamErase/data/illumina/intra/paper/results/'
files = os.listdir(path)
original_input_path = '/home/people/malhal/contamErase/data/illumina/intra/paper/subsets/'

for item in files:
    if item != 't.py':
        print (item)
        if os.path.exists(path + item + '/alignments.tar.gz'):
            os.system('tar -xzf {}/{}/alignments.tar.gz -C {}/{}'.format(path, item, path, item))
            os.system('rm {}/{}/alignments.tar.gz'.format(path, item))

        alignment_path = path + item + '/' + item + '/alignments/'
        output_path = path + item + '/output/'
        original_1_read = original_input_path + item + '/' + item + '_R1.fastq'
        original_2_read = original_input_path + item + '/' + item + '_R2.fastq'
        primary_1_read = original_input_path + item + '/primary_R1.fastq'
        primary_2_read = original_input_path + item + '/primary_R2.fastq'
        secondary_1_read = original_input_path + item + '/secondary_R1.fastq'
        secondary_2_read = original_input_path + item + '/secondary_R2.fastq'
        extra_1_read = original_input_path + item + '/extra_R1.fastq'
        extra_2_read = original_input_path + item + '/extra_R2.fastq'

        primary_1_read_set = set()
        primary_2_read_set = set()

        secondary_1_read_set = set()
        secondary_2_read_set = set()

        extra_1_read_set = set()
        extra_2_read_set = set()

        unmapped_primary_reads = set()
        unmapped_secondary_reads = set()
        unmapped_extra_reads = set()
        print (primary_1_read)

        t1 = 0
        with open(primary_1_read, 'r') as f:
            for line in f:
                if line.startswith('@N'):
                    primary_1_read_set.add(line[1:].strip())
                t1 += 1

        print (len(primary_1_read_set))

        t2 = 0
        with open(primary_2_read, 'r') as f:
            for line in f:
                if line.startswith('@N'):
                    primary_2_read_set.add(line[1:].strip())
                t2 += 1

        with open(secondary_1_read, 'r') as f:
            for line in f:
                if line.startswith('@N'):
                    secondary_1_read_set.add(line[1:].strip())

        with open(secondary_2_read, 'r') as f:
            for line in f:
                if line.startswith('@N'):
                    secondary_2_read_set.add(line[1:].strip())

        with open(extra_1_read, 'r') as f:
            for line in f:
                if line.startswith('@N'):
                    extra_1_read_set.add(line[1:].strip())

        with open(extra_2_read, 'r') as f:
            for line in f:
                if line.startswith('@N'):
                    extra_2_read_set.add(line[1:].strip())

        all_reads = primary_1_read_set.union(primary_2_read_set, secondary_1_read_set, secondary_2_read_set,
                                             extra_1_read_set, extra_2_read_set)
        print ('all reads')
        print (len(all_reads))
        print ('primary reads')
        print (len(primary_1_read_set.union(primary_2_read_set)))
        print ('secondary reads')
        print (len(secondary_1_read_set.union(secondary_2_read_set)))
        print ('extra reads')
        print (len(extra_1_read_set.union(extra_2_read_set)))


        mapped_reads = set()

        with open(alignment_path + 'bacteria_alignment.frag', 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                mapped_reads.add(line[-1])

        for item in all_reads:
            if item not in mapped_reads:
                if item.startswith('NC_008253'):
                    unmapped_primary_reads.add(item)
                if item.startswith('NZ_CP026027'):
                    unmapped_secondary_reads.add(item)
                if item.startswith('NB502092'):
                    unmapped_extra_reads.add(item)

        results_1_primary_read = output_path + '1_R1.fastq'
        results_2_primary_read = output_path + '1_R2.fastq'

        secodary_1_read = output_path + 'intra_species_contamination_R1.fastq'
        secodary_2_read = output_path + 'intra_species_contamination_R2.fastq'

        extra_1_read = output_path + 'contaminations/staphylococcus_aureus_read_pool_R1.fastq'
        extra_2_read = output_path + 'contaminations/staphylococcus_aureus_read_pool_R2.fastq'

        results_1_primary_read_set = set()
        results_2_primary_read_set = set()

        secodary_1_read_set = set()
        secodary_2_read_set = set()

        extra_1_read_set = set()
        extra_2_read_set = set()

        print ('Unmapped')
        print (len(unmapped_primary_reads))
        print (len(unmapped_secondary_reads))
        print (len(unmapped_extra_reads))

        #Here derive stats to tsv file



        header_list = ['sample', 'primary_reads', 'secondary_reads', 'extra_reads', 'unmapped_primary_reads', 'unmapped_secondary_reads', 'unmapped_extra_reads']
        stats_dict = {}
        stats_dict['sample'] = item
        stats_dict['primary_reads'] = primary_1_read_set.union(primary_2_read_set)
        stats_dict['secondary_reads'] = secondary_1_read_set.union(secondary_2_read_set)
        stats_dict['extra_reads'] = extra_1_read_set.union(extra_2_read_set)
        stats_dict['unmapped_primary_reads'] = unmapped_primary_reads
        stats_dict['unmapped_secondary_reads'] = unmapped_secondary_reads
        stats_dict['unmapped_extra_reads'] = unmapped_extra_reads

        stats_dict['primary_reads_hits'] = 0
        for read in mapped_reads:
            if read.startswith('NC_008253'):
                stats_dict['primary_reads_hits'] += 1

        print (len(stats_dict['primary_reads']))
        print (stats_dict['primary_reads_hits'])

        sys.exit()
