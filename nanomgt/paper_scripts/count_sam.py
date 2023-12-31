import os
import sys


def count_alignments_and_zero_pos(sam_file):
    total_alignments = 0
    zero_pos_count = 0
    total_cov = 0

    with open(sam_file, 'r') as f:
        for line in f:
            # Skip header lines
            if line.startswith('@'):
                continue
            total_alignments += 1

            # Splitting each line by tab to get the fields
            fields = line.strip().split('\t')

            # The 4th field in a SAM file is the leftmost mapping POSition
            pos = int(fields[3])
            tlen = int(fields[8])
            query_length = len(fields[9])
            print (pos)
            if pos == 1:
                zero_pos_count += 1
                if query_length > tlen:
                    total_cov += 1

    print(f'Total Alignments: {total_alignments}')
    print(f'Total Reads with POS=0: {zero_pos_count}')
    print(f'Total Coverage: {total_cov}')

# Example usage:
count_alignments_and_zero_pos('test/test.sam')
