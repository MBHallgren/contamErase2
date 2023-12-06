import os
import sys

csv_file = 'Metadata_multi_country.csv'
kmers = 'source_9mer.matrix'


nine_mers = None
extended_headers = list()

header_dict = dict()

true_headers = list()

with open(kmers, 'r') as f:
    for line in f:
        if line.startswith('genome'):
            line = line.strip().split(',')
            nine_mers = line[1:]
        else:
            line = line.strip().split(',')
            extended_headers.append(line[0])

print (len(extended_headers))

with open (csv_file, 'r') as f:
    for line in f:
        if not line.startswith('New_Sou'):
            line = line.strip().split(';')
            true_headers.append(line[1])

print (len(true_headers))


for true_header in true_headers:
    for extended_header in extended_headers:
        if true_header in extended_header:
            header_dict[true_header] = extended_header

missing = list()

for true_header in true_headers:
    if true_header not in header_dict:
        missing.append(true_header)

#for extended_header in extended_headers:
#    if extended_header not in header_dict.values():
#        missing.append(extended_header)

print (len(header_dict))

for item in missing:
    print (item)

#print (header_dict)
#for key, value in header_dict.items():
#    print (key + '\t' + value)