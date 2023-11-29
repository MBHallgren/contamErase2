import os
import sys

gene_names = [
    "BACT000001", "BACT000002", "BACT000003", "BACT000004", "BACT000005",
    "BACT000006", "BACT000007", "BACT000008", "BACT000009", "BACT000010",
    "BACT000011", "BACT000012", "BACT000013", "BACT000014", "BACT000015",
    "BACT000016", "BACT000017", "BACT000018", "BACT000019", "BACT000020",
    "BACT000021", "BACT000030", "BACT000031", "BACT000032", "BACT000033",
    "BACT000034", "BACT000035", "BACT000036", "BACT000038", "BACT000039",
    "BACT000040", "BACT000042", "BACT000043", "BACT000044", "BACT000045",
    "BACT000046", "BACT000047", "BACT000048", "BACT000049", "BACT000050",
    "BACT000051", "BACT000052", "BACT000053", "BACT000056", "BACT000057",
    "BACT000058", "BACT000059", "BACT000060", "BACT000061", "BACT000062",
    "BACT000063", "BACT000064", "BACT000065"
]

gene_set = set()

with open('test/NZ.res', 'r') as f:
    for line in f:
        if line.startswith('#'):
            if gene_set != set():
                for item in gene_names:
                    if item not in gene_set:
                        print (item)
            print (line.strip())
            gene_set = set()
        else:
            line = line.strip().split('_')
            gene_set.add(line[0])

for item in gene_names:
    if item not in gene_set:
        print(item)


