import os
import sys

def quick_determine_multiple_rmlst_alleles_from_res_files(res_file, cutoff):
    rmlst_dict = {}
    with open(res_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                gene = line.split()[0].split('_')[0]
                if gene not in rmlst_dict:
                    depth = float(line.split('\t')[-3])
                    if depth >= cutoff:
                        rmlst_dict[gene] = {}
                        allele = line.split()[0]
                        rmlst_dict[gene][allele] = depth
                else:
                    depth = float(line.split('\t')[-3])
                    if depth >= cutoff:
                        allele = line.split()[0]
                        rmlst_dict[gene][allele] = depth
    return rmlst_dict

mlst_hits = quick_determine_multiple_rmlst_alleles_from_res_files('rmlst_escherichia_coli.res', 1)
for item in mlst_hits:
    print(item, mlst_hits[item])