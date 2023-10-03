import os
import sys

secondary_true_positives = 0
other_random_hit = 0
secondary_hit_should_have_been_primary = 0

cont_1 = 'NZ_CP087110.1'
cont_2 = 'NZ_CP079993.1'
primary = 'NC_008253.1'
with open('intra_test/bacteria_alignment.frag' , 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        if line[-2].startswith(cont_1) or line[-2].startswith(cont_2) and not line[-1].startswith(primary):
            if line[-1].startswith('NZ_CP026027'):
                secondary_true_positives += 1
            else:
                other_random_hit += 1
        elif line[-2].startswith(cont_1) or line[-2].startswith(cont_2) and line[-1].startswith(primary):
            secondary_hit_should_have_been_primary +=1

print ('secondary true hits:', secondary_true_positives )
print ('secondary hits from a primary read:', secondary_hit_should_have_been_primary)
print ('klebs random hit:', other_random_hit)
print ('')
