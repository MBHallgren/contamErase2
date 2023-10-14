from itertools import combinations

def generate_combinations(input_set):
    all_combinations = set()
    input_list = list(input_set)

    # Generate singles
    all_combinations.update(input_set)

    # Generate pairs, triplets, and so on
    for r in range(2, len(input_list) + 1):
        for combo in combinations(input_list, r):
            all_combinations.add(tuple(sorted(combo)))

    return all_combinations

input_set = {'269_G_C', '334_G_T', '985_G_A', '669_T_C', '494_A_G', '894_T_C', '1103_T_G', '1076_G_C', '665_T_C', '492_A_G', '1105_G_T', '1027_G_A', '282_T_C', '1078_G_A', '263_C_T', '1075_A_G'}

combinations = generate_combinations(input_set)

# Print the result
for combo in combinations:
    print(combo)
