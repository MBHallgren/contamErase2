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
#i_list = ['78_T_C', '193_G_T', '195_C_A', '196_G_T', '197_A_G', '201_T_C', '264_T_G', '356_C_T', '410_T_C', '708_T_A', '722_A_G', '828_T_C', '1071_T_C', '1072_G_T', '1077_T_C', '1078_G_T', '1082_C_T', '1085_C_G', '1095_G_A', '1103_T_C', '1109_G_A']
#i_list = ['5_C_T', '15_C_A', '19_G_A', '20_A_G', '30_T_C', '64_G_A', '66_G_A', '117_A_C', '121_A_G', '125_G_C', '126_G_A', '356_C_T', '392_G_A', '397_G_A', '400_C_A', '402_G_T', '626_C_T', '628_T_G', '670_G_A', '676_T_A', '708_T_A', '828_T_C', '855_C_G', '1017_G_A', '1109_G_A']
i_list = ['30_T_C', '64_G_A', '66_G_A', '117_A_C', '121_A_G', '125_G_C', '126_G_A', '356_C_T', '392_G_A', '397_G_A', '400_C_A', '402_G_T', '626_C_T', '628_T_G', '670_G_A', '676_T_A', '708_T_A', '828_T_C', '855_C_G', '1017_G_A', '1109_G_A']

input_set = set(i_list)
combinations = generate_combinations(input_set)

# Print the result
for combo in combinations:
    print(combo)
