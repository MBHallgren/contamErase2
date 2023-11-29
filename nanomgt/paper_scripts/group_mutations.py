import pandas as pd

def determine_group(df):
    # Convert 'Position' to integers
    df['Position'] = pd.to_numeric(df['Position'], errors='coerce')

    # Convert 'MutationComment' to string
    df['MutationComment'] = df['MutationComment'].astype(str)

    # Convert 'Gene' to string and modify each gene name to only include the gene number
    df['Gene'] = df['Gene'].astype(str).apply(lambda gene: gene.split('_')[0][-2:])

    # Sort by Gene and Position to make the proximity check easier
    df = df.sort_values(by=['Gene', 'Position'])

    # Add a column for proximity (within 5 bases of another mutation)
    df['Proximity'] = False
    for i in range(1, len(df)):
        if df.iloc[i]['Gene'] == df.iloc[i-1]['Gene'] and abs(df.iloc[i]['Position'] - df.iloc[i-1]['Position']) <= 5:
            df.at[i, 'Proximity'] = True
            df.at[i-1, 'Proximity'] = True

    # Add a column for Novelty based on MutationComment
    df['Novel'] = df['MutationComment'].str.contains("Novel mutation")

    # Determine the group for each mutation
    df['Group'] = df.apply(lambda row: 1 if row['Novel'] and row['Proximity'] else
                                       2 if not row['Novel'] and row['Proximity'] else
                                       3 if row['Novel'] and not row['Proximity'] else
                                       4, axis=1)

    return df

# Load the data from a file
file_path = 'error1.csv'  # Replace with your file path
df = pd.read_csv(file_path)

# Apply the function to determine the group
grouped_df = determine_group(df)

# Count the number of mutations in each group
group_counts = grouped_df['Group'].value_counts().sort_index()

# Save the result
grouped_df.to_csv('output_modified.csv', index=False)  # This will save the output to a file named 'output_modified.csv'

# Print the count of each group
print(group_counts)
