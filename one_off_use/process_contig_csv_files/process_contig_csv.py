import pandas as pd
import os

# 1. Load the list of good samples from the text file
with open('Alistipes_finegoldii_56071_good_samples.txt', 'r') as f:
    # strip() removes newlines and whitespace
    good_samples = [line.strip() for line in f if line.strip()]

# Metadata columns to keep at the start of every file
info_cols = ['contig', 'gene_id', 'site_pos', 'site_type']

file_data = []

# 2. Process each of the 86 files
for i in range(1, 87):
    # Construct filename (e.g., AENZ01000001_haplotypes.csv)
    input_fname = f"AENZ01{i:06d}_haplotypes.csv"
    
    if not os.path.exists(input_fname):
        print(f"Warning: {input_fname} not found. Skipping...")
        continue

    # Read the CSV
    df = pd.read_csv(input_fname, low_memory=False)
    
    # Strictly include ONLY the metadata and the samples from your txt file
    # We check if the sample exists in the current dataframe columns
    target_cols = [col for col in info_cols if col in df.columns] + \
                  [col for col in good_samples if col in df.columns]
    
    filtered_df = df[target_cols].copy()
    
    # Store the number of sites (rows) and the filtered dataframe
    num_sites = len(filtered_df)
    file_data.append({
        'num_sites': num_sites,
        'df': filtered_df,
        'original_name': input_fname
    })

# 3. Sort by number of sites descending (largest/most sites first)
file_data.sort(key=lambda x: x['num_sites'], reverse=True)

# 4. Save the files with the new names
print("\nRenaming and saving files...")
for index, item in enumerate(file_data, start=1):
    # Pad index with a zero for 01, 02... 86
    new_name = f"A_finegoldii_{index:02d}.csv"
    
    item['df'].to_csv(new_name, index=False)
    print(f"{item['original_name']} -> {new_name} ({item['num_sites']} sites)")

print("\nProcessing complete.")