import numpy as np
import argparse
from tqdm import tqdm
import os

def sort_by_frequency(input_path, output_path):
    # Load with mmap_mode to save RAM
    print(f"Loading {input_path}...")
    data = np.load(input_path, mmap_mode='r')
    num_sims, num_haps, num_snps, num_channels = data.shape
    
    # Pre-allocate output array
    print(f"Allocating output: {output_path}...")
    output = np.lib.format.open_memmap(
        output_path, 
        dtype=np.int8, 
        mode='w+', 
        shape=(num_sims, num_haps, num_snps, num_channels)
    )

    for i in tqdm(range(num_sims), desc="Sorting Simulations"):
        # Load one simulation into memory
        # Format: (num_haps, num_snps, 2)
        sim = np.array(data[i])
        
        # We sort based on Channel 0 (Genotype)
        genotypes = sim[:, :, 0]
        
        # 1. Find unique haplotypes, their first index, and their counts
        # unique_indices points to the first time that specific hap appeared
        unique_haps, first_indices, counts = np.unique(
            genotypes, axis=0, return_index=True, return_counts=True
        )
        
        # 2. Sort the unique groups by frequency (descending)
        sorted_group_order = np.argsort(-counts)
        
        sorted_sim_data = []
        for group_idx in sorted_group_order:
            # Identify the actual haplotype pattern of this group
            pattern = unique_haps[group_idx]
            
            # Find every row in the original sim that matches this pattern
            # (Using genotypes[:,None] logic or simple equality)
            match_indices = np.where((genotypes == pattern).all(axis=1))[0]
            
            # Add all occurrences (full 2-channel data) to our list
            for idx in match_indices:
                sorted_sim_data.append(sim[idx])
        
        # Write the sorted sim to the output memmap
        output[i] = np.array(sorted_sim_data, dtype=np.int8)

    # Flush changes to disk
    del output
    print(f"\nSuccess! Sorted file saved to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sort SLiM numpy arrays by row frequency.")
    parser.add_argument("input", help="Input .npy file")
    parser.add_argument("output", help="Output .npy file")
    args = parser.parse_args()

    if os.path.exists(args.output):
        confirm = input(f"File {args.output} exists. Overwrite? (y/n): ")
        if confirm.lower() != 'y':
            exit()

    sort_by_frequency(args.input, args.output)