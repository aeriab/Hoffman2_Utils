import numpy as np
import argparse
from tqdm import tqdm
import os

def sort_by_frequency(input_path, output_path, window_size):
    # Load with mmap_mode to save RAM
    print(f"Loading {input_path}...")
    data = np.load(input_path, mmap_mode='r')
    num_sims, num_haps, num_snps, num_channels = data.shape
    
    # Calculate middle window indices
    if window_size > num_snps:
        print(f"Warning: Window size {window_size} exceeds SNP count {num_snps}. Sorting by all SNPs.")
        start, end = 0, num_snps
    else:
        start = (num_snps - window_size) // 2
        end = start + window_size
    
    print(f"Sorting based on middle {end-start} SNPs (indices {start} to {end}).")

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
        sim = np.array(data[i])
        
        # We define the "Identity" of the row by the middle slice of Channel 0
        genotypes = sim[:, :, 0]
        sorting_slice = genotypes[:, start:end]
        
        # 1. Find unique patterns in the WINDOW, their first index, and counts
        unique_patterns, first_indices, counts = np.unique(
            sorting_slice, axis=0, return_index=True, return_counts=True
        )
        
        # 2. Sort the unique groups by frequency (descending)
        sorted_group_order = np.argsort(-counts)
        
        sorted_sim_data = []
        for group_idx in sorted_group_order:
            pattern = unique_patterns[group_idx]
            
            # Find every row where the MIDDLE WINDOW matches this pattern
            match_indices = np.where((sorting_slice == pattern).all(axis=1))[0]
            
            # Add all occurrences (full 2-channel data) to our list
            for idx in match_indices:
                sorted_sim_data.append(sim[idx])
        
        # Write the sorted sim to the output memmap
        output[i] = np.array(sorted_sim_data, dtype=np.int8)

    # Flush changes to disk
    del output
    print(f"\nSuccess! Sorted file saved to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sort SLiM numpy arrays by middle-window row frequency.")
    parser.add_argument("input", help="Input .npy file")
    parser.add_argument("output", help="Output .npy file")
    parser.add_argument("--window", type=int, default=100, help="Number of middle SNPs to use for sorting identity (default: 100)")
    args = parser.parse_args()

    if os.path.exists(args.output):
        confirm = input(f"File {args.output} exists. Overwrite? (y/n): ")
        if confirm.lower() != 'y':
            exit()

    sort_by_frequency(args.input, args.output, args.window)