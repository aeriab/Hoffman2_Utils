import numpy as np
import argparse
from tqdm import tqdm
import sys
import os

def determine_maj(image):
    # (samples, sites, channels)
    geno = image[:, :, 0]
    color = image[:, :, 1]
    
    num_samples = geno.shape[0]
    num_sites = geno.shape[1]

    for col in range(num_sites):
        col_geno = geno[:, col]
        col_color = color[:, col]
        
        count_ones = np.sum(col_geno == 1)
        
        if count_ones > (num_samples / 2):
            # Identify Site Type (-1 = syn, 1 = nonsyn)
            site_type = 0
            if np.any(col_color == -1):
                site_type = -1
            elif np.any(col_color == 1):
                site_type = 1

            # Flip Genotype
            new_col_geno = np.where(col_geno == 1, -1, 1)
            image[:, col, 0] = new_col_geno

            # Update Color: Only label the new minor allele (1)
            new_col_color = np.zeros_like(col_color)
            new_col_color[new_col_geno == 1] = site_type
            image[:, col, 1] = new_col_color
            
    return image

def entire_arr_maj(all_images):
    for i in tqdm(range(all_images.shape[0]), desc="Flipping Alleles"):
        all_images[i] = determine_maj(all_images[i])
    return all_images

def sort_by_frequency(data_array, output_path, window_size):
    """
    Takes the processed NUMPY ARRAY directly instead of loading from disk.
    """
    num_sims, num_haps, num_snps, num_channels = data_array.shape
    
    # Safety check for window size vs total SNPs
    if window_size > num_snps:
        print(f"Warning: Window size {window_size} exceeds SNP count ({num_snps}). Sorting by all.")
        start, end = 0, num_snps
    else:
        start = (num_snps - window_size) // 2
        end = start + window_size
    
    print(f"Sorting based on middle {end-start} SNPs (indices {start} to {end}).")

    # Use memmap for the output to handle large files safely
    output = np.lib.format.open_memmap(
        output_path, 
        dtype=np.int8, 
        mode='w+', 
        shape=(num_sims, num_haps, num_snps, num_channels)
    )

    for i in tqdm(range(num_sims), desc="Sorting Simulations"):
        sim = data_array[i] # This is already the flipped data
        
        genotypes = sim[:, :, 0]
        sorting_slice = genotypes[:, start:end]
        
        unique_patterns, first_indices, counts = np.unique(
            sorting_slice, axis=0, return_index=True, return_counts=True
        )
        
        sorted_group_order = np.argsort(-counts)
        
        sorted_sim_data = []
        for group_idx in sorted_group_order:
            pattern = unique_patterns[group_idx]
            match_indices = np.where((sorting_slice == pattern).all(axis=1))[0]
            
            for idx in match_indices:
                sorted_sim_data.append(sim[idx])
        
        output[i] = np.array(sorted_sim_data, dtype=np.int8)

    output.flush() # Ensure data is written
    del output
    print(f"Success! Sorted file saved to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Flip majority alleles and sync color channels, then sort.")
    parser.add_argument("input_npy", help="Path to input .npy")
    parser.add_argument("output_npy", help="Path to output .npy")
    # Set default to None so we can calculate the dynamic default later
    parser.add_argument("--window", type=int, default=None, 
                        help="Middle SNP window size (default: 50 percent of SNPs, min 100)")
    args = parser.parse_args()

    if os.path.exists(args.output_npy):
        print(f"Warning: File {args.output_npy} already exists. Overwriting...")

    # if os.path.exists(args.output_npy):
    #     confirm = input(f"File {args.output_npy} exists. Overwrite? (y/n): ")
    #     if confirm.lower() != 'y':
    #         sys.exit("Operation cancelled.")

    print(f"Loading {args.input_npy} into memory...")
    data = np.load(args.input_npy)
    
    # Calculate the dynamic window size if not explicitly provided
    num_snps = data.shape[2]
    if args.window is None:
        # 50% of total SNPs, but at least 100 (capped at num_snps total)
        calculated_window = max(100, int(num_snps * 0.5))
        window_size = min(num_snps, calculated_window)
        print(f"No window specified. Using default: {window_size} (middle 50% or 100 min).")
    else:
        window_size = args.window

    print("Processing images (flipping major/minor)...")
    processed_data = entire_arr_maj(data)

    print(f"Sorting and saving to {args.output_npy}...")
    sort_by_frequency(processed_data, args.output_npy, window_size)

if __name__ == "__main__":
    main()