import numpy as np
import pandas as pd
import sys
import os
import argparse
from tqdm import tqdm

# How to use it:
# Standard (2 channels, sorted, downsampled to 100):
# python HMP_csv_to_numpy.py data.csv output.npy --window_h 201 --slide_step 10 --channels 2 --target_samples 100 --sort

# Minimalist (1 channel, no sorting, full samples):
# python HMP_csv_to_numpy.py data.csv output.npy --channels 1


# --- Internal Sorting Import ---
try:
    import helper_haplotypeSorter
except ImportError:
    helper_haplotypeSorter = None

def main():
    parser = argparse.ArgumentParser(description="Convert CSV genomic data to windowed NumPy arrays with downsampling and sorting.")
    parser.add_argument("input_csv", help="Path to input CSV")
    parser.add_argument("output_npy", help="Path for output .npy file")
    parser.add_argument("--window_h", type=int, default=201, help="Window size (sites)")
    parser.add_argument("--slide_step", type=int, default=10, help="Step size for sliding window")
    parser.add_argument("--channels", type=int, choices=[1, 2], default=2, help="1 for Genotype only, 2 for Genotype + Syn/Non-syn")
    parser.add_argument("--target_samples", type=int, default=None, help="Downsample to this many haplotypes (picks ones with least missing data)")
    parser.add_argument("--sort", action="store_true", help="Sort haplotypes by distance (requires haplotype_sorter.py)")
    parser.add_argument("--missing_val", type=int, default=-1, help="Value representing missing data")

    args = parser.parse_args()

    # Generate map filename
    output_prefix = args.output_npy.rsplit('.', 1)[0]
    map_filename = f"{output_prefix}_map.npy"

    print(f"Loading data from {args.input_csv}...")
    try:
        df = pd.read_csv(args.input_csv)
    except Exception as e:
        print(f"Error loading file: {e}")
        sys.exit(1)

    # 1. Extract Data
    site_pos_col = df['site_pos'].values.astype(np.int32)
    genotype_data = df.iloc[:, 2:].values.astype(np.int8) # (Total_Sites, Total_Samples)
    total_sites, total_samples = genotype_data.shape

    # Handle Sample Count
    final_num_samples = args.target_samples if args.target_samples else total_samples
    if args.target_samples and args.target_samples > total_samples:
        print(f"Warning: Target samples ({args.target_samples}) > available ({total_samples}). Using {total_samples}.")
        final_num_samples = total_samples

    # Site Types (Syn vs Non-Syn)
    site_type_map = {'syn': 0, 'nonsyn': 1}
    site_types = df['site_type'].map(site_type_map).fillna(0).values.astype(np.int8)

    # 2. Sliding Window Logic
    num_images = int(np.floor((total_sites - args.window_h) / args.slide_step) + 1)
    
    # Define Shape: (Images, Samples, Sites, [Channels])
    if args.channels == 2:
        final_shape = (num_images, final_num_samples, args.window_h, 2)
    else:
        final_shape = (num_images, final_num_samples, args.window_h)
    
    print(f"Preparing {num_images} windows. Output shape: {final_shape}")
    
    final_data = np.zeros(final_shape, dtype=np.int8)
    final_site_indices = np.zeros((num_images, args.window_h), dtype=np.int32)

    # 3. Processing Loop
    for i in tqdm(range(num_images)):
        start_idx = i * args.slide_step
        end_idx = start_idx + args.window_h
        
        final_site_indices[i] = site_pos_col[start_idx:end_idx]
        
        # Extract Block: (Sites, Samples) -> Transpose to (Samples, Sites)
        window_block = genotype_data[start_idx:end_idx, :].T 
        
        # --- A. DOWNSAMPLING (Quality Based) ---
        if args.target_samples and total_samples > args.target_samples:
            missing_counts = (window_block == args.missing_val).sum(axis=1)
            best_indices = np.argsort(missing_counts)[:args.target_samples]
            best_indices.sort() # Keep relative order
            window_block = window_block[best_indices]
        elif args.target_samples:
            # If target is same or larger, we still use the full block
            pass

        # --- B. CHANNEL CONSTRUCTION ---
        if args.channels == 2:
            current_image = np.zeros((final_num_samples, args.window_h, 2), dtype=np.int8)
            current_image[:, :, 0] = window_block
            current_site_types = site_types[start_idx:end_idx]
            current_image[:, :, 1] = np.tile(current_site_types, (final_num_samples, 1))
        else:
            current_image = window_block.copy()

        # --- C. SORTING ---
        if args.sort:
            if helper_haplotypeSorter is not None:
                # Sorts in-place; for 2 channels, sorter usually looks at channel 0
                helper_haplotypeSorter.sort_haplotypes(current_image, ordering='rows_dist')
            else:
                if i == 0: print("Warning: Sorting requested but haplotype_sorter.py not found.")

        final_data[i] = current_image

    # 4. Save
    np.save(args.output_npy, final_data)
    np.save(map_filename, final_site_indices)

    print(f"\nSuccess! Data saved to {args.output_npy}")

if __name__ == "__main__":
    main()