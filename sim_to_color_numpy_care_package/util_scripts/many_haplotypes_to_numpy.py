import os
import sys
import numpy as np
import pandas as pd
import argparse
from tqdm import tqdm
import glob
import re

# Import the helper if available
try:
    import helper_haplotypeSorter
except ImportError:
    helper_haplotypeSorter = None
    print("Warning: helper_haplotypeSorter.py not found. Sorting will be skipped.")

def recode_genotype(raw_block):
    """Recode: 0=major -> -1, 1=minor -> 1, -1=missing -> 0"""
    recoded = np.zeros_like(raw_block, dtype=np.int8)
    recoded[raw_block == 0] = -1
    recoded[raw_block == 1] = 1
    recoded[raw_block == -1] = 0
    return recoded

def build_color_channel(recoded_geno, site_types_window):
    """Channel 1: -1=syn, 0=major/missing, 1=nonsyn"""
    n_samples, n_sites = recoded_geno.shape
    color = np.zeros((n_samples, n_sites), dtype=np.int8)
    is_minor = (recoded_geno == 1)
    # Align site types to the sample rows
    site_types_tiled = np.tile(site_types_window, (n_samples, 1))
    color[(is_minor) & (site_types_tiled == 0)] = -1 # Syn
    color[(is_minor) & (site_types_tiled == 1)] = 1  # Nonsyn
    return color

def process_single_csv(file_path, args):
    """Loads a CSV and extracts the window centered on the physical chromosome midpoint."""
    if not os.path.exists(file_path):
        return None

    df = pd.read_csv(file_path)
    site_pos_col = df['site_pos'].values.astype(np.int32)
    raw_genotype_data = df.iloc[:, 2:].values.astype(np.int8) # (Sites, Samples)
    
    site_type_map = {'syn': 0, 'nonsyn': 1}
    site_types = df['site_type'].map(site_type_map).fillna(0).values.astype(np.int8)

    total_sites = raw_genotype_data.shape[0]
    
    if total_sites < args.window_h:
        print(f"Warning: {file_path} only has {total_sites} SNPs. Skipping.")
        return None

    # --- Identify Middle Window based on Physical Position ---
    # Find the SNP index closest to (chrom_len / 2)
    chrom_midpoint = args.chrom_len / 2
    mid_idx = np.argmin(np.abs(site_pos_col - chrom_midpoint))
    
    # Calculate start and end indices for the SNP window
    start_idx = mid_idx - (args.window_h // 2)
    end_idx = start_idx + args.window_h
    
    # Boundary checks: ensure the window doesn't slide off the array edges
    if start_idx < 0:
        start_idx = 0
        end_idx = args.window_h
    elif end_idx > total_sites:
        end_idx = total_sites
        start_idx = end_idx - args.window_h

    # Slice the data (Samples, Sites)
    current_window_raw = raw_genotype_data[start_idx:end_idx, :].T 

    # --- Filter & Downsample Samples ---
    missing_mask = (current_window_raw == args.raw_missing_val)
    missing_rate_per_sample = missing_mask.mean(axis=1)
    
    good_samples_mask = missing_rate_per_sample <= (1 - args.complete_threshold)
    good_indices = np.where(good_samples_mask)[0]

    if len(good_indices) < args.target_samples:
        print(f"Warning: {file_path} lacks enough good samples. Skipping.")
        return None

    selected_indices = np.random.choice(good_indices, args.target_samples, replace=False)
    if args.sort == "none":
        selected_indices.sort()
        
    downsampled_window = current_window_raw[selected_indices, :]

    ch0 = recode_genotype(downsampled_window)
    current_site_types = site_types[start_idx:end_idx]
    ch1 = build_color_channel(ch0, current_site_types)

    image = np.stack([ch0, ch1], axis=-1) # (Samples, Sites, 2)

    # --- Sorting ---
    if args.sort != "none" and helper_haplotypeSorter is not None:
        image = helper_haplotypeSorter.sort_haplotypes(image, ordering=args.sort)

    return image

def natural_key(string_):
    """Helper for natural sorting (1, 2, 10 instead of 1, 10, 2)"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def main():
    parser = argparse.ArgumentParser(description="Extract SNP window centered on physical genomic position.")
    parser.add_argument("input_dir", help="Folder containing the rep_*.csv files")
    parser.add_argument("output_npy", help="Path for final .npy file")
    parser.add_argument("--chrom_len", type=int, default=50000, help="Total length of simulated chromosome")
    parser.add_argument("--window_h", type=int, default=200, help="Number of SNPs in window")
    parser.add_argument("--target_samples", type=int, default=100, help="Number of haplotypes per image")
    parser.add_argument("--complete_threshold", type=float, default=0.5)
    parser.add_argument("--sort", type=str, choices=["rows_freq", "rows_dist", "none"], default="rows_freq")
    parser.add_argument("--raw_missing_val", type=int, default=-1)

    args = parser.parse_args()
    np.random.seed(2023)

    all_images = []

    print(f"Scanning {args.input_dir} for replicates...")
    
    # 1. Use glob to find all matching files dynamically
    search_pattern = os.path.join(args.input_dir, "rep_*.*_haplotypes.csv")
    file_list = glob.glob(search_pattern)
    
    # 2. Sort them naturally so they follow your 1.1, 1.2... sequence
    file_list.sort(key=natural_key)

    total_found = len(file_list)
    print(f"Found {total_found} files. Starting processing...")

    all_images = []
    
    # 3. Iterate directly over the files found on disk
    with tqdm(total=total_found, desc="Processing Replicates") as pbar:
        for file_path in file_list:
            img = process_single_csv(file_path, args)
            
            if img is not None:
                all_images.append(img)
            
            pbar.update(1)

    if not all_images:
        print("No valid images were processed. Check your directory and file names.")
        return

    final_array = np.array(all_images, dtype=np.int8)
    np.save(args.output_npy, final_array)

    print(f"\nProcessing Complete!")
    print(f"Final shape: {final_array.shape}")
    print(f"Centered on chrom position: {args.chrom_len // 2}")
    print(f"Saved to: {args.output_npy}")

if __name__ == "__main__":
    main()
