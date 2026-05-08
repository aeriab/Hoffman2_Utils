import numpy as np
import pandas as pd
import sys
import argparse
from tqdm import tqdm

# =============================================================================
# NEW ENCODING SCHEME
# =============================================================================
# Channel 0 (Genotype): -1=major, 0=missing, 1=minor
# Channel 1 (Mutation): -1=syn, 0=major/missing, 1=nonsyn
# =============================================================================

try:
    import helper_haplotypeSorter
except ImportError:
    helper_haplotypeSorter = None

def recode_genotype(raw_block):
    """Recode: 0=major -> -1, 1=minor -> 1, -1=missing -> 0"""
    recoded = np.zeros_like(raw_block, dtype=np.int8)
    recoded[raw_block == 0] = -1
    recoded[raw_block == 1] = 1
    recoded[raw_block == -1] = 0
    return recoded

def build_color_channel(recoded_geno, site_types_window):
    n_samples, n_sites = recoded_geno.shape
    color = np.zeros((n_samples, n_sites), dtype=np.int8)
    is_minor = (recoded_geno == 1)
    site_types_tiled = np.tile(site_types_window, (n_samples, 1))
    color[(is_minor) & (site_types_tiled == 0)] = -1 # Syn
    color[(is_minor) & (site_types_tiled == 1)] = 1  # Nonsyn
    return color

def main():
    # Set seed for reproducibility as requested
    np.random.seed(2023)
    
    parser = argparse.ArgumentParser(description="Convert HMP CSV to windowed NumPy arrays.")
    parser.add_argument("input_csv", help="Path to input CSV")
    parser.add_argument("output_npy", help="Path for output .npy file")
    parser.add_argument("--window_h", type=int, default=200, help="Window size (sites)")
    parser.add_argument("--slide_step", type=int, default=1, help="Step size")
    parser.add_argument("--target_samples", type=int, default=100, help="downsample_n")
    parser.add_argument("--snp_threshold", type=int, default=10, help="Min SNPs in window (S)")
    parser.add_argument("--complete_threshold", type=float, default=0.5, 
                        help="complete_data_threshold (fraction of non-missing required)")
    parser.add_argument("--sort", type=str, choices=["rows_freq", "rows_dist", "none"], default="rows_freq")
    parser.add_argument("--raw_missing_val", type=int, default=-1)

    args = parser.parse_args()
    output_prefix = args.output_npy.rsplit('.', 1)[0]

    # --- Load Data ---
    print(f"Loading data from {args.input_csv}...")
    df = pd.read_csv(args.input_csv)
    
    site_pos_col = df['site_pos'].values.astype(np.int32)
    raw_genotype_data = df.iloc[:, 2:].values.astype(np.int8) # (Sites, Samples)
    total_sites, total_samples = raw_genotype_data.shape

    site_type_map = {'syn': 0, 'nonsyn': 1}
    site_types = df['site_type'].map(site_type_map).fillna(0).values.astype(np.int8)

    # --- Sliding Window Setup ---
    num_possible_windows = int(np.floor((total_sites - args.window_h) / args.slide_step) + 1)
    
    # We collect valid windows in a list because we might throw some out
    valid_windows_list = []
    valid_site_indices_list = []

    print(f"Processing up to {num_possible_windows} windows...")

    for i in tqdm(range(num_possible_windows), desc="Filtering & Downsampling"):
        start_idx = i * args.slide_step
        end_idx = start_idx + args.window_h

        # 1. Check if window has minimum SNPs (S >= snp_threshold)
        # S is typically defined as the number of variant sites in the window
        # In this context, it's the number of rows in our current slice.
        # Since windows are fixed height, we check if the actual data exists.
        current_window_raw = raw_genotype_data[start_idx:end_idx, :].T # (Samples, Sites)
        
        # Calculate S: number of sites that aren't 100% missing or monomorphic
        # (Standard check: if the window is smaller than window_h, skip)
        if current_window_raw.shape[1] < args.window_h:
            continue
            
        # 2. Identify "Good Samples"
        # Missingness is defined as raw_missing_val (-1)
        missing_mask = (current_window_raw == args.raw_missing_val)
        missing_rate_per_sample = missing_mask.mean(axis=1)
        
        # good_hap_bool = haps_missing_data <= (1 - complete_data_threshold)
        good_samples_mask = missing_rate_per_sample <= (1 - args.complete_threshold)
        good_indices = np.where(good_samples_mask)[0]
        num_good_samples = len(good_indices)

        # 3. Check if downsample_n > good_n
        if args.target_samples > num_good_samples:
            # Skip this window (equivalent to throwing it out)
            continue

        # 4. Randomly choose [target_samples] from good_samples without replacement
        selected_indices = np.random.choice(good_indices, args.target_samples, replace=False)
        # Keep original relative order if sorting is 'none', otherwise sorting happens later
        if args.sort == "none":
            selected_indices.sort()
            
        downsampled_window = current_window_raw[selected_indices, :]

        # 5. Process Channels
        ch0 = recode_genotype(downsampled_window)
        current_site_types = site_types[start_idx:end_idx]
        ch1 = build_color_channel(ch0, current_site_types)

        current_image = np.stack([ch0, ch1], axis=-1) # (Samples, Sites, 2)

        # 6. Sorting
        if args.sort != "none" and helper_haplotypeSorter is not None:
            helper_haplotypeSorter.sort_haplotypes(current_image, ordering=args.sort)

        valid_windows_list.append(current_image)
        valid_site_indices_list.append(site_pos_col[start_idx:end_idx])

    # --- Save ---
    if not valid_windows_list:
        print("Error: No windows passed the SNP and missing data thresholds.")
        sys.exit(1)

    final_data = np.array(valid_windows_list, dtype=np.int8)
    final_site_indices = np.array(valid_site_indices_list, dtype=np.int32)

    np.save(args.output_npy, final_data)
    np.save(f"{output_prefix}_map.npy", final_site_indices)

    print(f"\nSuccess! Kept {len(final_data)} windows.")
    print(f"Final shape: {final_data.shape}")
    print(f"Saved to {args.output_npy} and {output_prefix}_map.npy")

if __name__ == "__main__":
    main()