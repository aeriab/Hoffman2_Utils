import numpy as np
import pandas as pd
import sys
import argparse
from tqdm import tqdm

# =============================================================================
# ENCODING SCHEME
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
    # Tile site types across all samples for vectorized masking
    site_types_tiled = np.tile(site_types_window, (n_samples, 1))
    color[(is_minor) & (site_types_tiled == 0)] = -1 # Syn
    color[(is_minor) & (site_types_tiled == 1)] = 1  # Nonsyn
    return color

def main():
    np.random.seed(2023)
    
    parser = argparse.ArgumentParser(description="Convert HMP CSV to fixed-width NumPy arrays (Peer-matching logic).")
    parser.add_argument("input_csv", help="Path to input CSV")
    parser.add_argument("output_npy", help="Path for output .npy file")
    parser.add_argument("--window_h", type=int, default=354, help="Window size (sites)")
    parser.add_argument("--slide_step", type=int, default=1, help="Step size")
    parser.add_argument("--target_samples", type=int, default=120, help="downsample_n")
    parser.add_argument("--snp_threshold", type=int, default=120, help="Min SNPs in window (pre-filter)")
    parser.add_argument("--complete_threshold", type=float, default=0.7, 
                        help="complete_data_threshold (fraction of non-missing required)")
    parser.add_argument("--sort", type=str, choices=["rows_freq", "rows_dist", "none"], default="rows_freq")
    parser.add_argument("--raw_missing_val", type=int, default=-1)

    args = parser.parse_args()
    output_prefix = args.output_npy.rsplit('.', 1)[0]

    # --- Load Data ---
    print(f"Loading data from {args.input_csv}...")
    df = pd.read_csv(args.input_csv)
    
    site_pos_col = df['site_pos'].values.astype(np.int32)
    # Genotype data starts at 3rd column: (Sites, Samples)
    raw_genotype_data = df.iloc[:, 2:].values.astype(np.int8) 
    total_sites, total_samples = raw_genotype_data.shape

    site_type_map = {'syn': 0, 'nonsyn': 1}
    site_types = df['site_type'].map(site_type_map).fillna(0).values.astype(np.int8)

    # --- Sliding Window Setup ---
    num_possible_windows = int(np.floor((total_sites - args.window_h) / args.slide_step) + 1)
    
    valid_windows_list = []
    valid_site_indices_list = []

    print(f"Processing up to {num_possible_windows} windows...")

    for i in tqdm(range(num_possible_windows), desc="Processing Windows"):
        start_idx = i * args.slide_step
        end_idx = start_idx + args.window_h

        # 1. Slice raw window
        current_window_raw = raw_genotype_data[start_idx:end_idx, :] # (Sites, Samples)
        current_site_types = site_types[start_idx:end_idx]
        current_site_pos = site_pos_col[start_idx:end_idx]

        # Safety check for end of file
        if current_window_raw.shape[0] < args.window_h:
            continue
            
        # 2. Identify "Good Samples" based on missingness (Peer Logic)
        missing_mask = (current_window_raw == args.raw_missing_val)
        missing_rate_per_sample = missing_mask.mean(axis=0) 
        
        good_samples_mask = missing_rate_per_sample <= (1 - args.complete_threshold)
        good_indices = np.where(good_samples_mask)[0]
        num_good_samples = len(good_indices)

        # 3. Downsampling Check
        # If we don't have enough 'complete' samples, skip this window entirely
        if num_good_samples < args.target_samples:
            continue

        selected_indices = np.random.choice(good_indices, args.target_samples, replace=False)
        if args.sort == "none":
            selected_indices.sort()
            
        # Transpose to (Samples, Sites) for CNN format
        downsampled_window = current_window_raw[:, selected_indices].T 

        # --- IMPORTANT CHANGE: WE NO LONGER REMOVE INVARIANT SITES HERE ---
        # This ensures the window width remains exactly args.window_h (e.g., 354)

        # 4. Process Channels
        ch0 = recode_genotype(downsampled_window)
        ch1 = build_color_channel(ch0, current_site_types)

        current_image = np.stack([ch0, ch1], axis=-1) # Shape: (target_samples, window_h, 2)

        # 5. Sorting (Crucial for CNN patterns)
        if args.sort != "none" and helper_haplotypeSorter is not None:
            helper_haplotypeSorter.sort_haplotypes(current_image, ordering=args.sort)

        valid_windows_list.append(current_image)
        valid_site_indices_list.append(current_site_pos)

    # --- Save ---
    if not valid_windows_list:
        print("Error: No windows passed the sample completeness threshold.")
        sys.exit(1)

    # Because windows are now fixed-width, np.array() will succeed without the 'object' warning
    final_data = np.array(valid_windows_list, dtype=np.int8)
    final_site_indices = np.array(valid_site_indices_list, dtype=np.int32)
    
    np.save(args.output_npy, final_data)
    np.save(f"{output_prefix}_map.npy", final_site_indices)
    
    print(f"\nSuccess! Kept {len(final_data)} windows.")
    print(f"Final shape: {final_data.shape}") # Should be (N, 120, 354, 2)

if __name__ == "__main__":
    main()