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
    np.random.seed(2023)
    
    parser = argparse.ArgumentParser(description="Convert HMP CSV to windowed NumPy arrays with Invariant Site Removal.")
    parser.add_argument("input_csv", help="Path to input CSV")
    parser.add_argument("output_npy", help="Path for output .npy file")
    parser.add_argument("--window_h", type=int, default=200, help="Window size (sites)")
    parser.add_argument("--slide_step", type=int, default=1, help="Step size")
    parser.add_argument("--target_samples", type=int, default=100, help="downsample_n")
    parser.add_argument("--snp_threshold", type=int, default=10, help="Min SNPs in window after filtering")
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
    
    valid_windows_list = []
    valid_site_indices_list = []

    print(f"Processing up to {num_possible_windows} windows...")

    for i in tqdm(range(num_possible_windows), desc="Filtering & Downsampling"):
        start_idx = i * args.slide_step
        end_idx = start_idx + args.window_h

        # 1. Slice raw window
        # We slice the genotype data and the metadata (types/positions) simultaneously
        current_window_raw = raw_genotype_data[start_idx:end_idx, :] # (Sites, Samples)
        current_site_types = site_types[start_idx:end_idx]
        current_site_pos = site_pos_col[start_idx:end_idx]

        if current_window_raw.shape[0] < args.window_h:
            continue
            
        # 2. Identify "Good Samples" (Pre-downsampling check)
        missing_mask = (current_window_raw == args.raw_missing_val)
        missing_rate_per_sample = missing_mask.mean(axis=0) # Mean across sites
        
        good_samples_mask = missing_rate_per_sample <= (1 - args.complete_threshold)
        good_indices = np.where(good_samples_mask)[0]
        num_good_samples = len(good_indices)

        # 3. Downsampling logic
        if args.target_samples > num_good_samples:
            continue

        selected_indices = np.random.choice(good_indices, args.target_samples, replace=False)
        if args.sort == "none":
            selected_indices.sort()
            
        # Create the sample-subset slice (Samples, Sites)
        downsampled_window = current_window_raw[:, selected_indices].T 

        # 4. REMOVE INVARIANT SITES
        # A site is invariant if there is no minor allele (1) present among the selected samples.
        # We check along the sample axis (axis 0).
        is_variant = np.any(downsampled_window == 1, axis=0)
        
        # Filter genotypes and metadata
        downsampled_window = downsampled_window[:, is_variant]
        filtered_site_types = current_site_types[is_variant]
        filtered_site_pos = current_site_pos[is_variant]

        # 5. POST-FILTER THRESHOLD CHECK
        # Now we check if the window still meets our size requirements.
        # If your model requires EXACTLY window_h sites, use: if downsampled_window.shape[1] != args.window_h:
        # If you just need a minimum number of SNPs (S), use args.snp_threshold:
        if downsampled_window.shape[1] < args.snp_threshold:
            continue

        # (Optional) If you need the window to be padded back to window_h or if you 
        # only want windows that stayed at window_h after invariant removal:
        if downsampled_window.shape[1] < args.window_h:
             # Logic choice: skip windows that lost sites, or keep them as variable-length?
             # Most CNNs need fixed width. If so, uncomment the 'continue' below:
             # continue 
             pass

        # 6. Process Channels
        ch0 = recode_genotype(downsampled_window)
        ch1 = build_color_channel(ch0, filtered_site_types)

        current_image = np.stack([ch0, ch1], axis=-1) # (Samples, Sites, 2)

        # 7. Sorting
        if args.sort != "none" and helper_haplotypeSorter is not None:
            helper_haplotypeSorter.sort_haplotypes(current_image, ordering=args.sort)

        valid_windows_list.append(current_image)
        valid_site_indices_list.append(filtered_site_pos)

    # --- Save ---
    if not valid_windows_list:
        print("Error: No windows passed the SNP and missing data thresholds.")
        sys.exit(1)

    # Note: If windows have different lengths, np.array() will fail or create an object array.
    # We ensure they are consistent or warn the user.
    try:
        final_data = np.array(valid_windows_list, dtype=np.int8)
        final_site_indices = np.array(valid_site_indices_list, dtype=np.int32)
        np.save(args.output_npy, final_data)
        np.save(f"{output_prefix}_map.npy", final_site_indices)
        print(f"\nSuccess! Kept {len(final_data)} windows.")
        print(f"Final shape: {final_data.shape}")
    except ValueError:
        print("\nWarning: Windows have variable widths after invariant site removal.")
        print("Saving as object array (pickled). Models may require fixed-width padding.")
        np.save(args.output_npy, np.array(valid_windows_list, dtype=object))
        np.save(f"{output_prefix}_map.npy", np.array(valid_site_indices_list, dtype=object))

if __name__ == "__main__":
    main()