import numpy as np
import pandas as pd
import sys
import argparse
from tqdm import tqdm

# =============================================================================
# NEW ENCODING SCHEME
# =============================================================================
# Channel 0 ("Black & White" - Genotype):
#   -1 = major allele     (CHARCOAL #2F3640)
#    0 = missing data      (GREY     #BDC3C7)
#    1 = minor allele      (BLUE     #3498DB)
#
# Channel 1 ("Color" - Mutation Type):
#   -1 = synonymous        (YELLOW   #E1B12C)
#    0 = major/missing      (CHARCOAL #2F3640)
#    1 = non-synonymous     (GREEN    #44BD32)
#
# CSV RAW VALUES:  0 = ref/major,  1 = alt/minor,  -1 = missing
# Ch0:   0 -> -1 (major),  -1 -> 0 (missing),  1 -> 1 (minor)
# Ch1:   Only encode syn/nonsyn where sample carries minor allele (ch0==1)
#                  syn site + minor allele  -> -1
#                  nonsyn site + minor allele -> 1
#                  otherwise                  -> 0
# =============================================================================

# How to use:
# Standard (2 channels, sorted by frequency, downsampled to 100):
#   python HMP_csv_to_numpy.py data.csv output.npy --window_h 200 --slide_step 10 --target_samples 100 --sort rows_freq
#
# Sort by distance:
#   python HMP_csv_to_numpy.py data.csv output.npy --sort rows_dist
#
# No sorting:
#   python HMP_csv_to_numpy.py data.csv output.npy --sort none
#
# Minimalist (defaults to rows_freq sorting, full samples):
#   python HMP_csv_to_numpy.py data.csv output.npy

try:
    import helper_haplotypeSorter
except ImportError:
    helper_haplotypeSorter = None


def recode_genotype(raw_block):
    """
    Recode raw CSV genotype values to the new encoding.
    Raw:  0=major, 1=minor, -1=missing
    New: -1=major,  1=minor,  0=missing
    """
    recoded = np.zeros_like(raw_block, dtype=np.int8)
    recoded[raw_block == 0] = -1   # major allele
    recoded[raw_block == 1] = 1    # minor allele
    recoded[raw_block == -1] = 0   # missing data
    return recoded


def build_color_channel(recoded_geno, site_types_window):
    """
    Build the color channel (channel 1).

    Args:
        recoded_geno: (n_samples, n_sites) with new encoding (-1=major, 0=missing, 1=minor)
        site_types_window: (n_sites,) with 0=syn, 1=nonsyn (from site_type column)

    Returns:
        color_channel: (n_samples, n_sites) int8
            -1 = synonymous mutation (minor allele at syn site)
             0 = major allele or missing data
             1 = non-synonymous mutation (minor allele at nonsyn site)
    """
    n_samples, n_sites = recoded_geno.shape
    color = np.zeros((n_samples, n_sites), dtype=np.int8)

    # Mask: only sites where sample carries minor allele
    is_minor = (recoded_geno == 1)

    # Tile site types across samples
    site_types_tiled = np.tile(site_types_window, (n_samples, 1))  # 0=syn, 1=nonsyn

    # Synonymous + minor allele -> -1
    color[(is_minor) & (site_types_tiled == 0)] = -1

    # Non-synonymous + minor allele -> 1
    color[(is_minor) & (site_types_tiled == 1)] = 1

    # Everything else stays 0 (major allele or missing)
    return color


def main():
    np.random.seed(2023)
    parser = argparse.ArgumentParser(
        description="Convert HMP CSV to windowed NumPy arrays with NEW encoding scheme."
    )
    parser.add_argument("input_csv", help="Path to input CSV")
    parser.add_argument("output_npy", help="Path for output .npy file")
    parser.add_argument("--window_h", type=int, default=200, help="Window size (sites)")
    parser.add_argument("--slide_step", type=int, default=1, help="Step size for sliding window")
    parser.add_argument("--target_samples", type=int, default=None,
                        help="Downsample to this many haplotypes (picks ones with least missing data)")
    parser.add_argument("--sort", type=str, choices=["rows_freq", "rows_dist", "none"],
                        default="rows_freq",
                        help="Haplotype sorting method (default: rows_freq)")
    parser.add_argument("--raw_missing_val", type=int, default=-1,
                        help="Value representing missing data in the RAW csv (default: -1)")

    args = parser.parse_args()

    output_prefix = args.output_npy.rsplit('.', 1)[0]
    map_filename = f"{output_prefix}_map.npy"

    # --- Load Data ---
    print(f"Loading data from {args.input_csv}...")
    try:
        df = pd.read_csv(args.input_csv)
    except Exception as e:
        print(f"Error loading file: {e}")
        sys.exit(1)

    site_pos_col = df['site_pos'].values.astype(np.int32)
    raw_genotype_data = df.iloc[:, 2:].values.astype(np.int8)  # (Total_Sites, Total_Samples)
    total_sites, total_samples = raw_genotype_data.shape

    # Site types: syn=0, nonsyn=1
    site_type_map = {'syn': 0, 'nonsyn': 1}
    site_types = df['site_type'].map(site_type_map).fillna(0).values.astype(np.int8)

    # Handle sample count
    final_num_samples = args.target_samples if args.target_samples else total_samples
    if args.target_samples and args.target_samples > total_samples:
        print(f"Warning: Target samples ({args.target_samples}) > available ({total_samples}). Using {total_samples}.")
        final_num_samples = total_samples

    # --- Sliding Window ---
    num_images = int(np.floor((total_sites - args.window_h) / args.slide_step) + 1)
    final_shape = (num_images, final_num_samples, args.window_h, 2)

    print(f"Preparing {num_images} windows. Output shape: {final_shape}")
    print(f"Encoding: Ch0[-1=major, 0=missing, 1=minor] Ch1[-1=syn, 0=major/missing, 1=nonsyn]")
    print(f"Sorting: {args.sort}")

    final_data = np.zeros(final_shape, dtype=np.int8)
    final_site_indices = np.zeros((num_images, args.window_h), dtype=np.int32)

    # --- Processing Loop ---
    for i in tqdm(range(num_images), desc="Processing windows"):
        start_idx = i * args.slide_step
        end_idx = start_idx + args.window_h

        final_site_indices[i] = site_pos_col[start_idx:end_idx]

        # Raw block: (Sites, Samples) -> Transpose to (Samples, Sites)
        raw_window = raw_genotype_data[start_idx:end_idx, :].T

        # --- A. DOWNSAMPLING (Quality Based) ---
        if args.target_samples and total_samples > args.target_samples:
            # Count missing in RAW encoding
            missing_counts = (raw_window == args.raw_missing_val).sum(axis=1)
            best_indices = np.argsort(missing_counts)[:args.target_samples]
            best_indices.sort()
            raw_window = raw_window[best_indices]

        # --- B. RECODE GENOTYPE (Channel 0) ---
        ch0 = recode_genotype(raw_window)

        # --- C. BUILD COLOR CHANNEL (Channel 1) ---
        current_site_types = site_types[start_idx:end_idx]
        ch1 = build_color_channel(ch0, current_site_types)

        # --- D. STACK CHANNELS ---
        current_image = np.stack([ch0, ch1], axis=-1)  # (samples, sites, 2)

        # --- E. SORTING ---
        if args.sort != "none":
            if helper_haplotypeSorter is not None:
                helper_haplotypeSorter.sort_haplotypes(current_image, ordering=args.sort)
            else:
                if i == 0:
                    print("Warning: Sorting requested but helper_haplotypeSorter.py not found.")

        final_data[i] = current_image

    # --- Save ---
    np.save(args.output_npy, final_data)
    np.save(map_filename, final_site_indices)

    print(f"\nSuccess! Data saved to {args.output_npy}")
    print(f"Site position map saved to {map_filename}")
    print(f"Final shape: {final_data.shape}")
    print(f"  Axis 0: {final_data.shape[0]} windows")
    print(f"  Axis 1: {final_data.shape[1]} samples")
    print(f"  Axis 2: {final_data.shape[2]} sites per window")
    print(f"  Axis 3: 2 channels [genotype, mutation_type]")


if __name__ == "__main__":
    main()