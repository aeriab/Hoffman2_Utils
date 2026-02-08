#!/usr/bin/env python3
"""
Filter out bad images from SLiM-processed numpy arrays.

Bad images are all-zero slices that result from failed simulations
(e.g., insufficient segregating sites after downsampling/filtering).

Usage:
    python filter_bad_images.py input.npy output.npy
    python filter_bad_images.py input.npy output.npy --verbose
    python filter_bad_images.py hard.npy hard_clean.npy soft.npy soft_clean.npy neutral.npy neutral_clean.npy
"""

import numpy as np
import argparse
import sys


def filter_bad_images(input_path, output_path, verbose=False):
    """
    Remove all-zero images from a numpy array.
    
    Parameters:
        input_path: Path to input .npy file (shape: [N, samples, snps, channels] or [N, samples, snps])
        output_path: Path to save filtered .npy file
        verbose: Print per-image diagnostics
    
    Returns:
        (n_original, n_kept, n_removed)
    """
    data = np.load(input_path, mmap_mode='r')
    n_original = data.shape[0]
    
    # An image is "bad" if every value in that slice is 0
    # Reshape each image to a flat vector and check if any value is nonzero
    reshaped = data.reshape(n_original, -1)
    good_mask = np.any(reshaped != 0, axis=1)
    
    n_kept = int(np.sum(good_mask))
    n_removed = n_original - n_kept
    
    if verbose:
        bad_indices = np.where(~good_mask)[0]
        if len(bad_indices) > 0:
            print(f"  Bad image indices: {bad_indices.tolist()}" 
                  if len(bad_indices) <= 20 
                  else f"  First 20 bad indices: {bad_indices[:20].tolist()}... ({len(bad_indices)} total)")
    
    # Save filtered array
    filtered = data[good_mask]
    np.save(output_path, filtered)
    
    return n_original, n_kept, n_removed


def main():
    parser = argparse.ArgumentParser(
        description="Remove bad (all-zero) images from SLiM-processed numpy arrays.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file
  python filter_bad_images.py hard_sorted_color.npy hard_clean.npy

  # Multiple files (pairs of input output)
  python filter_bad_images.py \\
      hard.npy hard_clean.npy \\
      soft.npy soft_clean.npy \\
      neutral.npy neutral_clean.npy

  # Verbose mode
  python filter_bad_images.py hard.npy hard_clean.npy --verbose
        """
    )
    parser.add_argument("files", nargs='+', 
                        help="Pairs of input/output .npy files: in1 out1 [in2 out2 ...]")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Print indices of bad images")
    
    args = parser.parse_args()
    
    if len(args.files) % 2 != 0:
        print("Error: Must provide pairs of input/output files.")
        sys.exit(1)
    
    pairs = [(args.files[i], args.files[i+1]) for i in range(0, len(args.files), 2)]
    
    for input_path, output_path in pairs:
        print(f"\n{'='*60}")
        print(f"Input:  {input_path}")
        print(f"Output: {output_path}")
        
        n_orig, n_kept, n_removed = filter_bad_images(input_path, output_path, args.verbose)
        
        print(f"  Original:  {n_orig} images")
        print(f"  Kept:      {n_kept} images")
        print(f"  Removed:   {n_removed} images ({100*n_removed/n_orig:.1f}%)")
        
        # Quick verification
        out = np.load(output_path, mmap_mode='r')
        print(f"  Output shape: {out.shape}, dtype: {out.dtype}")
    
    # Summary for multi-file runs
    if len(pairs) > 1:
        print(f"\n{'='*60}")
        print("Tip: If using for CNN training, make sure to truncate all")
        print("classes to the same count for balanced training.")


if __name__ == "__main__":
    main()