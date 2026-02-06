#!/usr/bin/env python3
"""
Batch process SLiM simulation files into numpy arrays.

Output Encoding:
  Channel 1 (Allele State):
    -1 = Major allele
     0 = Missing data
     1 = Minor allele

  Channel 2 (Mutation Type):
    -1 = Synonymous mutation
     0 = Major allele / Missing data
     1 = Non-synonymous mutation
"""

import numpy as np
import subprocess
import glob
import os
import sys
import argparse
from tqdm import tqdm


def main():
    parser = argparse.ArgumentParser(
        description="Batch process SLiM simulations into numpy arrays.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Encoding Scheme:
  Channel 1 (Allele State):  -1=major, 0=missing, 1=minor
  Channel 2 (Mutation Type): -1=synonymous, 0=major/missing, 1=non-synonymous
        """
    )
    parser.add_argument("output_npy", help="Path for the final output .npy file")
    parser.add_argument("input_dir", help="Directory containing SLiM .txt output files")
    parser.add_argument("--num_samps", type=int, default=100, 
                        help="Number of haplotypes to sample per simulation (default: 100)")
    parser.add_argument("--window_size", type=int, default=201, 
                        help="Window size in SNPs (default: 201)")
    parser.add_argument("--channels", type=int, choices=[1, 2], default=2, 
                        help="1=allele state only, 2=allele state + mutation type (default: 2)")
    parser.add_argument("--sort", action="store_true", 
                        help="Sort haplotypes by frequency (rows_freq)")
    parser.add_argument("--helper_path", 
                        default="/u/project/ngarud/Garud_lab/Brendan/Utils/helper_processSLiMsims.py",
                        help="Path to helper script")

    args = parser.parse_args()

    # Find input files
    input_files = sorted(glob.glob(os.path.join(args.input_dir, "*.txt")))
    num_sims = len(input_files)
    
    if num_sims == 0:
        print(f"Error: No .txt files found in {args.input_dir}")
        sys.exit(1)

    # Define output shape
    if args.channels == 2:
        sim_shape = (num_sims, args.num_samps, args.window_size, 2)
    else:
        sim_shape = (num_sims, args.num_samps, args.window_size)

    print(f"Processing {num_sims} simulations")
    print(f"Output shape: {sim_shape}")
    print(f"Output file: {args.output_npy}")
    print(f"Channels: {args.channels} ({'allele state + mutation type' if args.channels == 2 else 'allele state only'})")
    print(f"Sorting: {'rows_freq' if args.sort else 'none'}")
    print()

    # Pre-allocate output array with int8 (can store -1, 0, 1)
    print(f"Allocating {args.output_npy}...")
    big_array = np.lib.format.open_memmap(
        args.output_npy, dtype=np.int8, mode="w+", shape=sim_shape
    )
    del big_array  # Release to allow worker processes to write

    # Process each simulation
    sort_flag = "rows_freq" if args.sort else "none"
    failed = []
    
    for i, infile in enumerate(tqdm(input_files, desc="Processing")):
        try:
            result = subprocess.run([
                "python", args.helper_path,
                infile,                     # 1: input SLiM file
                args.output_npy,            # 2: output npy path
                str(args.num_samps),        # 3: number of samples
                str(args.window_size),      # 4: window size
                str(i),                     # 5: index in output array
                str(args.channels),         # 6: number of channels
                sort_flag                   # 7: sorting method
            ], check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            failed.append((infile, e.stderr))
            tqdm.write(f"Failed: {os.path.basename(infile)}")

    # Report results
    print()
    if failed:
        print(f"Completed with {len(failed)} failures:")
        for fname, err in failed[:5]:  # Show first 5 errors
            print(f"  {os.path.basename(fname)}: {err[:100]}...")
        if len(failed) > 5:
            print(f"  ... and {len(failed) - 5} more")
    else:
        print("Success! All simulations processed.")
    
    # Verify output
    print()
    print("Verifying output...")
    data = np.load(args.output_npy, mmap_mode='r')
    print(f"  Shape: {data.shape}")
    print(f"  Dtype: {data.dtype}")
    if args.channels == 2:
        print(f"  Channel 1 unique values: {np.unique(data[..., 0])}")
        print(f"  Channel 2 unique values: {np.unique(data[..., 1])}")
    else:
        print(f"  Unique values: {np.unique(data)}")


if __name__ == "__main__":
    main()