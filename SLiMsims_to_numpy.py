import numpy as np
import subprocess
import glob
import os
import sys
import argparse
from tqdm import tqdm

def main():
    parser = argparse.ArgumentParser(description="Batch process SLiM simulations (Direct Write Version).")
    parser.add_argument("output_npy", help="Path for the final output .npy file")
    parser.add_argument("input_dir", help="Directory containing SLiM .txt output files")
    parser.add_argument("--num_samps", type=int, default=100, help="Target haplotypes")
    parser.add_argument("--window_size", type=int, default=201, help="Window size (SNPs)")
    parser.add_argument("--channels", type=int, choices=[1, 2], default=2, help="1=BW, 2=Color")
    parser.add_argument("--sort", action="store_true", help="Apply rows_dist sorting")
    parser.add_argument("--helper_path", default="/u/project/ngarud/Garud_lab/Brendan/Utils/helper_processSLiMsims.py")

    args = parser.parse_args()

    input_files = sorted(glob.glob(os.path.join(args.input_dir, "*.txt")))
    num_sims = len(input_files)
    
    if num_sims == 0:
        print(f"Error: No .txt files in {args.input_dir}"); sys.exit(1)

    # Define Shape & Allocate
    sim_shape = (num_sims, args.num_samps, args.window_size)
    if args.channels == 2: sim_shape += (2,)

    print(f"Allocating {num_sims} sims to {args.output_npy}...")
    
    # Pre-create file as int8 (most efficient for 0/1 genotype data)
    big_array = np.lib.format.open_memmap(
        args.output_npy, dtype=np.int8, mode="w+", shape=sim_shape
    )
    del big_array # Hand off to workers

    # Execute Workers
    for i, infile in enumerate(tqdm(input_files)):
        sort_flag = "rows_dist" if args.sort else "none"
        try:
            subprocess.run([
                "python", args.helper_path,
                infile, args.output_npy,
                str(args.num_samps), str(args.window_size),
                str(i), str(num_sims), str(args.channels), sort_flag
            ], check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(f"\nFailed: {infile}\n{e.stderr.decode()}")

    print("Success!")

if __name__ == "__main__":
    main()