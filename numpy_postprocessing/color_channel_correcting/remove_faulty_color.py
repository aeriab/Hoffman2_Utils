import numpy as np
import argparse
import os
from tqdm import tqdm

def clean_labels(input_path, output_path):
    print(f"Loading input: {input_path}")
    # Load input as read-only
    input_data = np.load(input_path, mmap_mode='r')
    
    print(f"Allocating output: {output_path}")
    # Pre-allocate output memmap with the same shape/dtype as input
    output_data = np.lib.format.open_memmap(
        output_path, 
        mode='w+', 
        dtype=input_data.dtype, 
        shape=input_data.shape
    )

    # Process in chunks of simulations to keep progress tracking accurate 
    # and handle potentially massive files efficiently.
    num_sims = input_data.shape[0]
    print(f"Processing {num_sims} simulations...")

    for i in tqdm(range(num_sims), desc="Cleaning"):
        # Copy the simulation data to the output
        sim_data = np.array(input_data[i])
        
        # Logic: Channel 0 == -1 (Major) but Channel 1 != 0 (Ghost label)
        ghost_labels = (sim_data[..., 0] == -1) & (sim_data[..., 1] != 0)
        
        # Also clean missing data sites (Ch0 == 0) just in case
        missing_labels = (sim_data[..., 0] == 0) & (sim_data[..., 1] != 0)
        
        # Apply the fix
        sim_data[..., 1][ghost_labels] = 0
        sim_data[..., 1][missing_labels] = 0
        
        # Write back to the output memmap
        output_data[i] = sim_data

    # Ensure everything is written to disk
    output_data.flush()
    del input_data
    del output_data
    print(f"\nSuccess! Cleaned file saved to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove mutation labels from major allele sites.")
    parser.add_argument("input", help="Path to the source .npy file")
    parser.add_argument("output", help="Path for the cleaned .npy file")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: {args.input} not found.")
    elif os.path.exists(args.output):
        confirm = input(f"Warning: {args.output} already exists. Overwrite? (y/n): ")
        if confirm.lower() == 'y':
            clean_labels(args.input, args.output)
    else:
        clean_labels(args.input, args.output)