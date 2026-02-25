import numpy as np
import argparse
import os

def clean_labels(file_path):
    print(f"Opening {file_path} for in-place cleaning...")
    
    # Load with r+ mode to allow direct modification of the file on disk
    data = np.load(file_path, mmap_mode='r+')
    
    # Logic: 
    # If Channel 0 == -1 (Major) AND Channel 1 != 0 (Has a label)
    # Then set Channel 1 to 0.
    
    # We define the mask:
    # Channel 0 is at index 0, Channel 1 is at index 1
    major_allele_mask = (data[..., 0] == -1)
    has_mutation_label = (data[..., 1] != 0)
    
    # Combine conditions: Major allele sites that incorrectly have mutation info
    ghost_labels = major_allele_mask & has_mutation_label
    
    num_fixes = np.sum(ghost_labels)
    
    if num_fixes > 0:
        print(f"Found {num_fixes} ghost labels. Cleaning...")
        # Apply the fix directly to Channel 1
        data[..., 1][ghost_labels] = 0
        print("Cleaning complete.")
    else:
        print("No ghost labels found. Your file is already clean!")

    # Explicitly flush changes to disk and close
    data.flush()
    del data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove mutation labels from major allele sites.")
    parser.add_argument("input", help="Path to the .npy file to clean")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: {args.input} not found.")
    else:
        clean_labels(args.input)