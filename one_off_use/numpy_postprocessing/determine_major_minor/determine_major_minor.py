import numpy as np
import argparse
import sys

def determine_maj(image):
    """
    1. Flips -1 and 1 in Channel 0 if 1 is the majority.
    2. Updates Channel 1 so only the new 'minor' alleles are colored.
    """
    # geno: (samples, sites), color: (samples, sites)
    geno = image[:, :, 0]
    color = image[:, :, 1]
    
    num_samples = geno.shape[0]
    num_sites = geno.shape[1]

    for col in range(num_sites):
        col_geno = geno[:, col]
        col_color = color[:, col]
        
        count_ones = np.sum(col_geno == 1)
        
        # Check if 1 is the majority
        if count_ones > (num_samples / 2):
            
            # --- STEP 1: Identify Site Type ---
            # We need to know if this site is 'syn' or 'nonsyn' before we flip.
            # We look for ANY non-zero value in the current color column.
            # In your encoding: -1 = syn, 1 = nonsyn
            site_type = 0
            if np.any(col_color == -1):
                site_type = -1
            elif np.any(col_color == 1):
                site_type = 1

            # --- STEP 2: Flip Genotype (Channel 0) ---
            new_col_geno = col_geno.copy()
            new_col_geno[col_geno == 1] = -1
            new_col_geno[col_geno == -1] = 1
            image[:, col, 0] = new_col_geno

            # --- STEP 3: Update Color (Channel 1) ---
            # Only label the site if the NEW genotype is 1 (the new minor allele)
            new_col_color = np.zeros_like(col_color)
            new_col_color[new_col_geno == 1] = site_type
            image[:, col, 1] = new_col_color
            
    return image

def entire_arr_maj(all_images):
    """Loops through the 4D array and processes each image."""
    for i in range(all_images.shape[0]):
        all_images[i] = determine_maj(all_images[i])
    return all_images

def main():
    parser = argparse.ArgumentParser(description="Flip majority alleles and sync color channels.")
    parser.add_argument("input_npy", help="Path to input .npy")
    parser.add_argument("output_npy", help="Path to output .npy")
    
    args = parser.parse_args()

    print(f"Loading {args.input_npy}...")
    data = np.load(args.input_npy)
    
    print("Processing images and syncing channels...")
    processed_data = entire_arr_maj(data)
    
    print(f"Saving to {args.output_npy}...")
    np.save(args.output_npy, processed_data)
    print("Success!")

if __name__ == "__main__":
    main()