import pandas as pd
import sys
import os

def remove_invariant_sites(input_file, output_file):
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    df = pd.read_csv(input_file)
    
    # Isolate sample columns (everything after site_pos and site_type)
    allele_cols = df.iloc[:, 2:]

    # Keep rows where at least one minor allele (1.0) exists
    is_variant = (allele_cols == 1.0).any(axis=1)
    filtered_df = df[is_variant]

    filtered_df.to_csv(output_file, index=False)
    
    print(f"Processed {len(df)} sites. Retained {len(filtered_df)} variant sites.")

if __name__ == "__main__":
    # Check if arguments are provided; otherwise use defaults
    in_file = sys.argv[1] if len(sys.argv) > 1 else 'input.csv'
    out_file = sys.argv[2] if len(sys.argv) > 2 else 'output_variants.csv'
    
    remove_invariant_sites(in_file, out_file)