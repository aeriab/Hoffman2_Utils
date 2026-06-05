import pandas as pd
import argparse

def convert_csv(input_path, output_path):
    print(f"Loading data from {input_path}...")
    
    # Read the original CSV. Pandas automatically converts empty consecutive commas to NaN
    df = pd.read_csv(input_path)

    # 1. Drop unnecessary columns
    cols_to_drop = [col for col in ['contig', 'gene_id'] if col in df.columns]
    df = df.drop(columns=cols_to_drop, errors='ignore')

    # 2. Format the site_type column
    if 'site_type' in df.columns:
        # Update these dictionary keys to match whatever your raw data uses.
        # e.g., 'NC' to 'syn', and whatever designates non-synonymous to 'nonsyn'
        mapping_dict = {
            'NC': 'syn',
            'C': 'nonsyn', 
            'synonymous': 'syn',
            'non-synonymous': 'nonsyn'
        }
        df['site_type'] = df['site_type'].replace(mapping_dict)

    # 3. Separate metadata from genotype data
    metadata_cols = [col for col in ['site_pos', 'site_type'] if col in df.columns]
    sample_cols = [col for col in df.columns if col not in metadata_cols]

    metadata = df[metadata_cols]
    raw_genotypes = df[sample_cols]

    print("Formatting genotypes (0.0=major, 1.0=minor, NaN -> -1.0)...")
    
    # Since 0.0 and 1.0 are already correct, we just need to handle any missing data.
    # .fillna() is highly optimized in pandas and replaces all NaNs with -1.0.
    recoded_genotypes = raw_genotypes.fillna(-1.0)

    # 4. Combine metadata and recoded genotypes
    final_df = pd.concat([metadata, recoded_genotypes], axis=1)
    
    print(f"Saving formatted data to {output_path}...")
    
    # Save to CSV, ensuring floats are formatted with one decimal place (e.g. -1.0)
    final_df.to_csv(output_path, index=False, float_format='%.1f')
    print("Conversion complete!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format haplotype CSV to target format (-1.0/0.0/1.0).")
    parser.add_argument("input_csv", help="Path to the original raw CSV file")
    parser.add_argument("output_csv", help="Path to save the newly formatted CSV file")
    
    args = parser.parse_args()
    convert_csv(args.input_csv, args.output_csv)