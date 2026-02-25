import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import re

# 1. Load the good samples list
with open('Alistipes_finegoldii_56071_good_samples.txt', 'r') as f:
    good_samples = {line.strip() for line in f if line.strip()}

# 2. Find and sort the processed files
files = glob.glob("A_finegoldii_*.csv")
files.sort(key=lambda x: int(re.search(r'(\d+)', x).group()))

# Dictionary to store how many sites have exactly N samples
# Key: Number of samples, Value: Total number of sites
occupancy_counts = {}

for fname in files:
    df = pd.read_csv(fname, low_memory=False)
    
    # Identify sample columns present in this file that are also in our "good" list
    sample_cols = [c for c in df.columns if c in good_samples]
    if not sample_cols:
        continue
        
    # --- Usability Filter ---
    # Calculate missingness per sample within THIS file
    # A sample is usable if it has <= 10% missing data (NaNs) in this contig
    missing_pct = df[sample_cols].isna().mean()
    usable_samples = missing_pct[missing_pct <= 0.10].index.tolist()
    
    if not usable_samples:
        continue

    # For each site (row), count how many of these "usable" samples have data
    # (notna() returns True if the site is present)
    counts_per_site = df[usable_samples].notna().sum(axis=1)
    
    # Update global counts
    for count in counts_per_site:
        occupancy_counts[count] = occupancy_counts.get(count, 0) + 1

# 3. Prepare data for plotting
x_values = sorted(occupancy_counts.keys())
y_values = [occupancy_counts[x] for x in x_values]

# 4. Generate the Plot
plt.figure(figsize=(10, 6))
plt.bar(x_values, y_values, color='skyblue', edgecolor='navy')

plt.title(r'Site Occupancy for $A. finegoldii$ (Samples with $\leq 10\%$ Missing Data)')
plt.xlabel('Number of Samples')
plt.ylabel('Number of Usable Sites (Images)')
plt.grid(axis='y', linestyle='--', alpha=0.7)

plt.savefig('A_finegoldii_occupancy_plot.png', dpi=300)
print("Plot saved as A_finegoldii_occupancy_plot.png")