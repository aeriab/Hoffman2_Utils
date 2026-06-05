import matplotlib.pyplot as plt
import pandas as pd
import argparse
import numpy as np
import sys
import os
import glob
import re

parser = argparse.ArgumentParser(description="Plot genomic scan results for multiple contigs side-by-side")
parser.add_argument('pred_dir', type=str, help="Directory containing prediction .txt files")
parser.add_argument('--bin_size', type=int, default=3, help="Number of consecutive windows to average (default: 3)")
parser.add_argument('--title', type=str, default='Genomic Scan (Multiple Contigs)', help="Plot title")
parser.add_argument('--output', type=str, default='multi_genomic_scan.png', help="Output plot filename")
parser.add_argument('--dpi', type=int, default=300, help="Plot resolution (default: 300)")
parser.add_argument('--y_max', type=float, default=10.0, help="Maximum value for the y-axis (default: 10.0)")
args = parser.parse_args()

if not os.path.isdir(args.pred_dir):
    print(f"ERROR: Directory not found: {args.pred_dir}")
    sys.exit(1)

# --- Find and Sort Files ---
# Matches files ending in numbers before .txt (e.g., predictions_01.txt, example_predictions2.txt)
all_files = glob.glob(os.path.join(args.pred_dir, "*.txt"))

def extract_num(filepath):
    basename = os.path.basename(filepath)
    match = re.search(r'(\d+)\.txt$', basename)
    return int(match.group(1)) if match else -1

# Filter out files that don't have a numeric suffix and sort them
valid_files = [f for f in all_files if extract_num(f) != -1]
valid_files.sort(key=extract_num)

if not valid_files:
    print(f"ERROR: No valid prediction files with numerical suffixes found in {args.pred_dir}")
    sys.exit(1)

print(f"Found {len(valid_files)} prediction files. Processing...")

# --- Load, Bin, and Concatenate Data ---
all_binned_dfs = []
contig_boundaries = []
current_x_offset = 0

for filepath in valid_files:
    df = pd.read_csv(filepath)
    
    expected_cols = {'Image_Index', 'P_Neutral', 'P_Hard', 'P_Soft'}
    if not expected_cols.issubset(df.columns):
        print(f"WARNING: Skipping {os.path.basename(filepath)} - missing expected columns.")
        continue

    # Use Image_Index as the base x-axis
    df['Center'] = df['Image_Index']
    df = df.sort_values('Center').reset_index(drop=True)

    # Apply rolling average strictly within this contig
    df_binned = df.rolling(window=args.bin_size, center=True, min_periods=1).mean(numeric_only=True)
    
    # Offset the x-axis so it plots sequentially after the previous contig
    df_binned['Continuous_Center'] = df_binned['Center'] + current_x_offset
    
    # Record boundaries for the alternating background colors
    start_pos = current_x_offset
    end_pos = df_binned['Continuous_Center'].max()
    contig_boundaries.append((start_pos, end_pos))
    
    # Update offset for the next file (add a small gap of 1 unit)
    current_x_offset = end_pos + 1 
    
    all_binned_dfs.append(df_binned)

if not all_binned_dfs:
    print("ERROR: No data loaded.")
    sys.exit(1)

# Combine into one master dataframe
master_df = pd.concat(all_binned_dfs, ignore_index=True)

# --- Assign colors by top class ---
highest_prob = master_df[['P_Neutral', 'P_Soft', 'P_Hard']].idxmax(axis=1)

# P_Neutral > 0.01 corresponds to -log10(P_Neutral) < 2
final_labels = np.where(master_df['P_Neutral'] > 0.01, 'P_Neutral', highest_prob)

color_map = {
    'P_Neutral': 'grey',
    'P_Soft': 'skyblue',
    'P_Hard': 'red'
}
point_colors = pd.Series(final_labels).map(color_map)

# --- Plot ---
fig, ax = plt.subplots(figsize=(16, 6)) # Made slightly wider to accommodate multiple contigs

# Draw alternating background spans
for i, (start, end) in enumerate(contig_boundaries):
    # Alternate between two subtle grey hex colors
    bg_color = '#EAEAEA' if i % 2 == 0 else '#F7F7F7'
    ax.axvspan(start, end, facecolor=bg_color, alpha=1.0, zorder=0, edgecolor='none')

y_values = -np.log10(master_df['P_Neutral'].clip(lower=1e-10))

# Plot scatter points on top of the background (zorder=2)
ax.scatter(master_df['Continuous_Center'], y_values, c=point_colors, s=15, alpha=0.8, zorder=2)

ax.set_title(f'{args.title} ({args.bin_size}-window average)')
ax.set_xlabel('Continuous Window Index across Contigs')
ax.set_ylabel('-log10(P_Neutral)')
ax.grid(True, alpha=0.2, zorder=1)

ax.set_ylim(0, args.y_max)
ax.set_xlim(0, current_x_offset)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='grey', markersize=8, label='Neutral'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label='Hard Sweep'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='skyblue', markersize=8, label='Soft Sweep'),
]
ax.legend(handles=legend_elements, loc='upper right')

plt.tight_layout()
plt.savefig(args.output, dpi=args.dpi)
print(f"Plot saved to {args.output}")