import matplotlib.pyplot as plt
import pandas as pd
import argparse
import numpy as np
import sys
from matplotlib.lines import Line2D

# ─────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Overlay Brendan's CNN scan with Ricky's LD-based results."
)
parser.add_argument('predictions', type=str,
                    help="Path to CNN predictions CSV (Image_Index, P_Neutral, P_Hard, P_Soft, ...)")
parser.add_argument('ricky', type=str,
                    help="Path to Ricky's results CSV (must contain 'site_pos' column)")

# Position mapping for Brendan's windows (same options as original script)
pos_group = parser.add_mutually_exclusive_group()
pos_group.add_argument('--positions', type=str, default=None,
                       help="CSV mapping Image_Index -> Center genomic position")
pos_group.add_argument('--site_indices', type=str, default=None,
                       help="Path to _site_indices.npy from HMP_csv_to_numpy.py")

# Ricky's column to plot on the right axis
parser.add_argument('--ricky_col', type=str, default='iLDS_pval',
                    help="Column from Ricky's file to plot as -log10(p) on right axis "
                         "(default: iLDS_pval). Other options: rNrS_pval, r2G_pval, iLDS, rNrS ...")
parser.add_argument('--ricky_pval', action='store_true',
                    help="Treat --ricky_col as a p-value and plot -log10(value). "
                         "If not set, the raw value is plotted on the right axis.")
parser.add_argument('--ricky_contig', type=str, default=None,
                    help="Filter Ricky's data to a single contig (e.g. FP929051). "
                         "If omitted, all rows are used.")

parser.add_argument('--bin_size', type=int, default=3,
                    help="Rolling window size for CNN predictions (default: 3)")
parser.add_argument('--title', type=str, default='CNN vs LD Scan Comparison')
parser.add_argument('--output', type=str, default='comparison_scan.png')
parser.add_argument('--dpi', type=int, default=300)
args = parser.parse_args()

# ─────────────────────────────────────────────
# Load Brendan's predictions
# ─────────────────────────────────────────────
df = pd.read_csv(args.predictions)
expected = {'Image_Index', 'P_Neutral', 'P_Hard', 'P_Soft'}
if not expected.issubset(df.columns):
    print(f"ERROR: Predictions file missing columns. Need: {expected}, got: {set(df.columns)}")
    sys.exit(1)
print(f"Loaded {len(df)} CNN predictions.")

# ── Assign genomic positions ──
if args.site_indices:
    site_map = np.load(args.site_indices, allow_pickle=True)
    n = min(len(site_map), len(df))
    site_map, df = site_map[:n], df.iloc[:n].copy()
    df['Center'] = [(np.array(w).min() + np.array(w).max()) / 2.0 for w in site_map]
    x_label = 'Genomic Position (Site Index)'
elif args.positions:
    pos_df = pd.read_csv(args.positions, sep=None, engine='python')
    if not {'Image_Index', 'Center'}.issubset(pos_df.columns):
        print("ERROR: --positions file must have 'Image_Index' and 'Center' columns.")
        sys.exit(1)
    df = df.merge(pos_df[['Image_Index', 'Center']], on='Image_Index', how='left')
    df = df.dropna(subset=['Center'])
    x_label = 'Genomic Position (Center BP)'
else:
    df['Center'] = df['Image_Index']
    x_label = 'Window Index'
    print("No position mapping provided — using Image_Index as x-axis.")

df = df.sort_values('Center').reset_index(drop=True)

# ── Rolling average ──
df_binned = df.rolling(window=args.bin_size, center=True, min_periods=1).mean(numeric_only=True)
print(f"Rolling average (bin_size={args.bin_size}) → {len(df_binned)} positions")

# ── Color logic (same as original) ──
highest_prob  = df_binned[['P_Neutral', 'P_Soft', 'P_Hard']].idxmax(axis=1)
final_labels  = np.where(df_binned['P_Neutral'] > 0.01, 'P_Neutral', highest_prob)
color_map     = {'P_Neutral': 'grey', 'P_Soft': 'skyblue', 'P_Hard': 'red'}
point_colors  = pd.Series(final_labels).map(color_map)
y_brendan     = -np.log10(df_binned['P_Neutral'].clip(lower=1e-10))

# ─────────────────────────────────────────────
# Load Ricky's results
# ─────────────────────────────────────────────
rdf = pd.read_csv(args.ricky)
if 'site_pos' not in rdf.columns:
    print(f"ERROR: Ricky's file must contain a 'site_pos' column. Found: {list(rdf.columns)}")
    sys.exit(1)
if args.ricky_col not in rdf.columns:
    print(f"ERROR: --ricky_col '{args.ricky_col}' not found in Ricky's file.")
    print(f"  Available columns: {list(rdf.columns)}")
    sys.exit(1)

if args.ricky_contig:
    rdf = rdf[rdf['contig'] == args.ricky_contig].copy()
    print(f"Filtered Ricky's data to contig '{args.ricky_contig}': {len(rdf)} rows")
else:
    print(f"Loaded {len(rdf)} rows from Ricky's results.")

rdf = rdf.dropna(subset=['site_pos', args.ricky_col]).sort_values('site_pos').reset_index(drop=True)

if args.ricky_pval:
    y_ricky     = -np.log10(rdf[args.ricky_col].clip(lower=1e-10))
    ricky_label = f'-log10({args.ricky_col})'
else:
    y_ricky     = rdf[args.ricky_col]
    ricky_label = args.ricky_col

print(f"Ricky's y-axis: '{ricky_label}', range [{y_ricky.min():.3f}, {y_ricky.max():.3f}]")

# ─────────────────────────────────────────────
# Plot — shared x, two y-axes
# ─────────────────────────────────────────────
fig, ax1 = plt.subplots(figsize=(14, 6))

# Left axis: Brendan's CNN signal
sc = ax1.scatter(df_binned['Center'], y_brendan,
                 c=point_colors, s=15, alpha=0.8, zorder=3)
ax1.set_ylabel('-log10(P_Neutral)  [CNN]', color='black')
ax1.set_ylim(0, 10)
ax1.tick_params(axis='y', labelcolor='black')
ax1.set_xlabel(x_label)
ax1.grid(True, alpha=0.2)

# Right axis: Ricky's statistic
ax2 = ax1.twinx()
ax2.plot(rdf['site_pos'], y_ricky,
         color='darkorange', linewidth=1.2, alpha=0.85,
         label=f"Ricky: {ricky_label}", zorder=2)
ax2.scatter(rdf['site_pos'], y_ricky,
            color='darkorange', s=18, alpha=0.6, zorder=4)
ax2.set_ylabel(ricky_label + '  [Ricky]', color='darkorange')
ax2.tick_params(axis='y', labelcolor='darkorange')

# Optional: shade Ricky's "significant" sites if a 'significance' column exists
if 'significance' in rdf.columns:
    sig_sites = rdf[rdf['significance'] == True]
    if not sig_sites.empty:
        for xpos in sig_sites['site_pos']:
            ax1.axvline(xpos, color='darkorange', alpha=0.15, linewidth=1, zorder=1)
        print(f"Shaded {len(sig_sites)} significant sites from Ricky's data.")

ax1.set_title(f"{args.title}  ({args.bin_size}-window CNN average)")

# Legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',       markersize=8, label='CNN: Neutral'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red',        markersize=8, label='CNN: Hard Sweep'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='skyblue',    markersize=8, label='CNN: Soft Sweep'),
    Line2D([0], [0], color='darkorange', linewidth=2,                                   label=f'Ricky: {ricky_label}'),
]
ax1.legend(handles=legend_elements, loc='upper left', framealpha=0.8)

plt.tight_layout()
plt.savefig(args.output, dpi=args.dpi)
print(f"Plot saved to {args.output}")