"""
plot_downsample.py

Reads a bacterial SNP CSV (e.g. A_finegoldii_01.csv) and plots:
  x-axis: downsample number (# of high-missing samples we're willing to remove per image)
  y-axis: number of usable images

Images are defined as non-overlapping windows of exactly 90 SNPs, tiled per contig.
Windows with fewer than 90 SNPs (remainders at end of a contig) are discarded.
  rows = samples, columns = SNP positions (always exactly 90)

A sample is "high-missing" within an image if >= 10% of its 90 SNPs are missing (>= 9 sites).
An image is USABLE at downsample level D if it has <= D such high-missing samples.

Usage:
    python plot_downsample.py A_finegoldii_01.csv
    python plot_downsample.py A_finegoldii_01.csv --missing-threshold 0.10 --output plot.png
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

IMAGE_SIZE = 90  # fixed SNP window size


def load_data(csv_path: str):
    df = pd.read_csv(csv_path, low_memory=False)
    # First 4 columns are metadata; the rest are samples
    sample_cols = list(df.columns[4:])
    return df, sample_cols


def build_images(df, sample_cols):
    """
    Tile each contig's SNPs into non-overlapping windows of exactly IMAGE_SIZE.
    Returns list of bool arrays of shape (n_samples, IMAGE_SIZE): True = MISSING.
    """
    images = []
    for contig, grp in df.groupby("contig", sort=False):
        # Convert sample columns to numeric; non-numeric -> NaN (missing)
        snp_matrix = grp[sample_cols].apply(pd.to_numeric, errors="coerce").values
        # snp_matrix shape: (n_snps_in_contig, n_samples)
        n_snps = snp_matrix.shape[0]

        for start in range(n_snps - IMAGE_SIZE + 1):
            window = snp_matrix[start:start + IMAGE_SIZE, :]  # (90, n_samples)
            missing_mask = np.isnan(window).T                  # (n_samples, 90)
            images.append(missing_mask)

    return images


def count_usable_images(images, missing_threshold: float = 0.10):
    """
    For each sample count S in [1, n_samples], count how many images are usable
    when we keep exactly S samples (removing the worst n_samples - S samples).
    An image is usable at sample count S if it has <= (n_samples - S) bad samples,
    i.e. after removing up to (n_samples - S) high-missing samples, all remaining
    samples have < missing_threshold fraction missing.
    """
    if not images:
        raise ValueError("No images found. Check your CSV -- are there enough SNPs per contig?")

    n_samples = images[0].shape[0]

    bad_sample_counts = []
    for mask in images:
        per_sample_missing = mask.mean(axis=1)  # fraction missing per sample (over 90 SNPs)
        n_bad = int((per_sample_missing >= missing_threshold).sum())
        bad_sample_counts.append(n_bad)

    bad_counts = np.array(bad_sample_counts)

    # S = number of samples kept; D = n_samples - S = number removed
    # image usable at S iff bad_count <= D = n_samples - S
    S_values = np.arange(1, n_samples + 1)
    usable_counts = np.array([(bad_counts <= (n_samples - S)).sum() for S in S_values])

    return S_values, usable_counts, len(bad_counts), n_samples


def plot(S_values, usable_counts, total_images, n_samples, missing_threshold, output_path):
    fig, ax = plt.subplots(figsize=(9, 5))

    # Plot the continuous orange line
    ax.plot(S_values, usable_counts, color="darkorange", linewidth=2.5, zorder=2)

    # Black dots only every 5 sample numbers
    dot_mask = (S_values % 5 == 0)
    ax.scatter(S_values[dot_mask], usable_counts[dot_mask],
               color="black", s=40, zorder=3)

    ax.set_xlabel("Sample Number (samples retained per image)", fontsize=12)
    ax.set_ylabel("# Usable Images", fontsize=12)
    ax.set_title(
        f"Usable Images vs. Sample Number\n"
        f"({IMAGE_SIZE}-SNP windows, missing threshold = {missing_threshold:.0%}, "
        f"total images = {total_images})",
        fontsize=13
    )
    ax.set_xlim(left=1, right=n_samples)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    ax.axhline(total_images, color="gray", linestyle="--", linewidth=1, alpha=0.7,
               label=f"Total images ({total_images})")
    ax.legend(fontsize=10)
    ax.grid(True, linestyle="--", alpha=0.4)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    print(f"Plot saved to: {output_path}")
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Plot downsample number vs usable 90-SNP haplotype images."
    )
    parser.add_argument("csv", help="Path to input CSV (e.g. A_finegoldii_01.csv)")
    parser.add_argument("--missing-threshold", type=float, default=0.10,
                        help="Fraction of missing sites to flag a sample (default: 0.10 = 10%%)")
    parser.add_argument("--output", default="downsample_plot.png",
                        help="Output plot filename (default: downsample_plot.png)")
    args = parser.parse_args()

    print(f"Loading {args.csv}...")
    df, sample_cols = load_data(args.csv)
    print(f"  {len(df)} SNP sites, {len(sample_cols)} samples")

    print(f"Building {IMAGE_SIZE}-SNP images (non-overlapping windows per contig)...")
    images = build_images(df, sample_cols)
    print(f"  {len(images)} images found")

    if len(images) == 0:
        print("No images produced. Your contigs may all have fewer than 90 SNPs.")
        return

    S_values, usable_counts, total_images, n_samples = count_usable_images(images, args.missing_threshold)

    print(f"\nSample of results (sample count -> usable images):")
    step = max(1, len(S_values) // 10)
    for i in range(0, len(S_values), step):
        print(f"  S={S_values[i]:3d}: {usable_counts[i]:4d} / {total_images} usable")

    plot(S_values, usable_counts, total_images, n_samples, args.missing_threshold, args.output)


if __name__ == "__main__":
    main()