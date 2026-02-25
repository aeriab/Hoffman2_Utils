#!/usr/bin/env python3
"""
Visualize haplotype images from numpy arrays.
Generates PNG files that can be downloaded and viewed locally.

Data Encoding:
  Channel 1 (Allele State):
    -1 = Major allele (Red)
     0 = Missing data (Grey)
     1 = Minor allele (Blue)

  Channel 2 (Mutation Type):
    -1 = Synonymous mutation (Yellow)
     0 = Major allele / Missing (Muted Red)
     1 = Non-synonymous mutation (Green)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


def plot_haplotype_image(data, idx, output_dir, prefix="haplotype"):
    """
    Plot a single haplotype image showing both channels.
    
    Parameters:
        data: numpy array of shape (n_images, n_samples, n_sites, 2)
        idx: index of image to plot
        output_dir: directory to save output
        prefix: filename prefix
    """
    img = data[idx]  # shape: (n_samples, n_sites, 2)
    
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Channel 1: Allele State
    # -1 = major (red), 0 = missing (grey), 1 = minor (blue)
    # colors_ch1 = ['tab:red', 'grey', 'tab:blue']
    colors_ch1 = ['#2F3640', '#BDC3C7', '#3498DB']
    cmap_ch1 = ListedColormap(colors_ch1)
    
    im0 = axes[0].imshow(img[:, :, 0], aspect='auto', cmap=cmap_ch1, 
                          vmin=-1.5, vmax=1.5, interpolation='nearest')
    axes[0].set_title('Channel 1: Allele State')
    axes[0].set_xlabel('Sites')
    axes[0].set_ylabel('Samples')

    legend_ch1 = [
        Patch(facecolor='#2F3640', edgecolor='black', label='-1 (major allele)'),
        Patch(facecolor='#BDC3C7', edgecolor='black', label='0 (missing data)'),
        Patch(facecolor='#3498DB', edgecolor='black', label='1 (minor allele)')
    ]
    axes[0].legend(handles=legend_ch1, loc='upper right', fontsize=9, framealpha=0.9)
    
    # Channel 2: Mutation Type
    # -1 = synonymous (yellow), 0 = major/missing (muted red), 1 = non-syn (green)
    # colors_ch2 = ['gold', 'rosybrown', 'tab:green']
    colors_ch2 = ['#E1B12C', '#2F3640', '#44BD32']
    cmap_ch2 = ListedColormap(colors_ch2)

    im1 = axes[1].imshow(img[:, :, 1], aspect='auto', cmap=cmap_ch2, 
                          vmin=-1.5, vmax=1.5, interpolation='nearest')
    axes[1].set_title('Channel 2: Mutation Type')
    axes[1].set_xlabel('Sites')
    axes[1].set_ylabel('Samples')

    legend_ch2 = [
        Patch(facecolor='#E1B12C', edgecolor='black', label='-1 (synonymous)'),
        Patch(facecolor='#2F3640', edgecolor='black', label='0 (major/missing)'),
        Patch(facecolor='#44BD32', edgecolor='black', label='1 (non-synonymous)')
    ]
    axes[1].legend(handles=legend_ch2, loc='upper right', fontsize=9, framealpha=0.9)
    
    plt.tight_layout()
    outpath = Path(output_dir) / f"{prefix}_{idx:04d}.png"
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    
    return outpath


def plot_grid(data, indices, output_path, title="Haplotype Images"):
    """
    Plot multiple haplotype images in a grid (Channel 1 view).
    Uses the allele state channel for grid visualization.
    """
    n = len(indices)
    cols = min(5, n)
    rows = (n + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(3*cols, 3*rows))
    if n == 1:
        axes = np.array([[axes]])
    axes = axes.reshape(rows, cols) if rows > 1 else axes.reshape(1, -1)
    
    for i, idx in enumerate(indices):
        row, col = i // cols, i % cols
        ax = axes[row, col]
        
        img = data[idx]
        # Create RGB image from Channel 1 (allele state)
        combined = np.zeros((*img.shape[:2], 3))
        
        # # Channel 1 encoding: -1=major (red), 0=missing (grey), 1=minor (blue)
        # combined[img[:, :, 0] == -1] = [0.9, 0.2, 0.2]   # major → red
        # combined[img[:, :, 0] == 0] = [0.5, 0.5, 0.5]    # missing → grey
        # combined[img[:, :, 0] == 1] = [0.2, 0.4, 0.8]    # minor → blue

        # Convert Hex to RGB [0, 1] for the grid view
        combined[img[:, :, 0] == -1] = [0.17, 0.24, 0.31]   # #2C3E50 (Major)
        combined[img[:, :, 0] == 0]  = [0.74, 0.76, 0.78]   # #BDC3C7 (Missing)
        combined[img[:, :, 0] == 1]  = [0.20, 0.60, 0.86]   # #3498DB (Minor)
        
        ax.imshow(combined, aspect='auto', interpolation='nearest')
        ax.set_title(f'Image {idx}', fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    
    # Hide empty subplots
    for i in range(n, rows * cols):
        row, col = i // cols, i % cols
        axes[row, col].axis('off')
    
    plt.suptitle(title, fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    return output_path


def validate_encoding(data, n_samples=100):
    """
    Validate that the encoding is consistent.
    
    Expected encoding:
      Ch1=-1, Ch2=0  → Major allele (no mutation)
      Ch1=0,  Ch2=0  → Missing data
      Ch1=1,  Ch2=-1 → Minor allele, synonymous
      Ch1=1,  Ch2=1  → Minor allele, non-synonymous
    """
    indices = np.random.choice(data.shape[0], min(n_samples, data.shape[0]), replace=False)
    
    print("Encoding Validation:")
    print(f"  Shape: {data.shape}")
    print(f"  Dtype: {data.dtype}")
    print(f"  Channel 1 unique values: {np.unique(data[..., 0])}")
    print(f"  Channel 2 unique values: {np.unique(data[..., 1])}")
    
    ch1 = data[indices, :, :, 0]
    ch2 = data[indices, :, :, 1]
    
    # Check: Where ch1=-1 (major), ch2 should be 0
    major_sites = ch1 == -1
    ch2_at_major = ch2[major_sites]
    invalid_major = np.sum(ch2_at_major != 0)
    print(f"  Major allele sites with ch2≠0: {invalid_major} (should be 0)")
    
    # Check: Where ch1=0 (missing), ch2 should be 0
    missing_sites = ch1 == 0
    ch2_at_missing = ch2[missing_sites]
    invalid_missing = np.sum(ch2_at_missing != 0)
    print(f"  Missing data sites with ch2≠0: {invalid_missing} (should be 0)")
    
    # Check: Where ch1=1 (minor), ch2 should be -1 or 1
    minor_sites = ch1 == 1
    ch2_at_minor = ch2[minor_sites]
    invalid_minor = np.sum(ch2_at_minor == 0)
    print(f"  Minor allele sites with ch2=0: {invalid_minor} (should be 0)")
    
    all_valid = (invalid_major == 0) and (invalid_missing == 0) and (invalid_minor == 0)
    print(f"  Overall: {'VALID ✓' if all_valid else 'INVALID ✗'}")
    
    return all_valid


def main():
    parser = argparse.ArgumentParser(description='Visualize haplotype images')
    parser.add_argument('input', type=str, help='Path to .npy file')
    parser.add_argument('--indices', type=int, nargs='+', default=None,
                        help='Specific indices to plot (default: first N images)')
    parser.add_argument('--n', type=int, default=5,
                        help='Number of images to plot if --indices not specified')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Output directory for images')
    parser.add_argument('--grid', action='store_true',
                        help='Create a single grid image instead of individual files')
    parser.add_argument('--validate', action='store_true',
                        help='Validate encoding consistency')
    parser.add_argument('--prefix', type=str, default='haplotype',
                        help='Filename prefix')
    
    args = parser.parse_args()
    
    # Load data
    print(f"Loading {args.input}...")
    data = np.load(args.input, mmap_mode='r')
    print(f"Shape: {data.shape}")
    print(f"  {data.shape[0]} images")
    print(f"  {data.shape[1]} samples per image")
    print(f"  {data.shape[2]} sites per image")
    print(f"  {data.shape[3]} channels")
    
    # Validate encoding if requested
    if args.validate:
        validate_encoding(data)
        print()
    
    # Determine which indices to plot

    if args.indices:
        indices = args.indices
    else:
        # REVISED: Generate n indices evenly spaced across the entire dataset
        # This prevents only seeing the first few (often similar) simulations
        num_to_plot = min(args.n, data.shape[0])
        indices = np.linspace(0, data.shape[0] - 1, num_to_plot, dtype=int).tolist()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    if args.grid:
        outpath = output_dir / f"{args.prefix}_grid.png"
        plot_grid(data, indices, outpath, title=Path(args.input).stem)
        print(f"Saved grid: {outpath}")
    else:
        for idx in indices:
            outpath = plot_haplotype_image(data, idx, output_dir, args.prefix)
            print(f"Saved: {outpath}")


if __name__ == '__main__':
    main()