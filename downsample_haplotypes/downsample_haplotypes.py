#!/usr/bin/env python3
"""
Downsample haplotype images by removing samples with the most missing data.

Takes a numpy array of shape (N, samples, sites, 2) and reduces the samples
dimension by preferentially removing samples with the highest missing data count
(channel 0 == 0). Ties are broken randomly.

Usage:
    python downsample_haplotypes.py input.npy 120 output.npy
"""

import numpy as np
import argparse
import sys


def downsample(data, target_samples, seed=None):
    """
    Downsample along axis 1 by removing samples with most missing data.
    
    Args:
        data: np.ndarray of shape (N, samples, sites, 2)
        target_samples: int, number of samples to keep
        seed: optional random seed for reproducibility
        
    Returns:
        np.ndarray of shape (N, target_samples, sites, 2)
    """
    n_images, n_samples, n_sites, n_channels = data.shape
    
    if target_samples >= n_samples:
        print(f"Warning: target ({target_samples}) >= current samples ({n_samples}). No downsampling needed.")
        return data
    
    rng = np.random.default_rng(seed)
    out = np.empty((n_images, target_samples, n_sites, n_channels), dtype=data.dtype)
    
    for i in range(n_images):
        # Count missing data per sample (channel 0 == 0 means missing)
        missing_counts = np.sum(data[i, :, :, 0] == 0, axis=1)  # (n_samples,)
        
        # Break ties randomly by adding small random jitter to sort key
        jitter = rng.random(n_samples)
        sort_key = missing_counts + jitter  # lower = less missing = keep
        
        # Keep the target_samples with lowest missing count
        keep_indices = np.argpartition(sort_key, target_samples)[:target_samples]
        keep_indices.sort()  # preserve original row order
        
        out[i] = data[i, keep_indices]
    
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Downsample haplotype images by removing samples with most missing data."
    )
    parser.add_argument("input_npy", help="Input .npy file (shape: N x samples x sites x 2)")
    parser.add_argument("target_samples", type=int, help="Number of samples to keep per image")
    parser.add_argument("output_npy", help="Output .npy file path")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for tie-breaking")
    
    args = parser.parse_args()
    
    print(f"Loading {args.input_npy}...")
    data = np.load(args.input_npy)
    print(f"  Input shape: {data.shape}, dtype: {data.dtype}")
    
    if data.ndim != 4:
        print(f"Error: Expected 4D array (N, samples, sites, channels), got {data.ndim}D")
        sys.exit(1)
    
    if args.target_samples >= data.shape[1]:
        print(f"Error: target_samples ({args.target_samples}) must be less than current samples ({data.shape[1]})")
        sys.exit(1)
    
    print(f"Downsampling {data.shape[1]} -> {args.target_samples} samples per image...")
    result = downsample(data, args.target_samples, seed=args.seed)
    
    np.save(args.output_npy, result)
    print(f"  Output shape: {result.shape}, dtype: {result.dtype}")
    print(f"Saved to {args.output_npy}")


if __name__ == "__main__":
    main()