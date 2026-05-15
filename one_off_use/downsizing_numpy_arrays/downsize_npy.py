import numpy as np
import sys
import argparse
from tqdm import tqdm

def symmetric_crop(data, site_map, target_width):
    """Crops the site dimension (axis 1) and map (if it exists) symmetrically."""
    curr_width = data.shape[1]
    
    if curr_width < target_width:
        return None, None 

    start_idx = (curr_width - target_width) // 2
    end_idx = start_idx + target_width
    
    cropped_data = data[:, start_idx:end_idx, :]
    
    cropped_map = None
    if site_map is not None:
        cropped_map = site_map[start_idx:end_idx]
    
    return cropped_data, cropped_map

def main():
    parser = argparse.ArgumentParser(description="Downsample, filter invariants, and crop haplotype images (Map optional).")
    parser.add_argument("input_npy", help="Path to the .npy image file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--input_map", help="Optional: Path to the corresponding _map.npy file", default=None)
    parser.add_argument("--new_samples", type=int, required=True, help="Number of samples to select")
    parser.add_argument("--target_width", type=int, required=True, help="Final number of sites")
    
    args = parser.parse_args()

    # --- Load Data ---
    print(f"Loading {args.input_npy}...")
    data = np.load(args.input_npy, allow_pickle=True)
    
    site_maps = None
    if args.input_map:
        print(f"Loading map: {args.input_map}")
        site_maps = np.load(args.input_map, allow_pickle=True)
    
    num_windows, total_samples, original_sites, channels = data.shape
    print(f"Original shape: {data.shape}")

    if args.new_samples > total_samples:
        print(f"Error: Requested {args.new_samples} samples, but only {total_samples} available.")
        sys.exit(1)

    processed_images = []
    processed_maps = []

    # --- Processing ---
    for i in tqdm(range(num_windows), desc="Processing Windows"):
        window = data[i] 
        smap = site_maps[i] if site_maps is not None else None

        # 1. Randomly select new samples
        sample_indices = np.random.choice(total_samples, args.new_samples, replace=False)
        window = window[sample_indices, :, :]

        # 2. Remove Invariant Sites
        genotypes = window[:, :, 0]
        is_variant = np.any(genotypes == 1, axis=0)
        
        window = window[:, is_variant, :]
        if smap is not None:
            smap = smap[is_variant]

        # 3. Symmetric Cropping
        cropped_window, cropped_map = symmetric_crop(window, smap, args.target_width)

        if cropped_window is not None:
            processed_images.append(cropped_window)
            if cropped_map is not None:
                processed_maps.append(cropped_map)

    # --- Save Results ---
    if not processed_images:
        print("Error: No windows survived.")
        sys.exit(1)

    final_images = np.array(processed_images, dtype=np.int8)
    img_out = f"{args.output_prefix}.npy"
    np.save(img_out, final_images)
    print(f"\nSaved images to: {img_out} (Shape: {final_images.shape})")

    if processed_maps:
        final_maps = np.array(processed_maps, dtype=np.int32)
        map_out = f"{args.output_prefix}_map.npy"
        np.save(map_out, final_maps)
        print(f"Saved maps to: {map_out}")

if __name__ == "__main__":
    main()