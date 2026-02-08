# downsample_haplotypes.py

Reduces the number of samples (rows) per haplotype image by removing samples with the most missing data first. Ties are broken randomly.

## Why

When combining datasets with different sample counts (e.g., SLiM simulations with 154 samples and HMP data with 200 samples), you may need to downsample to a common size for CNN training. This script preferentially drops the lowest-quality samples (most missing data in channel 0).

## Location

```
/u/project/ngarud/Garud_lab/Brendan/Utils/downsample_haplotypes/downsample_haplotypes.py
```

### Basic (154 → 120 samples)

```bash
python "/u/project/ngarud/Garud_lab/Brendan/Utils/downsample_haplotypes/downsample_haplotypes.py" hard_sorted_color.npy 120 hard_120.npy
```

### With reproducible tie-breaking

```bash
python "/u/project/ngarud/Garud_lab/Brendan/Utils/downsample_haplotypes/downsample_haplotypes.py" hard_sorted_color.npy 120 hard_120.npy --seed 42
```

### All options

| Argument           | Default | Description                                          |
|--------------------|---------|------------------------------------------------------|
| `input_npy`        | —       | Input `.npy` file (shape: `N × samples × sites × 2`)|
| `target_samples`   | —       | Number of samples to keep per image                  |
| `output_npy`       | —       | Output `.npy` file path                              |
| `--seed`           | None    | Random seed for reproducible tie-breaking            |

## Example output

```
Loading hard_sorted_color.npy...
  Input shape: (30000, 154, 200, 2), dtype: int8
Downsampling 154 -> 120 samples per image...
  Output shape: (30000, 120, 200, 2), dtype: int8
Saved to hard_120.npy
```

## How it works

For each image, the script counts missing data per sample (channel 0 == 0), then keeps the `target_samples` with the fewest missing sites. Samples with equal missing counts are broken randomly. Original row order is preserved among kept samples.

## Note on sorting

If your images were sorted (by frequency or distance) before downsampling, the sort order will be partially disrupted since rows are removed. You may want to re-sort after downsampling.