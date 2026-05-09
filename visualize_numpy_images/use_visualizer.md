# Haplotype Image Visualizer

Visualize numpy haplotype arrays as PNG images directly on Hoffman2.

## Location
```
/u/project/ngarud/Garud_lab/Brendan/Utils/visualize_numpy_images/visualize_haplotypes.py
```

## Data Encoding

The visualizer expects numpy arrays with shape `(n_images, n_samples, n_sites, 2)`:

### Encoding Scheme
**Channel 0 ("Black & White" - Genotype):**
| Value | Meaning       | Color    | Hex       |
|-------|---------------|----------|-----------|
| -1    | Major allele  | CHARCOAL | `#2F3640` |
|  0    | Missing data  | GREY     | `#BDC3C7` |
|  1    | Minor allele  | BLUE     | `#3498DB` |

**Channel 1 ("Color" - Mutation Type):**
| Value | Meaning                | Color    | Hex       |
|-------|------------------------|----------|-----------|
| -1    | Synonymous mutation    | YELLOW   | `#E1B12C` |
|  0    | Major allele / missing | CHARCOAL | `#2F3640` |
|  1    | Non-synonymous mutation| GREEN    | `#44BD32` |

## Usage

### Basic: Plot first 10 images individually
```bash
python /u/project/ngarud/Garud_lab/Brendan/Utils/visualize_numpy_images/visualize_haplotypes.py Neu_sims.npy --output-dir ./plots
```

### Generates 10 images sampled at equal intervals across the entire dataset.
```bash
python /u/project/ngarud/Garud_lab/Brendan/Utils/visualize_numpy_images/visualize_haplotypes.py input.npy --n 10 --grid --prefix "spaced_samples"
```

### Grid view of first 10 images
```bash
python /u/project/ngarud/Garud_lab/Brendan/Utils/visualize_numpy_images/visualize_haplotypes.py \
    Neu_sims.npy --grid --output-dir ./plots
```

### Plot specific image indices
```bash
python /u/project/ngarud/Garud_lab/Brendan/Utils/visualize_numpy_images/visualize_haplotypes.py hard_sorted_color.npy --n 50 --grid --prefix "hard_sweep_sample"
```

### Validate encoding + visualize
```bash
python /u/project/ngarud/Garud_lab/Brendan/Utils/visualize_numpy_images/visualize_haplotypes.py \
    Neu_sims.npy --validate --n 5 --grid
```

### Compare all three selection classes
```bash
# Neutral
python /u/project/ngarud/Garud_lab/Brendan/Utils/visualize_numpy_images/visualize_haplotypes.py \
    Neu_sims.npy --indices 0 1 2 --grid --prefix neutral --output-dir ./plots

# Hard sweep
python /u/project/ngarud/Garud_lab/Brendan/Utils/visualize_numpy_images/visualize_haplotypes.py \
    HS_sims.npy --indices 0 1 2 --grid --prefix hard_sweep --output-dir ./plots

# Soft sweep
python /u/project/ngarud/Garud_lab/Brendan/Utils/visualize_numpy_images/visualize_haplotypes.py \
    SS_sims.npy --indices 0 1 2 --grid --prefix soft_sweep --output-dir ./plots
```

## Options

| Flag | Description | Default |
|------|-------------|---------|
| `--output-dir` | Directory to save output PNGs | `.` (current dir) |
| `--indices` | Specific image indices to plot | First N images |
| `--n` | Number of random images (if `--indices` not set) | 10 |
| `--grid` | Combine images into single grid PNG | Individual files |
| `--prefix` | Filename prefix | `haplotype` |
| `--validate` | Check encoding consistency | Off |

## Output

- **Individual mode:** `{prefix}_{index:04d}.png` for each image
- **Grid mode:** `{prefix}_grid.png` with all images in one file

## Downloading Results

After generating plots on Hoffman2:
```bash
scp username@hoffman2.idre.ucla.edu:/path/to/plots/*.png ./local_plots/
```

## Requirements

- numpy
- matplotlib

Both should be available in your `tf_A100_clean` conda environment.