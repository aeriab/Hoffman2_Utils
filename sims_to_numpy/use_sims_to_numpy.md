# SLiM Simulations to Numpy Converter

Convert SLiM `.txt` output files into numpy arrays for CNN training.

## Location
```
/u/project/ngarud/Garud_lab/Brendan/Utils/
├── sims_to_numpy/
│   ├── SLiMsims_to_numpy.py          # Manager script
│   ├── helper_processSLiMsims.py     # Helper script (called by manager)
│   ├── submit_process_sims.sh        # Job submission script
│   └── use_sims_to_numpy.md          # This documentation
└── numpy_image_generator/
    └── visualize_haplotypes.py       # Visualization script
```

## Output Encoding

The converter produces numpy arrays with shape `(n_images, n_samples, n_sites, 2)`:

**Channel 1 (Allele State):**
| Value | Meaning |
|-------|---------|
| -1 | Major allele |
| 0 | Missing data |
| 1 | Minor allele |

**Channel 2 (Mutation Type):**
| Value | Meaning |
|-------|---------|
| -1 | Synonymous mutation (selection coefficient = 0) |
| 0 | Major allele / Missing data |
| 1 | Non-synonymous mutation (selection coefficient ≠ 0) |

## Usage

### Command Line
```bash
python /u/project/ngarud/Garud_lab/Brendan/Utils/sims_to_numpy/SLiMsims_to_numpy.py \
    <output.npy> \
    <input_directory> \
    --num_samps <N> \
    --window_size <W> \
    --channels <1|2> \
    --sort \
    --helper_path /u/project/ngarud/Garud_lab/Brendan/Utils/sims_to_numpy/helper_processSLiMsims.py
```

### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `output.npy` | Path for output numpy file | (required) |
| `input_directory` | Directory containing SLiM `.txt` files | (required) |
| `--num_samps` | Number of haplotypes to sample per simulation | 100 |
| `--window_size` | Window size in SNPs | 201 |
| `--channels` | 1 = allele state only, 2 = allele state + mutation type | 2 |
| `--sort` | Sort haplotypes by frequency | off |
| `--helper_path` | Path to helper script | (see above) |

### Examples
```bash
# Basic usage with defaults
python /u/project/ngarud/Garud_lab/Brendan/Utils/sims_to_numpy/SLiMsims_to_numpy.py \
    neutral_sims.npy \
    /path/to/neutral_slim_outputs/

# Full options
python /u/project/ngarud/Garud_lab/Brendan/Utils/sims_to_numpy/SLiMsims_to_numpy.py \
    hard_sweep_sims.npy \
    /path/to/hard_sweep_outputs/ \
    --num_samps 154 \
    --window_size 200 \
    --channels 2 \
    --sort \
    --helper_path /u/project/ngarud/Garud_lab/Brendan/Utils/sims_to_numpy/helper_processSLiMsims.py
```

## Submitting as a Job on Hoffman2

1. Copy the submission script to your working directory:
```bash
cp /u/project/ngarud/Garud_lab/Brendan/Utils/sims_to_numpy/submit_process_sims.sh ./submit_process_sims.sh
```

2. Edit the configuration section in `submit_process_sims.sh`:
```bash
# --- Configuration ---
INPUT_FOLDER="/path/to/your/slim/outputs/"
OUTPUT_NAME="your_output.npy"
NUM_SAMPS=154
WINDOW_SIZE=200
CHANNELS=2
```

3. Submit the job:
```bash
qsub submit_process_sims.sh
```

## Validation and Visualization

After conversion, validate the output and generate sample plots:
```bash
python /u/project/ngarud/Garud_lab/Brendan/Utils/numpy_image_generator/visualize_haplotypes.py \
    your_output.npy \
    --validate \
    --n 5 \
    --grid \
    --output-dir ./plots
```

See [visualize_haplotypes documentation](../numpy_image_generator/use_visualize_haplotypes.md) for more options.

## Environment Setup

Ensure the conda environment is activated before running:
```bash
module load anaconda3
conda activate base
conda activate tf_A100_clean
```

**Note:** The `tf_A100_clean` environment is specific to Brendan's setup. You may need to create your own environment with `numpy`, `matplotlib`, and `tqdm` installed.

## Troubleshooting

**"No .txt files found"**
- Check that `input_directory` contains SLiM output files with `.txt` extension

**"Window size > available SNPs"**
- The simulation has fewer SNPs than requested window size
- Reduce `--window_size` or regenerate simulations with more mutations

**"Requested N samples but only M available"**
- The simulation has fewer haplotypes than requested
- Reduce `--num_samps` to match available haplotypes

**Encoding validation fails**
- Run the visualizer with `--validate` to diagnose
- Check that SLiM output format matches expected structure