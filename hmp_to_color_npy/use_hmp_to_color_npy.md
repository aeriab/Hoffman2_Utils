## How to convert HMP csv files to numpy arrays (Color Encoding v2):

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

> **Note:** Channel 1 only encodes syn/nonsyn at sites where the sample carries the minor allele. All other sites are 0.

---

### Standard (sorted, downsampled to 120):
```bash
python /u/project/ngarud/Garud_lab/Brendan/Utils/hmp_to_color_npy/HMP_csv_to_numpy.py data.csv output.npy --window_h 200 --slide_step 10 --target_samples 120 --sort
```

### Minimalist (no sorting, full samples):
```bash
python /u/project/ngarud/Garud_lab/Brendan/Utils/hmp_to_color_npy/HMP_csv_to_numpy.py data.csv output.npy
```

### All Options:
| Argument            | Default | Description                                              |
|---------------------|---------|----------------------------------------------------------|
| `input_csv`         | —       | Path to input CSV                                        |
| `output_npy`        | —       | Path for output `.npy` file                              |
| `--window_h`        | 200     | Window size (number of sites)                            |
| `--slide_step`      | 10      | Step size for sliding window                             |
| `--target_samples`  | None    | Downsample to N haplotypes (picks least missing data)    |
| `--sort`            | False   | Sort haplotypes by distance (requires `helper_haplotypeSorter.py`) |
| `--raw_missing_val` | -1      | Value representing missing data in the raw CSV           |

### Output
- `output.npy` — shape `(num_windows, num_samples, window_h, 2)`, dtype `int8`
- `output_map.npy` — shape `(num_windows, window_h)`, site positions for each window