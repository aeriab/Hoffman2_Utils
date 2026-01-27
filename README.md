# Brendan's Genomics & CNN Utils

Shared utilities for SLiMulations and CNN sweep detection.

Anywhere on hoffman2, these data processing scripts can be used like so:

--

## How to convert SLiM simulation files to numpy arrays:
python /u/project/ngarud/Garud_lab/Brendan/Utils/SLiMsims_to_numpy.py <out.npy> <n_samps> <window> <in_dir> <channels>

--

## How to convert HMP csv files to numpy arrays:
### Standard (2 channels, sorted, downsampled to 100):
python HMP_csv_to_numpy.py data.csv output.npy --window_h 201 --slide_step 10 --channels 2 --target_samples 100 --sort

### Minimalist (1 channel, no sorting, full samples):
python HMP_csv_to_numpy.py data.csv output.npy --channels 1
