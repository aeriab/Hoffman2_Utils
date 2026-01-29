# Brendan's Genomics & CNN Utils

Shared utilities for SLiMulations and CNN sweep detection.

Anywhere on hoffman2, these data processing scripts can be used like so:

--

## How to convert SLiM simulation files to numpy arrays:
python /u/project/ngarud/Garud_lab/Brendan/Utils/SLiMsims_to_numpy.py <out.npy> <n_samps> <window> <in_dir> <channels>

Or, you can modify the contents of submit_process_sims.sh and run a job on highp:
cp /u/project/ngarud/Garud_lab/Brendan/Utils/submit_process_sims.sh submit_process_sims.sh
qsub submit_process_sims.sh

Ensure that the conda environment is set up properly. This line is specific to Brendan:
conda activate tf_A100_clean

--

## How to convert HMP csv files to numpy arrays:
### Standard (2 channels, sorted, downsampled to 120):
python /u/project/ngarud/Garud_lab/Brendan/Utils/HMP_csv_to_numpy.py data.csv output.npy --window_h 201 --slide_step 10 --channels 2 --target_samples 120 --sort

### Minimalist (1 channel, no sorting, full samples):
python /u/project/ngarud/Garud_lab/Brendan/Utils/HMP_csv_to_numpy.py data.csv output.npy --channels 1


--


## Get a subsample of a numpy of shape (batch#, sample#, site#, 2)

### This example gets a snapshot of the bw channel:
/u/project/ngarud/Garud_lab/Brendan/Utils/hap_sample.sh sims.npy my_test

### This example gets a snapshot of the color channel:
/u/project/ngarud/Garud_lab/Brendan/Utils/hap_sample.sh sims.npy my_test 1

