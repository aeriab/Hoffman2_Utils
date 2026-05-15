#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=08:00:00,h_data=30G
#$ -pe shared 1

# --- Load Environment ---
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

echo "Job $JOB_ID started on: `hostname -s` at `date`"

# --- Hard-Coded Configurations ---
# INPUT_NPY should be the path to your simulated data
INPUT_NPY="/u/project/ngarud/Garud_lab/Brendan/path/to/your/simulated_data.npy"
OUTPUT_PREFIX="downsized_simulated_output"

# Reprocessing Parameters
NEW_SAMPLES=50
TARGET_WIDTH=90

# Path to your downsizing utility script
DOWNSIZE_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/one_off_use/downsizing_numpy_arrays/downsize_npy.py"

# --- Run Downsizing ---
echo "Starting downsizing process for simulated data..."
echo "Input file: $INPUT_NPY"

# Note: We are NOT passing --input_map here
python "$DOWNSIZE_SCRIPT" \
    "$INPUT_NPY" \
    "$OUTPUT_PREFIX" \
    --new_samples "$NEW_SAMPLES" \
    --target_width "$TARGET_WIDTH"

echo "Job $JOB_ID ended on: `hostname -s` at `date`"