#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=24:00:00,h_data=30G,highp,gpu
#$ -pe shared 1

# --- Load Environment ---
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# Echo job info
echo "Job $JOB_ID started on: `hostname -s` at `date`"
echo "Arguments provided: Input=$1, Output=$2, Samps=$3, Window=$4"

# --- Dynamic Configuration ---
# Internal paths (stable)
HELPER_PATH="/u/project/ngarud/Garud_lab/Brendan/Utils/helper_processSLiMsims.py"
MANAGER_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/SLiMsims_to_numpy.py"

# User-provided arguments (dynamic)
# $1 = INPUT_FOLDER (e.g., /u/home/b/baeria/.../hard/)
# $2 = OUTPUT_NAME  (e.g., hard_sorted_color.npy)
# $3 = NUM_SAMPS    (e.g., 154)
# $4 = WINDOW_SIZE  (e.g., 200)

INPUT_FOLDER="$1"
OUTPUT_NAME="$2"
NUM_SAMPS="${3:-154}" # Defaults to 154 if $3 is not provided
WINDOW_SIZE="${4:-200}" # Defaults to 200 if $4 is not provided


# --- Run Code ---
# Note: I'm keeping the 200 window size and 2 channels as defaults based on your example
python "$MANAGER_SCRIPT" \
    "$OUTPUT_NAME" \
    "$INPUT_FOLDER" \
    --num_samps "$NUM_SAMPS" \
    --window_size "$WINDOW_SIZE" \
    --channels 2 \
    --sort \
    --helper_path "$HELPER_PATH"

echo " "
echo "Job $JOB_ID ended on: `hostname -s` at `date`"