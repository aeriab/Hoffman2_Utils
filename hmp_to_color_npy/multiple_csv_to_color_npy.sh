#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=200:00:00,h_data=50G,highp
#$ -pe shared 1

# Converts multiple HMP csv files to multiple numpy arrays for CNN training.
#
# Output Encoding:
#   Channel 1 (Allele State):  -1=major, 0=missing, 1=minor
#   Channel 2 (Mutation Type): -1=synonymous, 0=major/missing, 1=non-synonymous
# =============================================================================

# Example use:
# qsub script.sh /path/to/csvs A_finegoldii 55

# Load the job environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# Echo job info on joblog
echo "=============================================="
echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started on: $(date)"
echo "=============================================="
echo ""

# Input/Output
INPUT_FOLDER=$1
PREFIX=$2
MAX_NUM=$3

# Ensure the user provided arguments
if [ -z "$INPUT_FOLDER" ] || [ -z "$PREFIX" ] || [ -z "$MAX_NUM" ]; then
    echo "Usage: qsub script.sh <input_folder> <prefix_term> <max_number>"
    exit 1
fi

# Iterate from 01 to the max number
for i in $(seq -f "%02g" 1 "$MAX_NUM"); do
    FILE_NAME="${PREFIX}_${i}"
    INPUT_PATH="${INPUT_FOLDER}/${FILE_NAME}.csv"
    OUTPUT_PATH="${INPUT_FOLDER}/${FILE_NAME}.npy"

    if [ -f "$INPUT_PATH" ]; then
        echo "Processing: $FILE_NAME"
        python /u/project/ngarud/Garud_lab/Brendan/Utils/hmp_to_color_npy/HMP_csv_to_numpy.py \
            "$INPUT_PATH" "$OUTPUT_PATH" --window_h 200 --slide_step 1 --target_samples 120
    else
        echo "File $INPUT_PATH not found, skipping..."
    fi
done


# Echo job info on joblog
echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(hostname -s)"
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="