#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=200:00:00,h_data=50G,highp
#$ -pe shared 1

# Converts multiple HMP csv files to multiple numpy arrays for CNN training.
# Updated to support SNP thresholding and missing data filtering.
# =============================================================================

# Load the job environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# --- Configuration Parameters ---
WINDOW_H=200
SLIDE_STEP=1
TARGET_SAMPLES=120
SNP_THRESH=120        # Minimum SNPs in a window to keep it
COMPLETE_THRESH=0.5  # % of non-missing data required per sample (e.g., 0.5 = 50%)

# Input/Output Arguments
INPUT_PATH="example.csv"
OUTPUT_PATH="example.npy"

# Ensure the user provided arguments
if [ -z "$INPUT_FOLDER" ] || [ -z "$PREFIX" ]; then
    echo "Usage: qsub script.sh <input_folder> <prefix_term>"
    exit 1
fi

echo "=============================================="
echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started on: $(date)"
echo "Target Samples: $TARGET_SAMPLES"
echo "SNP Threshold: $SNP_THRESH"
echo "Complete Data Threshold: $COMPLETE_THRESH"
echo "=============================================="
echo ""

# Iterate from 01 to the max number
for i in $(seq -f "%02g" 1 "$MAX_NUM"); do

    if [ -f "$INPUT_PATH" ]; then
        echo "----------------------------------------------"
        
        # Calling the updated Python script with new threshold flags
        python /u/project/ngarud/Garud_lab/Brendan/Utils/hmp_to_color_npy/HMP_csv_to_numpy.py \
            "$INPUT_PATH" \
            "$OUTPUT_PATH" \
            --window_h "$WINDOW_H" \
            --slide_step "$SLIDE_STEP" \
            --target_samples "$TARGET_SAMPLES" \
            --snp_threshold "$SNP_THRESH" \
            --complete_threshold "$COMPLETE_THRESH" \
            --sort rows_freq
            
    else
        echo "File $INPUT_PATH not found, skipping..."
    fi
done

echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(hostname -s)"
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="