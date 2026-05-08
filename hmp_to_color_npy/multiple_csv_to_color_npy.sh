#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=200:00:00,h_data=50G,highp
#$ -pe shared 1

# Converts a specific HMP csv file to a numpy array for CNN training.
# =============================================================================

# Load the job environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# --- Configuration Parameters ---
# Set to 354 per your mentors' recent recommendation
WINDOW_H=354
SLIDE_STEP=1
TARGET_SAMPLES=120
SNP_THRESH=120
COMPLETE_THRESH=0.5

# --- Specific File Paths ---
# Replace these with the exact filenames you need to process
INPUT_FILE="example_input_data.csv"
OUTPUT_FILE="processed_haplotype_image.npy"
PYTHON_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/hmp_to_color_npy/HMP_csv_to_numpy.py"

echo "=============================================="
echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started on: $(date)"
echo "Input File:     $INPUT_FILE"
echo "SNP Window:     $WINDOW_H"
echo "Target Samples: $TARGET_SAMPLES"
echo "=============================================="
echo ""

# Verify the file exists before running the Python script
if [ -f "$INPUT_FILE" ]; then
    echo "Starting conversion..."
    echo "----------------------------------------------"
    
    python "$PYTHON_SCRIPT" \
        "$INPUT_FILE" \
        "$OUTPUT_FILE" \
        --window_h "$WINDOW_H" \
        --slide_step "$SLIDE_STEP" \
        --target_samples "$TARGET_SAMPLES" \
        --snp_threshold "$SNP_THRESH" \
        --complete_threshold "$COMPLETE_THRESH" \
        --sort rows_freq
        
    echo "----------------------------------------------"
    echo "Conversion complete."
else
    echo "ERROR: File '$INPUT_FILE' not found in the current directory."
    exit 1
fi

echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="