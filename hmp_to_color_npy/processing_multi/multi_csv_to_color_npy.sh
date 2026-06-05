#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=200:00:00,h_data=50G,highp
#$ -pe shared 1

# Converts a directory of HMP csv files to numpy arrays for CNN training.
# =============================================================================

# Load the job environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# --- Configuration Parameters ---
WINDOW_H=90
SLIDE_STEP=1
TARGET_SAMPLES=50
SNP_THRESH=50
COMPLETE_THRESH=0.7

# --- Specific File Paths ---
INPUT_DIRECTORY="/u/home/b/baeria/project-ngarud/Research/WideSimsCNN/A_finegoldii/sorted_csv"
INPUT_PREFIX="A_finegoldii_"

OUTPUT_DIRECTORY="output_numpy"
OUTPUT_PREFIX="A_fine_numpy_"

PYTHON_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/hmp_to_color_npy/HMP_csv_to_numpy.py"

echo "=============================================="
echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started on: $(date)"
echo "Input Directory: $INPUT_DIRECTORY"
echo "SNP Window:      $WINDOW_H"
echo "Target Samples:  $TARGET_SAMPLES"
echo "=============================================="
echo ""

# Verify the input directory exists
if [ -d "$INPUT_DIRECTORY" ]; then
    echo "Input directory found. Starting batch conversion..."
    echo "----------------------------------------------"
    
    # Create the output directory if it doesn't exist
    mkdir -p "$OUTPUT_DIRECTORY"
    
    # Loop through all matching .csv files in the input directory
    for INPUT_CSV in "$INPUT_DIRECTORY"/"$INPUT_PREFIX"*.csv; do
        
        # Guard clause in case the directory is empty or no files match
        if [ ! -e "$INPUT_CSV" ]; then
            echo "No files matching '${INPUT_PREFIX}*.csv' found."
            break
        fi

        # Extract the base filename without the path and .csv extension
        BASENAME=$(basename "$INPUT_CSV" .csv)
        
        # Strip the base prefix (e.g., "A_finegoldii_") to isolate the suffix
        SUFFIX=${BASENAME#$INPUT_PREFIX}
        
        # Construct the output file path using the new prefix and the extracted suffix
        CURRENT_OUTPUT_FILE="$OUTPUT_DIRECTORY/${OUTPUT_PREFIX}${SUFFIX}.npy"
        
        echo "Converting: $BASENAME.csv -> ${OUTPUT_PREFIX}${SUFFIX}.npy"
        
        # Run conversion
        python "$PYTHON_SCRIPT" \
            "$INPUT_CSV" \
            "$CURRENT_OUTPUT_FILE" \
            --window_h "$WINDOW_H" \
            --slide_step "$SLIDE_STEP" \
            --target_samples "$TARGET_SAMPLES" \
            --snp_threshold "$SNP_THRESH" \
            --complete_threshold "$COMPLETE_THRESH" \
            --sort rows_freq
            
    done
        
    echo "----------------------------------------------"
    echo "Batch conversion complete."
else
    echo "ERROR: Directory '$INPUT_DIRECTORY' not found."
    exit 1
fi

echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="