#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=10:00:00,h_data=40G,highp
#$ -pe shared 1

# Process 1000 SLiM replicate CSVs into a single centered middle-window numpy array.
# =============================================================================

# 1. Load the job environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# 2. Configuration Parameters
WINDOW_H=354
TARGET_SAMPLES=120
COMPLETE_THRESH=0.7
SORT_ORDER="rows_freq"

# 3. Specific File Paths
# Update INPUT_DIR to the folder containing your rep_1.1...rep_500.2 files
INPUT_DIR="ms_slimulations_results/soft"
OUTPUT_FILE="soft_sims.npy"
PYTHON_SCRIPT="util_scripts/many_haplotypes_to_numpy.py"

echo "=============================================="
echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started on: $(date)"
echo "Input Directory:  $INPUT_DIR"
echo "Window Size:      $WINDOW_H (Centered)"
echo "Target Samples:   $TARGET_SAMPLES"
echo "Sort Order:       $SORT_ORDER"
echo "=============================================="
echo ""

# 4. Verify the input directory exists
if [ -d "$INPUT_DIR" ]; then
    echo "Starting batch processing of replicates..."
    echo "----------------------------------------------"
    
    python "$PYTHON_SCRIPT" \
        "$INPUT_DIR" \
        "$OUTPUT_FILE" \
        --window_h "$WINDOW_H" \
        --target_samples "$TARGET_SAMPLES" \
        --complete_threshold "$COMPLETE_THRESH" \
        --sort "$SORT_ORDER"
        
    echo "----------------------------------------------"
    echo "Processing complete."
else
    echo "ERROR: Directory '$INPUT_DIR' not found."
    exit 1
fi

echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="
