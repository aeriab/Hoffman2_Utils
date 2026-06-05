#!/bin/bash

#$ -cwd
#$ -o plot_joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=20:00:00,h_data=40G
#$ -pe shared 1

# =============================================================================
# Takes a directory of sequential prediction text files and plots them
# =============================================================================

# 1. Load the job environment (adjust if you need your tf conda env, though standard anaconda3 usually has matplotlib/pandas)
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# 2. Specific File Paths
INPUT_PREDICTIONS_DIR="multi_prediction_directory"
PYTHON_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/plot_inference_results/multi_plot_inference_results.py"
OUTPUT_PLOT_NAME="genome_wide_scan_A_finegoldii.png"

echo "=============================================="
echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started on: $(date)"
echo "=============================================="
echo ""

# 3. Verify directory exists and run
if [ -d "$INPUT_PREDICTIONS_DIR" ]; then
    echo "Directory found. Starting plotting job..."
    
    python "$PYTHON_SCRIPT" \
        "$INPUT_PREDICTIONS_DIR" \
        --title "A. finegoldii Genome Wide Scan" \
        --output "$OUTPUT_PLOT_NAME" \
        --bin_size 3 \
        --y_max 10.0
        
    echo "Plotting complete."
else
    echo "ERROR: Directory '$INPUT_PREDICTIONS_DIR' not found."
    exit 1
fi

echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="