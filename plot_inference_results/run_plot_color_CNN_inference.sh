#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=10:00:00,h_data=40G,highp
#$ -pe shared 1

# 1. Load the job environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# 2. Specific File Paths
# Make sure this points to where you saved the python code you showed me!
PYTHON_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/plot_inference_results/plot_inference_results.py" 

INPUT_PREDICTIONS="example_preds.txt"
INPUT_SITES_INDICES="example_freq_map.npy"
TITLE="R. bromii (HMP, n=120) Color CNN Genomic Scan"
BIN_SIZE=10
OUTPUT_FILE="example_scan.png"

echo "=============================================="
echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started on: $(date)"
echo "=============================================="

# 3. Verify the input exists
if [ -f "$INPUT_PREDICTIONS" ]; then
    echo "Starting Plotting..."
    echo "----------------------------------------------"
    
    # FIXED: Added the --site_indices flag and used the defined PYTHON_SCRIPT variable
    python "$PYTHON_SCRIPT" \
        "$INPUT_PREDICTIONS" \
        --site_indices "$INPUT_SITES_INDICES" \
        --title "$TITLE" \
        --bin_size $BIN_SIZE \
        --output "$OUTPUT_FILE"
        
    echo "----------------------------------------------"
    echo "Processing complete."
else
    echo "ERROR: File '$INPUT_PREDICTIONS' not found."
    exit 1
fi

echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="