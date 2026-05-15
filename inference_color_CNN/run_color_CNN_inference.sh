#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=10:00:00,h_data=40G,highp
#$ -pe shared 1

# Takes a trained color CNN and a filtered+sorted haplotype file and outputs txt results
# =============================================================================

# 1. Load the job environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean


# 3. Specific File Paths
INPUT_MODEL_JSON="example_model.json"
INPUT_MODEL_WEIGHTS="example_model.weights.h5"
INPUT_SORTED_FILTERED_NPY="example_color.npy"

OUTPUT_FILE="example_predictions.txt"
PYTHON_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/inference_color_CNN/color_CNN_inference.py"

echo "=============================================="
echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started on: $(date)"
echo "=============================================="
echo ""

# 4. Verify the input directory exists
if [ -f "$INPUT_SORTED_FILTERED_NPY" ]; then
    echo "Starting prediction..."
    echo "----------------------------------------------"
    
    python "$PYTHON_SCRIPT" \
        "$INPUT_MODEL_JSON" \
        "$INPUT_MODEL_WEIGHTS" \
        "$INPUT_SORTED_FILTERED_NPY" \
        --output "$OUTPUT_FILE"
        
    echo "----------------------------------------------"
    echo "Processing complete."
else
    echo "ERROR: Directory '$INPUT_SORTED_FILTERED_NPY' not found."
    exit 1
fi

echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="