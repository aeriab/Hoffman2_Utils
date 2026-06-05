#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=20:00:00,h_data=64G,highp
#$ -pe shared 1

# Takes a trained color CNN and a directory of filtered+sorted haplotype files and outputs txt results
# =============================================================================

# 1. Load the job environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# 3. Specific File Paths
INPUT_MODEL_JSON="example_model.json"
INPUT_MODEL_WEIGHTS="example_model.weights.h5"
INPUT_SORTED_FILTERED_NPY_DIRECTORY="example_color_npy_directory"
NUMPY_FILE_NAME="example_color"

OUTPUT_DIRECTORY="example_prediction_directory"
PYTHON_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/inference_color_CNN/color_CNN_inference.py"

echo "=============================================="
echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started on: $(date)"
echo "=============================================="
echo ""

# 4. Verify the input directory exists (-d checks for directory instead of file)
if [ -d "$INPUT_SORTED_FILTERED_NPY_DIRECTORY" ]; then
    echo "Input directory found. Starting batch predictions..."
    echo "----------------------------------------------"
    
    # Create the output directory if it doesn't already exist
    mkdir -p "$OUTPUT_DIRECTORY"
    
    # Loop through all matching .npy files in the input directory
    for INPUT_FILE in "$INPUT_SORTED_FILTERED_NPY_DIRECTORY"/"$NUMPY_FILE_NAME"*.npy; do
        
        # Guard clause in case the directory is empty or no files match
        if [ ! -e "$INPUT_FILE" ]; then
            echo "No files matching '${NUMPY_FILE_NAME}*.npy' found."
            break
        fi

        # Extract the base filename without the path and without the .npy extension (e.g., "example_color01")
        BASENAME=$(basename "$INPUT_FILE" .npy)
        
        # Strip the base "example_color" part to isolate just the suffix (e.g., "01")
        SUFFIX=${BASENAME#$NUMPY_FILE_NAME}
        
        # Construct the output file path using the extracted suffix
        CURRENT_OUTPUT_FILE="$OUTPUT_DIRECTORY/example_predictions${SUFFIX}.txt"
        
        echo "Processing: $BASENAME.npy -> example_predictions${SUFFIX}.txt"
        
        # Run inference
        python "$PYTHON_SCRIPT" \
            "$INPUT_MODEL_JSON" \
            "$INPUT_MODEL_WEIGHTS" \
            "$INPUT_FILE" \
            --output "$CURRENT_OUTPUT_FILE"
            
    done
        
    echo "----------------------------------------------"
    echo "Processing complete."
else
    echo "ERROR: Directory '$INPUT_SORTED_FILTERED_NPY_DIRECTORY' not found."
    exit 1
fi

echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="