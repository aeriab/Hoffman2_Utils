#!/bin/bash
#$ -cwd
#$ -o joblog_sort.$JOB_ID
#$ -j y
#$ -l h_rt=10:00:00,h_data=40G
#$ -m abe

# 1. Load the environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate tf_A100_clean

# 2. Define Paths
SORT_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/sim_to_color_numpy_care_package/util_scripts/remajor_and_middle_sort.py"
INPUT_DIRECTORY="output_numpy"
INPUT_DATA="A_fine_numpy_"  # Replace with your actual input filename

OUTPUT_DIRECTORY="flipped_numpy"
OUTPUT_DATA="flipped_A_fine_"

echo "=============================================="
echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started on: $(date)"
echo "=============================================="
echo ""

# 4. Verify the input directory exists (-d checks for directory instead of file)
if [ -d "$INPUT_DIRECTORY" ]; then
    echo "Input directory found. Starting batch predictions..."
    echo "----------------------------------------------"
    
    # Create the output directory if it doesn't already exist
    mkdir -p "$OUTPUT_DIRECTORY"
    
    # Loop through all matching .npy files in the input directory
    for INPUT_FILE in "$INPUT_DIRECTORY"/"$INPUT_DATA"*.npy; do
        
        # Guard clause in case the directory is empty or no files match
        if [ ! -e "$INPUT_FILE" ]; then
            echo "No files matching '${INPUT_DATA}*.npy' found."
            break
        fi

        # Extract the base filename without the path and without the .npy extension (e.g., "example_color01")
        BASENAME=$(basename "$INPUT_FILE" .npy)
        
        # Strip the base "example_color" part to isolate just the suffix (e.g., "01")
        SUFFIX=${BASENAME#$INPUT_DATA}
        
        # Construct the output file path using the extracted suffix
        CURRENT_OUTPUT_FILE="$OUTPUT_DIRECTORY/${OUTPUT_DATA}${SUFFIX}.npy"
        
        echo "Processing: $BASENAME.npy -> example_predictions${SUFFIX}.txt"
        
        # Run inference
        python "$SORT_SCRIPT" "$INPUT_FILE" "$CURRENT_OUTPUT_FILE"
            
    done
        
    echo "----------------------------------------------"
    echo "Processing complete."
else
    echo "ERROR: Directory '$INPUT_DIRECTORY' not found."
    exit 1
fi

echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="