#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=200:00:00,h_data=50G,highp
#$ -pe shared 1

# Converts a directory of raw HMP csv files into formatted/recoded csv files.
# =============================================================================

# Load the job environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# --- Specific File Paths ---
INPUT_DIRECTORY="/u/home/b/baeria/project-ngarud/Research/WideSimsCNN/A_finegoldii/sorted_csv"
INPUT_PREFIX="A_finegoldii_"

OUTPUT_DIRECTORY="output_formatted_csv"
OUTPUT_PREFIX="A_fine_formatted_"

PYTHON_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/one_off_use/process_contig_csv_files/recoding_and_column_formatting.py"

echo "=============================================="
echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started on: $(date)"
echo "Input Directory:  $INPUT_DIRECTORY"
echo "Output Directory: $OUTPUT_DIRECTORY"
echo "=============================================="
echo ""

# Verify the input directory exists
if [ -d "$INPUT_DIRECTORY" ]; then
    echo "Input directory found. Starting batch formatting..."
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
        CURRENT_OUTPUT_FILE="$OUTPUT_DIRECTORY/${OUTPUT_PREFIX}${SUFFIX}.csv"
        
        echo "Formatting: $BASENAME.csv -> ${OUTPUT_PREFIX}${SUFFIX}.csv"
        
        # Run conversion script (only requires input and output arguments)
        python "$PYTHON_SCRIPT" \
            "$INPUT_CSV" \
            "$CURRENT_OUTPUT_FILE"
            
    done
        
    echo "----------------------------------------------"
    echo "Batch formatting complete."
else
    echo "ERROR: Directory '$INPUT_DIRECTORY' not found."
    exit 1
fi

echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="