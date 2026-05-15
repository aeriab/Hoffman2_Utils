#!/bin/bash
#$ -cwd
#$ -o joblog_sort.$JOB_ID
#$ -j y
#$ -l h_rt=08:00:00,h_data=16G
#$ -m abe

# 1. Load the environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate tf_A100_clean

# 2. Define Paths
SORT_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/sim_to_color_numpy_care_package/util_scripts/remajor_and_middle_sort.py"
INPUT_DATA="your_input_here.npy"  # Replace with your actual input filename
OUTPUT_DATA="example.npy"

# 3. Process the dataset
if [ -f "$INPUT_DATA" ]; then
    echo "Processing: $INPUT_DATA"
    python "$SORT_SCRIPT" "$INPUT_DATA" "$OUTPUT_DATA"
    echo "Processing complete. Saved to $OUTPUT_DATA"
else
    echo "Error: $INPUT_DATA not found."
    exit 1
fi