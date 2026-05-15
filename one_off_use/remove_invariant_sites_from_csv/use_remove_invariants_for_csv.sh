#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=08:00:00,h_data=16G
#$ -pe shared 1

# --- Load Environment ---
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
# Note: Keeping tf_A100_clean if it contains your pandas installation, 
# though 'base' usually suffices for this task.
conda activate tf_A100_clean

echo "Job $JOB_ID started on: `hostname -s` at `date`"

# --- CONFIGURATION ---
# Filenames (assumed to be in the current directory)
INPUT_NAME="cropped_r_bromii_data.csv"
OUTPUT_NAME="filtered_cropped_r_bromii.csv"

# Path to the utility script
PY_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/one_off_use/remove_invariant_sites_from_csv/remove_invariants_for_csv.py"
# ---------------------

echo "Starting invariant site removal..."
echo "Input: ./$INPUT_NAME"

# Execute the Python script
python3 "$PY_SCRIPT" "$INPUT_NAME" "$OUTPUT_NAME"

echo "Done. Cleaned file saved as: $OUTPUT_NAME"
echo "Job $JOB_ID ended on: `hostname -s` at `date`"