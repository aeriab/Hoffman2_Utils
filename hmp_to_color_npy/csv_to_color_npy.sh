#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=200:00:00,h_data=50G,highp
#$ -pe shared 1

# =============================================================================
# SLiM Simulations to Numpy Converter
# =============================================================================
# Converts SLiM .txt output files to numpy arrays for CNN training.
#
# Output Encoding:
#   Channel 1 (Allele State):  -1=major, 0=missing, 1=minor
#   Channel 2 (Mutation Type): -1=synonymous, 0=major/missing, 1=non-synonymous
# =============================================================================

# Load the job environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# Echo job info on joblog
echo "=============================================="
echo "Job $JOB_ID started on: $(hostname -s)"
echo "Job $JOB_ID started on: $(date)"
echo "=============================================="
echo ""

# Input/Output
INPUT_FOLDER=$1
OUTPUT_NAME=$2

# Ensure the user provided arguments
if [ -z "$INPUT_FOLDER" ] || [ -z "$OUTPUT_NAME" ]; then
    echo "Usage: qsub script.sh <input_csv> <output_npy_name>"
    exit 1
fi

python /u/project/ngarud/Garud_lab/Brendan/Utils/hmp_to_color_npy/HMP_csv_to_numpy.py "$INPUT_FOLDER" "$OUTPUT_NAME" --window_h 200 --slide_step 1 --target_samples 120

# Echo job info on joblog
echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(hostname -s)"
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="