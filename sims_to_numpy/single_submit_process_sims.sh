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

# --- Configuration ---
HELPER_PATH="/u/project/ngarud/Garud_lab/Brendan/Utils/sims_to_numpy/helper_processSLiMsims.py"
MANAGER_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/sims_to_numpy/SLiMsims_to_numpy.py"
VISUALIZER="/u/project/ngarud/Garud_lab/Brendan/Utils/numpy_image_generator/visualize_haplotypes.py"

# Input/Output
INPUT_FOLDER="/u/home/b/baeria/project-ngarud/hmp_SLiMulations/dann_slimulations_12080244/hard/"
OUTPUT_NAME="hard_sorted_color.npy"

# Parameters
NUM_SAMPS=154
WINDOW_SIZE=200
CHANNELS=2      # 1=allele state only, 2=allele state + mutation type

# --- Run Conversion ---
echo "Converting SLiM simulations to numpy..."
echo "  Input:  $INPUT_FOLDER"
echo "  Output: $OUTPUT_NAME"
echo "  Samples: $NUM_SAMPS, Window: $WINDOW_SIZE, Channels: $CHANNELS"
echo ""

python "$MANAGER_SCRIPT" \
    "$OUTPUT_NAME" \
    "$INPUT_FOLDER" \
    --num_samps $NUM_SAMPS \
    --window_size $WINDOW_SIZE \
    --channels $CHANNELS \
    --sort \
    --helper_path "$HELPER_PATH"

# --- Optional: Validate and Visualize ---
if [ -f "$OUTPUT_NAME" ]; then
    echo ""
    echo "Generating validation plots..."
    mkdir -p plots
    python "$VISUALIZER" \
        "$OUTPUT_NAME" \
        --validate \
        --n 5 \
        --grid \
        --prefix "hard_sweep" \
        --output-dir ./plots
    echo "Plots saved to ./plots/"
fi

# Echo job info on joblog
echo ""
echo "=============================================="
echo "Job $JOB_ID ended on: $(hostname -s)"
echo "Job $JOB_ID ended on: $(date)"
echo "=============================================="