#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=200:00:00,h_data=50G,gpu,highp
#$ -pe shared 1

# --- Load Environment ---
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

echo "Job $JOB_ID started on: $(hostname -s) at $(date)"

# --- Configuration: Paths ---
HELPER_PATH="/u/project/ngarud/Garud_lab/Brendan/Utils/sims_to_numpy/helper_processSLiMsims.py"
MANAGER_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/sims_to_numpy/SLiMsims_to_numpy.py"
VISUALIZER="/u/project/ngarud/Garud_lab/Brendan/Utils/visualize_numpy_images/visualize_haplotypes.py"
TRAIN_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/train_color_CNN/color_CNN_train.py"

BASE_INPUT="/u/home/b/baeria/project-ngarud/wide_hmp_SLiMulations/dann_slimulations_results"

# --- Configuration: Parameters ---
NUM_SAMPS=120
WINDOW_SIZE=200
CHANNELS=2
MODEL_NAME="wide_hmp_color_cnn_model" # Update this to your preferred model name
BATCH_SIZE=32
TRAIN_PROP=1.0
TEST_PROP=0.2

# --- Part 1: Generate Numpy Arrays (Hard, Soft, Neutral) ---
for TYPE in hard soft neutral
do
    INPUT_FOLDER="${BASE_INPUT}/${TYPE}"
    OUTPUT_NAME="${TYPE}_sorted_color.npy"

    echo "--- Processing $TYPE ---"
    python "$MANAGER_SCRIPT" \
        "$OUTPUT_NAME" \
        "$INPUT_FOLDER" \
        --num_samps $NUM_SAMPS \
        --window_size $WINDOW_SIZE \
        --channels $CHANNELS \
        --sort \
        --helper_path "$HELPER_PATH"

    # Optional: Brief validation check
    if [ ! -f "$OUTPUT_NAME" ]; then
        echo "ERROR: $OUTPUT_NAME failed to generate. Exiting."
        exit 1
    fi
done

# --- Part 2: Train the CNN ---
echo "--- Starting CNN Training ---"
# We explicitly pass the files generated in the loop above
python "$TRAIN_SCRIPT" \
    "hard_sorted_color.npy" \
    "soft_sorted_color.npy" \
    "neutral_sorted_color.npy" \
    "$MODEL_NAME" \
    --batch_size "$BATCH_SIZE" \
    --train_prop "$TRAIN_PROP" \
    --test_prop "$TEST_PROP"

echo "Job $JOB_ID ended on: $(hostname -s) at $(date)"