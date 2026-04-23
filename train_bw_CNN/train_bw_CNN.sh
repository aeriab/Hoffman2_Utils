#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=24:00:00,h_data=30G,gpu
#$ -pe shared 1

# --- Load Environment ---
# Ensure these match your specific cluster setup
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

echo "Job $JOB_ID started on: `hostname -s` at `date`"

# --- Hard-Coded Configurations ---
# Replace these strings with your actual absolute paths
HARD_NPY="/example_hard.npy"
NEUTRAL_NPY="/example_neutral.npy"
SOFT_NPY="/example_soft.npy"
MODEL_NAME="CNN_Example_Title"
BATCH_SIZE=32
TRAIN_PROP=1.0
TEST_PROP=0.2

# Path to your actual training logic script
TRAIN_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/train_bw_CNN/bw_CNN_train.py"

# --- Run Training ---
echo "Beginning training for model: $MODEL_NAME"
echo "Using data: $HARD_NPY, $SOFT_NPY, $NEUTRAL_NPY"

python "$TRAIN_SCRIPT" \
    "$HARD_NPY" \
    "$SOFT_NPY" \
    "$NEUTRAL_NPY" \
    "$MODEL_NAME" \
    --batch_size "$BATCH_SIZE" \
    --train_prop "$TRAIN_PROP" \
    --test_prop "$TEST_PROP"

echo "Job $JOB_ID ended on: `hostname -s` at `date`"