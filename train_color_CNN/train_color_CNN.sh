#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=24:00:00,h_data=30G,gpu
#$ -pe shared 1

# --- Load Environment ---
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

echo "Job $JOB_ID started on: `hostname -s` at `date`"

# --- User Arguments ---
# $1 = hard sweep npy
# $2 = soft sweep npy
# $3 = neutral npy
# $4 = output model name
# $5 = batch size (optional, default 32)
# $6 = train proportion (optional, default 1.0)
# $7 = test proportion (optional, default 0.2)

HARD_NPY="$1"
SOFT_NPY="$2"
NEUTRAL_NPY="$3"
MODEL_NAME="$4"
BATCH_SIZE="${5:-32}"
TRAIN_PROP="${6:-1.0}"
TEST_PROP="${7:-0.2}"

TRAIN_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/train_color_CNN/color_CNN_train.py"

python "$TRAIN_SCRIPT" \
    "$HARD_NPY" \
    "$SOFT_NPY" \
    "$NEUTRAL_NPY" \
    "$MODEL_NAME" \
    --batch_size "$BATCH_SIZE" \
    --train_prop "$TRAIN_PROP" \
    --test_prop "$TEST_PROP"

echo "Job $JOB_ID ended on: `hostname -s` at `date`"

# You'd submit it like:
#### bash:
# qsub train_cnn.sh hard_sorted_color.npy soft_sorted_color.npy neutral_sorted_color.npy my_color_model
# qsub train_cnn.sh hard.npy soft.npy neutral.npy my_model 32 1.0 0.2