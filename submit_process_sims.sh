#!/bin/bash

#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -m a
#$ -l h_rt=200:00:00,h_data=50G,highp
#$ -pe shared 1

# Load the job environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
# Ensuring the base is active before the specific clean env
conda activate base
conda activate tf_A100_clean

# Echo job info on joblog
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# --- Configuration ---
# Update these paths to point to your actual utility scripts and data
HELPER_PATH="/u/project/ngarud/Garud_lab/Brendan/Utils/helper_processSLiMsims.py"
MANAGER_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/SLiMsims_to_numpy.py"
INPUT_FOLDER="/u/home/b/baeria/project-ngarud/hmp_SLiMulations/dann_slimulations_12080244/hard/"
OUTPUT_NAME="hard_sorted_color.npy"

# --- Run Code ---
# We pass the manager script arguments: Output, Input, Samps, Window, Channels, Sort
python "$MANAGER_SCRIPT" \
    "$OUTPUT_NAME" \
    "$INPUT_FOLDER" \
    --num_samps 154 \
    --window_size 200 \
    --channels 2 \
    --sort \
    --helper_path "$HELPER_PATH"

# Echo job info on joblog
echo " "
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "