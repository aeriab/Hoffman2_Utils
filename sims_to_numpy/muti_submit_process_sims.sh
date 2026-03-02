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
conda activate base
conda activate tf_A100_clean

# --- Configuration ---
HELPER_PATH="/u/project/ngarud/Garud_lab/Brendan/Utils/sims_to_numpy/helper_processSLiMsims.py"
MANAGER_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/sims_to_numpy/SLiMsims_to_numpy.py"
VISUALIZER="/u/project/ngarud/Garud_lab/Brendan/Utils/numpy_image_generator/visualize_haplotypes.py"

BASE_INPUT="/u/home/b/baeria/project-ngarud/wide_hmp_SLiMulations/dann_slimulations_results"
NUM_SAMPS=120
WINDOW_SIZE=200
CHANNELS=2

# --- Processing Loop ---
for TYPE in hard soft neutral
do
    INPUT_FOLDER="${BASE_INPUT}/${TYPE}"
    OUTPUT_NAME="${TYPE}_sorted_color.npy"

    echo "----------------------------------------------"
    echo "Processing Sweep Type: $TYPE"
    echo "  Input:  $INPUT_FOLDER"
    echo "  Output: $OUTPUT_NAME"
    echo "----------------------------------------------"

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
        echo "Generating validation plots for $TYPE..."
        mkdir -p "plots_${TYPE}"
        python "$VISUALIZER" \
            "$OUTPUT_NAME" \
            --validate \
            --n 5 \
            --grid \
            --prefix "${TYPE}_sweep" \
            --output-dir "./plots_${TYPE}"
    fi
done

echo "=============================================="
echo "All sweep types processed."
echo "=============================================="