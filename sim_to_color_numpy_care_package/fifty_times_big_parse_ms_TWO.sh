#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -l h_rt=08:00:00,h_data=8G
#$ -t 1-3000
#$ -tc 100

# =============================================================================
# MS to CSV Conversion Job Array (50x Scale)
# Processes 3 regimes (neutral, hard, soft) across 1000 bins each.
# Each bin contains 50 replicates.
# Total tasks: 3000
# =============================================================================

# 1. Load the environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# 2. Define Paths 
# Adjust BASE_DIR if the "fifty_times_bigger" folder is in a different parent directory
CONVERT_SCRIPT="util_scripts/parse_ms_mallelic_minimal.py"
BASE_DIR="ms_slimulations_results"
WINDOW_SIZE=50000
NREP=50

# 3. Determine Regime and Rep Bin based on SGE_TASK_ID
# Tasks 1-1000: neutral, 1001-2000: hard, 2001-3000: soft
if [ $SGE_TASK_ID -le 1000 ]; then
    REGIME="neutral"
    REP_BIN=$SGE_TASK_ID
elif [ $SGE_TASK_ID -le 2000 ]; then
    REGIME="hard"
    REP_BIN=$(($SGE_TASK_ID - 1000))
else
    REGIME="soft"
    REP_BIN=$(($SGE_TASK_ID - 2000))
fi

INPUT_DIR="$BASE_DIR/$REGIME"

# 4. Loop through the 50 replicates for the specific bin
for i in $(seq 1 $NREP); do
    IN_FILE="$INPUT_DIR/rep_$REP_BIN.$i.MS"
    OUT_FILE="$INPUT_DIR/rep_$REP_BIN.${i}_haplotypes.csv"

    # Check if the MS file exists before attempting conversion
    if [ -f "$IN_FILE" ]; then
        echo "Processing: $IN_FILE"

        python "$CONVERT_SCRIPT" \
            -i "$IN_FILE" \
            -o "$OUT_FILE" \
            -w $WINDOW_SIZE \
            --fitnessEffects

    else
        echo "Warning: File $IN_FILE not found."
    fi
done

echo "Task $SGE_TASK_ID for $REGIME bin $REP_BIN complete."
