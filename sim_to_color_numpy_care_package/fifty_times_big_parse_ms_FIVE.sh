#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -l h_rt=08:00:00,h_data=8G
#$ -t 1-3000
#$ -tc 100

# =============================================================================
# MS to CSV Conversion Job Array (50x Scale)
# Processes 3 regimes (neutral, hard, soft) across REP_BIN_MAX-1 bins each.
# Each bin contains NREP replicates.
#
# NOTE on the #$ -t 1-3000 default above: that's only used if this script is
# qsub'd standalone, where it matches the default REP_BIN_MAX=1001 (1000 bins
# x 3 regimes = 3000). When launched from run_full_pipeline.sh, the orchestrator
# overrides this with `qsub -t 1-$TOTAL_TASKS ...` computed from the actual
# REP_BIN_MAX, so the two always stay in sync.
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
WINDOW_SIZE=${WINDOW_SIZE:-50000}
NREP=${NREP:-50}
REP_BIN_MAX=${REP_BIN_MAX:-1001}

# REP_BIN_MAX is exclusive, matching three_bash_slim_scripts_generator_ONE.py's
# `range(1, REP_BIN_MAX)`, so each regime has (REP_BIN_MAX - 1) rep bins.
BINS_PER_REGIME=$((REP_BIN_MAX - 1))

# 3. Determine Regime and Rep Bin based on SGE_TASK_ID
# Tasks 1..BINS_PER_REGIME: neutral
# Tasks BINS_PER_REGIME+1..2*BINS_PER_REGIME: hard
# Tasks 2*BINS_PER_REGIME+1..3*BINS_PER_REGIME: soft
if [ $SGE_TASK_ID -le $BINS_PER_REGIME ]; then
    REGIME="neutral"
    REP_BIN=$SGE_TASK_ID
elif [ $SGE_TASK_ID -le $((BINS_PER_REGIME * 2)) ]; then
    REGIME="hard"
    REP_BIN=$(($SGE_TASK_ID - BINS_PER_REGIME))
else
    REGIME="soft"
    REP_BIN=$(($SGE_TASK_ID - BINS_PER_REGIME * 2))
fi

INPUT_DIR="$BASE_DIR/$REGIME"

# 4. Loop through the replicates for the specific bin
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
