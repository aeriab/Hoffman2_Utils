#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -l h_rt=04:00:00,h_data=8G
#$ -t 1-1500

# =============================================================================
# MS to CSV Conversion Job Array
# Processes 3 regimes (neutral, hard, soft) across 500 bins each.
# Total tasks: 1500
# =============================================================================

# 1. Load the environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate base
conda activate tf_A100_clean

# 2. Define Hard-coded Paths
CONVERT_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/parse_ms_mallelic_minimal.py"
BASE_DIR="/u/home/b/baeria/project-ngarud/wide_hmp_SLiMulations/middle_mutation_ms_output/ms_slimulations_results"
WINDOW_SIZE=50000

# 3. Determine Regime and Rep Bin based on SGE_TASK_ID
# Tasks 1-500: neutral, 501-1000: hard, 1001-1500: soft
if [ $SGE_TASK_ID -le 500 ]; then
    REGIME="neutral"
    REP_BIN=$SGE_TASK_ID
elif [ $SGE_TASK_ID -le 1000 ]; then
    REGIME="hard"
    REP_BIN=$(($SGE_TASK_ID - 500))
else
    REGIME="soft"
    REP_BIN=$(($SGE_TASK_ID - 1000))
fi

INPUT_DIR="$BASE_DIR/$REGIME"

# 4. Loop through the two replicates (.1 and .2) for the specific bin
for i in 1 2; do
    IN_FILE="$INPUT_DIR/rep_$REP_BIN.$i"
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