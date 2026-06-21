#!/bin/bash
#$ -cwd
#$ -o joblog_all_sort.$JOB_ID
#$ -j y
#$ -l h_rt=08:00:00,h_data=16G
#$ -m abe

# =============================================================================
# Haplotype Processing: Major Allele Flipping & Middle-Window Frequency Sorting
# Explicitly processing Hard, Neutral, and Soft datasets.
# =============================================================================

# 1. Load the environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate tf_A100_clean

# 2. Define Script Path (Now assumed to be in the current directory)
SORT_SCRIPT="util_scripts/remajor_and_middle_sort.py"


# IN_* files:
IN_HARD="hard_sims.npy"
IN_NEUT="neutral_sims.npy"
IN_SOFT="soft_sims.npy"

# -----------------------------------------------------------------------------
# 3. Process HARD Regime
# -----------------------------------------------------------------------------
OUT_HARD="hard_flipped_mid50sorted.npy"

if [ -f "$IN_HARD" ]; then
    echo "Starting Hard Regime: $IN_HARD"
    python "$SORT_SCRIPT" "$IN_HARD" "$OUT_HARD"
else
    echo "Warning: $IN_HARD not found."
fi

# -----------------------------------------------------------------------------
# 4. Process NEUTRAL Regime
# -----------------------------------------------------------------------------
OUT_NEUT="neutral_flipped_mid50sorted.npy"

if [ -f "$IN_NEUT" ]; then
    echo "Starting Neutral Regime: $IN_NEUT"
    python "$SORT_SCRIPT" "$IN_NEUT" "$OUT_NEUT"
else
    echo "Warning: $IN_NEUT not found."
fi

# -----------------------------------------------------------------------------
# 5. Process SOFT Regime
# -----------------------------------------------------------------------------
OUT_SOFT="soft_flipped_mid50sorted.npy"

if [ -f "$IN_SOFT" ]; then
    echo "Starting Soft Regime: $IN_SOFT"
    python "$SORT_SCRIPT" "$IN_SOFT" "$OUT_SOFT"
else
    echo "Warning: $IN_SOFT not found."
fi

echo "----------------------------------------------------------------"
echo "All processing tasks complete."
