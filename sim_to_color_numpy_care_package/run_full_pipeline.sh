#!/bin/bash
# =============================================================================
# run_full_pipeline.sh
#
# Orchestrates the full SLiM -> CSV -> numpy -> flipped/sorted pipeline.
# Each stage blocks until the previous one is FULLY complete (including all
# SGE array tasks) before starting the next, using `qsub -sync y`.
#
# Within a stage, the three regimes (neutral/hard/soft) run in parallel.
#
# IMPORTANT: This script itself takes a long time to run (potentially many
# hours to days, since it blocks on qsub -sync y for each stage). Do NOT run
# it in a bare SSH session that might disconnect. Use one of:
#   tmux new -s pipeline        (then run the script inside the tmux session)
#   screen -S pipeline
#   nohup ./run_full_pipeline.sh > pipeline_nohup.log 2>&1 &
# =============================================================================

set -uo pipefail

# =============================================================================
# PIPELINE PARAMETERS
# Edit these values to control every stage from one place. They're exported
# so Stage 1's Python script can read them directly, and passed into the
# qsub jobs for Stages 3 and 4 via `qsub -v`.
# =============================================================================

# Adjust these to adjust the haplotype image dimensions
export WINDOW_H=1500          # Number of SNPs (image height/rows of the SNP axis) kept from the centered middle window of each haplotype image.
export SAMPLE=420         # Number of haploid genomes sampled out of the SLiM population per replicate. This is the ceiling on how many haplotypes are available downstream — TARGET_SAMPLES (below) must be <= SAMPLE.
export TARGET_SAMPLES=220    # Number of haplotypes each image is downsampled to (image width). Must be <= SAMPLE above, since you can't downsample to more haplotypes than were originally sampled.


# --- Used by STAGE 1 (ONE: SLiM simulation generator) ---
export CHR_LENGTH=300000   # Length, in base pairs, of the simulated chromosome/locus. Bigger = more SNPs per sim, longer SLiM runtime.
export N=10000            # Diploid population size in the SLiM forward simulation. Bigger = more realistic drift/coalescent behavior, slower runtime.

# --- Shared by STAGE 1 (ONE) and STAGE 3 (FIVE) - must match between the two ---
export NREP=50            # Number of replicate simulations generated per (regime, bin) in Stage 1, and the number of replicates looped over per bin when parsing in Stage 3. Changing this changes how many total simulations exist per bin.
export REP_BIN_MAX=11     # Ussually 1001, but set lower for testing
# NOTE: REP_BIN_MAX is exclusive (Stage 1 uses range(1, REP_BIN_MAX)), so the
# actual number of rep bins generated PER REGIME is (REP_BIN_MAX - 1).
# Stage 3 below derives its array size and regime boundaries from this same
# value, so the two stages always agree on how many bins exist.

# --- Used by STAGE 3 (FIVE: MS -> CSV parsing) ---
export WINDOW_SIZE=50000  # Genomic window size, in base pairs, used when slicing MS output into per-window haplotype CSVs.

# --- Used by STAGE 4 (SIX/SEVEN/EIGHT: CSV -> numpy conversion) ---
export COMPLETE_THRESH=0.7   # Minimum fraction (0-1) of "complete" (non-missing) SNPs a window needs to be kept; windows below this threshold are discarded.

# Sanity check: catch an obviously bad config before burning cluster time
if (( TARGET_SAMPLES > SAMPLE )); then
    echo "CONFIG ERROR: TARGET_SAMPLES ($TARGET_SAMPLES) cannot be greater than SAMPLE ($SAMPLE)."
    echo "TARGET_SAMPLES is a downsample of SAMPLE, so it must be less than or equal to it."
    exit 1
fi

if (( REP_BIN_MAX < 2 )); then
    echo "CONFIG ERROR: REP_BIN_MAX ($REP_BIN_MAX) must be >= 2 (need at least 1 rep bin per regime)."
    exit 1
fi

# Number of rep bins per regime, derived once here so every stage that needs
# it (currently Stage 3) stays consistent with what Stage 1 generates.
BINS_PER_REGIME=$((REP_BIN_MAX - 1))

LOG="pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "=============================================================="
echo "PIPELINE PARAMETERS FOR THIS RUN"
echo "=============================================================="
echo "CHR_LENGTH      = $CHR_LENGTH"
echo "N               = $N"
echo "SAMPLE          = $SAMPLE"
echo "NREP            = $NREP"
echo "REP_BIN_MAX     = $REP_BIN_MAX  (=> $BINS_PER_REGIME rep bins per regime)"
echo "WINDOW_SIZE     = $WINDOW_SIZE"
echo "WINDOW_H        = $WINDOW_H"
echo "TARGET_SAMPLES  = $TARGET_SAMPLES"
echo "COMPLETE_THRESH = $COMPLETE_THRESH"
echo "=============================================================="

fail() {
    echo ""
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "PIPELINE FAILED: $1"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    exit 1
}

stage_header() {
    echo ""
    echo "=============================================================="
    echo "$1"
    echo "Started: $(date)"
    echo "=============================================================="
}

# -----------------------------------------------------------------------------
# STAGE 1 (ONE): Generate the per-regime SLiM array-job scripts
# -----------------------------------------------------------------------------
stage_header "STAGE 1: Generating SLiM SGE array scripts"

python three_bash_slim_scripts_generator_ONE.py || fail "Stage 1 (script generation) failed"

for f in run_slim_neutral.sh run_slim_hard.sh run_slim_soft.sh; do
    [ -f "$f" ] || fail "Expected generated script $f not found after Stage 1"
done
echo "Stage 1 complete: run_slim_neutral.sh, run_slim_hard.sh, run_slim_soft.sh generated."

# -----------------------------------------------------------------------------
# STAGE 2 (TWO/THREE/FOUR): Run SLiM simulations, all 3 regimes in parallel
# -----------------------------------------------------------------------------
stage_header "STAGE 2: SLiM simulations (neutral, hard, soft) - parallel"

qsub -sync y run_slim_neutral.sh & PID_NEUTRAL=$!
qsub -sync y run_slim_hard.sh    & PID_HARD=$!
qsub -sync y run_slim_soft.sh    & PID_SOFT=$!

STAGE2_FAIL=0
wait $PID_NEUTRAL || { echo "WARNING: neutral SLiM array job (qsub) reported nonzero exit"; STAGE2_FAIL=1; }
wait $PID_HARD    || { echo "WARNING: hard SLiM array job (qsub) reported nonzero exit"; STAGE2_FAIL=1; }
wait $PID_SOFT    || { echo "WARNING: soft SLiM array job (qsub) reported nonzero exit"; STAGE2_FAIL=1; }

for regime in neutral hard soft; do
    DIR="ms_slimulations_results/$regime"
    [ -d "$DIR" ] && [ -n "$(ls -A "$DIR" 2>/dev/null)" ] || fail "Stage 2: $DIR is missing or empty - SLiM run for $regime did not produce output"
done

[ "$STAGE2_FAIL" -eq 0 ] || echo "NOTE: qsub reported a nonzero exit for at least one regime, but output files exist. Spot-check joblog.* files before trusting these results."
echo "Stage 2 complete: SLiM simulations finished for all regimes."

# -----------------------------------------------------------------------------
# STAGE 3 (FIVE): Parse MS output to CSV (all regimes, single array job)
# -----------------------------------------------------------------------------
stage_header "STAGE 3: Parsing MS output to CSV"

# Stage 3's array job is laid out as [neutral bins][hard bins][soft bins], each
# block of size BINS_PER_REGIME (derived from REP_BIN_MAX above). The -t range
# is passed on the qsub command line because that's the only way to override
# the #$ -t default baked into the .sh file at submission time; REP_BIN_MAX
# itself is passed via -v so the script can compute the same regime boundaries
# internally.
TOTAL_TASKS=$((BINS_PER_REGIME * 3))

qsub -sync y -t 1-$TOTAL_TASKS \
    -v WINDOW_SIZE="$WINDOW_SIZE",NREP="$NREP",REP_BIN_MAX="$REP_BIN_MAX" \
    fifty_times_big_parse_ms_FIVE.sh || echo "WARNING: qsub reported nonzero exit for MS parsing job"

# Spot-check: at least one CSV should exist per regime
for regime in neutral hard soft; do
    CSV_COUNT=$(ls ms_slimulations_results/"$regime"/*_haplotypes.csv 2>/dev/null | wc -l)
    [ "$CSV_COUNT" -gt 0 ] || fail "Stage 3: no *_haplotypes.csv files found for $regime - MS parsing did not produce output"
done
echo "Stage 3 complete: MS files parsed to CSV for all regimes."

# -----------------------------------------------------------------------------
# STAGE 4 (SIX/SEVEN/EIGHT): Convert haplotype CSVs to numpy, parallel
# -----------------------------------------------------------------------------
stage_header "STAGE 4: Converting haplotype CSVs to numpy (hard, neutral, soft) - parallel"

QV="WINDOW_H=$WINDOW_H,TARGET_SAMPLES=$TARGET_SAMPLES,COMPLETE_THRESH=$COMPLETE_THRESH"
qsub -sync y -v "$QV" hard_haplotypes_to_numpy_SIX.sh      & PID_H=$!
qsub -sync y -v "$QV" neutral_haplotypes_to_numpy_SEVEN.sh & PID_N=$!
qsub -sync y -v "$QV" soft_haplotypes_to_numpy_EIGHT.sh    & PID_S=$!

STAGE4_FAIL=0
wait $PID_H || { echo "WARNING: hard numpy conversion job reported nonzero exit"; STAGE4_FAIL=1; }
wait $PID_N || { echo "WARNING: neutral numpy conversion job reported nonzero exit"; STAGE4_FAIL=1; }
wait $PID_S || { echo "WARNING: soft numpy conversion job reported nonzero exit"; STAGE4_FAIL=1; }

for f in hard_sims.npy neutral_sims.npy soft_sims.npy; do
    [ -f "$f" ] || fail "Stage 4: expected output $f not found"
done

[ "$STAGE4_FAIL" -eq 0 ] || echo "NOTE: qsub reported a nonzero exit for at least one regime, but output .npy files exist. Spot-check joblog.* files before trusting these results."
echo "Stage 4 complete: hard_sims.npy, neutral_sims.npy, soft_sims.npy created."

# -----------------------------------------------------------------------------
# STAGE 5 (NINE): Flip major allele + middle-window sort
# -----------------------------------------------------------------------------
stage_header "STAGE 5: Flipping major allele + middle-window sorting"

qsub -sync y raw_numpy_to_flipped_sorted_NINE.sh || echo "WARNING: qsub reported nonzero exit for flip/sort job"

for f in hard_flipped_mid50sorted.npy neutral_flipped_mid50sorted.npy soft_flipped_mid50sorted.npy; do
    [ -f "$f" ] || fail "Stage 5: expected output $f not found - check that raw_numpy_to_flipped_sorted_NINE.sh's input filenames (IN_HARD/IN_NEUT/IN_SOFT) match the actual files produced in Stage 4 (hard_sims.npy, neutral_sims.npy, soft_sims.npy)"
done

stage_header "PIPELINE COMPLETE"
echo "Final outputs: hard_flipped_mid50sorted.npy, neutral_flipped_mid50sorted.npy, soft_flipped_mid50sorted.npy"
