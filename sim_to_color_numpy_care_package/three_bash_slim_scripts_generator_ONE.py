import itertools
import os

# Parameters for the jobs
all_regime = ["neutral", "hard", "soft"]
# all_rep_bin = range(1, 1001)
all_rep_bin = range(1, int(os.environ.get("REP_BIN_MAX", 1001)))


# Simulation parameters - read from environment (set by run_full_pipeline.sh),
# falling back to these defaults if run standalone.
CHR_LENGTH = int(os.environ.get("CHR_LENGTH", 50000))
N = int(os.environ.get("N", 10000))
NREP = int(os.environ.get("NREP", 50))
SAMPLE = int(os.environ.get("SAMPLE", 200))

# The template for the .sh files
# Note: I changed -tc to 100 to respect your per-job cap
job_script_template = """#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -l h_rt=23:00:00,h_data=8G
#$ -tc 100
#$ -pe shared 1
#$ -t 1-{task_limit}

declare -a all_regime_array=({regime_list})
declare -a all_rep_bin_array=({rep_list})

# Get the parameters for this job
rand_uniform() {{
  awk -v min="$1" -v max="$2" 'BEGIN {{srand(); print min + (max - min) * rand()}}'
}}

current_index=$(($SGE_TASK_ID-1))
REGIME=${{all_regime_array[$current_index]}}
REP_BIN=${{all_rep_bin_array[$current_index]}}

CHR_LENGTH={chr_length}
N={n}
NREP={nrep}
SAMPLE={sample}

SLIM_SCRIPT="util_scripts/ms_dann_slim.txt"
OUTPUT_DIR="ms_slimulations_results/$REGIME"
mkdir -p "$OUTPUT_DIR"

# Load modules
. /u/local/Modules/default/init/modules.sh
module load python
module load slim/4.0.1
module load conda
module load R

rm -f "$OUTPUT_DIR/rep_$REP_BIN."*

# Run replicates
for i in $(seq 1 $NREP); do
    if [ "$REGIME" == "neutral" ]; then
        S=0
        THETA=$(rand_uniform 0.01 5)
    elif [ "$REGIME" == "soft" ]; then
        S=$(rand_uniform 0.01 0.1)
        THETA=$(rand_uniform 1 5)
    else
        S=$(rand_uniform 0.01 0.1)
        THETA=$(rand_uniform 0.001 0.01)
    fi

    RHO=$(rand_uniform 1e-8 1e-7)
    MU=$(rand_uniform 1e-7 1e-6)
    SWEEP_END=$(rand_uniform 0.2 0.9)
    SWEEP_LOC=0.5
    FILE_PATH="$OUTPUT_DIR/rep_$REP_BIN.$i"
    
    slim -d S=$S -d THETA=$THETA -d RHO=$RHO -d MU=$MU -d SWEEP_END=$SWEEP_END -d CHR_LENGTH=$CHR_LENGTH -d N=$N -d SWEEP_LOC=$SWEEP_LOC -d SAMPLE=$SAMPLE -d "FILE_PATH='$FILE_PATH'" $SLIM_SCRIPT
done
"""

# Iterate through each regime to create a specific script for it
for regime in all_regime:
    # Create the lists for just this specific regime
    regime_list = " ".join([regime] * len(all_rep_bin))
    rep_list = " ".join([str(r) for r in all_rep_bin])
    
    # Format the template
    formatted_script = job_script_template.format(
        task_limit=len(all_rep_bin),
        regime_list=regime_list,
        rep_list=rep_list,
        chr_length=CHR_LENGTH,
        n=N,
        nrep=NREP,
        sample=SAMPLE
    )
    
    # Write to a unique filename
    filename = f"run_slim_{regime}.sh"
    with open(filename, 'w') as f:
        f.write(formatted_script)
    
    print(f"Generated {filename}")
