import itertools
import os

# // Parameters for the jobs
all_regime = ["neutral", "hard", "soft"]
all_rep_bin = range(1, 501)
job_name = 'dann_slimulations'

# // SGE Shell Script Template
# // Note: We use 8G of data because we are sampling 200 individuals inside SLiM,
# // which is much more memory-efficient than loading 30,000 later.
job_script = """#### {job_name}.sh START ####
#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -l h_rt=23:00:00,h_data=8G
#$ -tc 499
#$ -pe shared 1s
#$ -t 1-{job_array_size}:1

declare -a all_regime_array=({regimes})
declare -a all_rep_bin_array=({reps})

# Helper for random parameter generation
rand_uniform() {{
  awk -v min="$1" -v max="$2" 'BEGIN {{srand(); print min + (max - min) * rand()}}'
}}

current_index=$(($SGE_TASK_ID-1))
REGIME=${{all_regime_array[$current_index]}}
REP_BIN=${{all_rep_bin_array[$current_index]}}

# Global Simulation Constants
CHR_LENGTH=100000
N=10000
NREP=101
SAMPLE=200 

SLIM_SCRIPT="000001_dann_slim.txt"

# Set up output directory
OUTPUT_DIR="dann_slimulations_$JOB_ID/$REGIME"
mkdir -p $OUTPUT_DIR

# Load necessary modules
. /u/local/Modules/default/init/modules.sh
module load python/3.9.6
module load slim/4.0.1

# Run the SLiM script and process into NumPy
for i in $(seq 1 $NREP); do
    # Define parameters based on regime
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
    SWEEP_LOC=$(rand_uniform 0.25 0.75)
    
    FILE_BASE="$OUTPUT_DIR/rep_$REP_BIN.$i"
    
    # 1. Run SLiM - Outputs a .ms file using Peter's writeMSwithS function
    slim -d S=$S -d THETA=$THETA -d RHO=$RHO -d MU=$MU \\
         -d SWEEP_END=$SWEEP_END -d CHR_LENGTH=$CHR_LENGTH \\
         -d N=$N -d SWEEP_LOC=$SWEEP_LOC -d SAMPLE=$SAMPLE \\
         -d "FILE_PATH='$FILE_BASE.ms'" $SLIM_SCRIPT

    # 2. Process MS into 2-Channel, Frequency-Sorted NumPy
    # We call the conversion logic appended to this generator
    python -c "
import sys
import numpy as np
import pandas as pd

def parse_and_save(infile, outfile, window_size):
    # Parse MS with selection coefficients (Skip 3 rows for Peter's header)
    MS = pd.read_csv(infile, skiprows=3, header=None)
    
    # Get Selection Coeffs for Channel 2
    sel_coeffs = list(MS.iloc[1, :])[0].split()[1:]
    sel_coeffs = np.array(sel_coeffs, dtype=float)
    
    # Get Positions and resolve multi-allelic collisions
    pos = list(MS.iloc[0, :])[0].split()[1:]
    pos = np.round(np.array(pos, dtype=float) * window_size)
    
    # Process Haplotypes (Transposed: Rows=Samples, Cols=SNPs)
    haps = MS.iloc[2:, 0].apply(lambda x: pd.Series(list(x))).astype(int)
    
    # Frequency Sorting (Rows)
    hap_strings = haps.astype(str).apply(lambda x: ''.join(x), axis=1)
    counts = hap_strings.value_counts()
    new_order = []
    for h_str in counts.index:
        new_order.extend(hap_strings[hap_strings == h_str].index.tolist())
    haps = haps.loc[new_order].reset_index(drop=True)

    # Encode Channel 1: 0 -> -1 (Major), 1 -> 1 (Minor)
    # (Note: In MS, 0 is the 'reference' which we treat as Major)
    ch1 = haps.replace({{0: -1}}).values
    
    # Encode Channel 2: Syn/Nonsyn
    # -1 = Syn (s=0), 1 = Nonsyn (s!=0), 0 = Major/Missing
    ch2 = np.zeros_like(ch1)
    for col_idx in range(ch2.shape[1]):
        s = sel_coeffs[col_idx]
        val = 1 if s != 0 else -1
        # Only color where the minor allele (1) exists
        ch2[ch1[:, col_idx] == 1, col_idx] = val
        
    # Stack into (Samples, SNPs, 2)
    final_img = np.stack([ch1, ch2], axis=-1).astype(np.int8)
    np.save(outfile, final_img)

parse_and_save('$FILE_BASE.ms', '$FILE_BASE.npy', $CHR_LENGTH)
"
    # Cleanup the intermediate MS file
    rm "$FILE_BASE.ms"
done
"""

# // Generate Combinations
combinations = list(itertools.product(all_regime, all_rep_bin))
all_regime_array = [str(c[0]) for c in combinations]
all_rep_bin_array = [str(c[1]) for c in combinations]

# // Format and Write File
formatted_script = job_script.format(
    job_name=job_name,
    job_array_size=len(combinations),
    regimes=' '.join(all_regime_array),
    reps=' '.join(all_rep_bin_array)
)

with open(f"{job_name}.sh", 'w') as f:
    f.write(formatted_script)

print(f"Successfully generated {job_name}.sh with {len(combinations)} tasks.")