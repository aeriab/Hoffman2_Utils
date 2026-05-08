
import itertools
# // Args: S_IND, THETA, FILE_PATH
## parameters for the jobs
all_regime = ["neutral", "hard", "soft"]
all_rep_bin = range(1, 501)

size = len(all_regime) * len(all_rep_bin)
job_name = 'all_dann_slimulations'

job_script = """#### all_dann_slimulations.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=23:00:00,h_data=8G
#$ -tc 499
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 1s
# Email address to notify
# Notify when
#$ -t 1-{}:1

declare -a all_regime_array=({})
declare -a all_rep_bin_array=({})

# Get the parameters for this job
rand_uniform() {{
  awk -v min="$1" -v max="$2" 'BEGIN {{srand(); print min + (max - min) * rand()}}'
}}

current_index=$(($SGE_TASK_ID-1))
REGIME=${{all_regime_array[$current_index]}}
REP_BIN=${{all_rep_bin_array[$current_index]}}

CHR_LENGTH=50000
N=10000
NREP=2
SAMPLE=200

SLIM_SCRIPT="ms_dann_slim.txt"

# Set slim parameters

OUTPUT_DIR="ms_slimulations_results/$REGIME"
mkdir -p "$OUTPUT_DIR"

#load modules
. /u/local/Modules/default/init/modules.sh
module load python
module load slim/4.0.1
module load conda
module load R

rm -f "$OUTPUT_DIR/rep_$REP_BIN."*

# run the SLiM script multiple times for the specified number of replicates
for i in $(seq 1 $NREP); do
    # get params
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
    # run slim
    FILE_PATH="$OUTPUT_DIR/rep_$REP_BIN.$i"
    slim -d S=$S -d THETA=$THETA -d RHO=$RHO -d MU=$MU -d SWEEP_END=$SWEEP_END -d CHR_LENGTH=$CHR_LENGTH -d N=$N -d SWEEP_LOC=$SWEEP_LOC -d SAMPLE=$SAMPLE -d "FILE_PATH='$FILE_PATH'" $SLIM_SCRIPT

    # we can add more processing here if needed
done
"""

# Generate all combinations of the parameters
combinations = list(itertools.product(all_regime, all_rep_bin))

# failed_jobs = [("neutral", 290), ("neutral", 296)]
# failed_jobs.extend([("hard", i) for i in range(401, 500)])
# failed_jobs.extend([("soft", i) for i in range(1, 378)])

# combinations = [("hard", 1), ("hard", 2), 
#                 ("neutral", 1), ("neutral", 2),
#                 ("soft", 1), ("soft", 2)]


# Generate the parameter arrays for each combination
all_regime_array = [str(combination[0]) for combination in combinations]
all_rep_bin_array = [str(combination[1]) for combination in combinations]

# Generate the job script content
job_array_size = len(combinations)
job_script_content = job_script.format(
    job_array_size,
    ' '.join(all_regime_array),
    ' '.join(all_rep_bin_array),
    script=job_name + '.py',
    job_name=job_name
)

# Write the job script to a file
with open(job_name + '.sh', 'w') as file:
    file.write(job_script_content)


