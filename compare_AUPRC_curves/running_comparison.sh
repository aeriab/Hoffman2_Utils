#!/bin/bash
#$ -cwd
#$ -o joblog_compare.$JOB_ID
#$ -j y
#$ -l h_rt=08:00:00,h_data=20G
#$ -pe shared 1

. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate tf_A100_clean

# Hard-coded Paths
HARD="/u/home/b/baeria/project-ngarud/Research/WideSimsCNN/har_remajored_middle_sorted.npy"
NEU="/u/home/b/baeria/project-ngarud/Research/WideSimsCNN/neu_remajored_middle_sorted.npy"
SOFT="/u/home/b/baeria/project-ngarud/Research/WideSimsCNN/sof_remajored_middle_sorted.npy"


# Color Model Files (epoch 12 had lowest validation loss)
C_JSON="/u/home/b/baeria/project-ngarud/Research/WideSimsCNN/Apr_25_adjusted_code_color_CPU_CNN/Color_CNN_354Win_120Samp_200MiddleSort_model.json"
C_WGHTS="/u/home/b/baeria/project-ngarud/Research/WideSimsCNN/Apr_25_adjusted_code_color_CPU_CNN/Color_CNN_354Win_120Samp_200MiddleSort.12.weights.h5"

# B&W Model Files (epoch 18 had lowest validation loss)
BW_JSON="/u/home/b/baeria/project-ngarud/Research/WideSimsCNN/middle_sort_bw_no_gpu_CNN/BW_CNN_354Win_120Samp_200MiddleSort_model.json"
BW_WGHTS="/u/home/b/baeria/project-ngarud/Research/WideSimsCNN/middle_sort_bw_no_gpu_CNN/BW_CNN_354Win_120Samp_200MiddleSort.18.weights.h5"

python /u/project/ngarud/Garud_lab/Brendan/Utils/compare_AUPRC_curves/compare_color_to_bw_AUPRC.py \
    --color_json "$C_JSON" --color_weights "$C_WGHTS" \
    --bw_json "$BW_JSON" --bw_weights "$BW_WGHTS" \
    --hard_npy "$HARD" --soft_npy "$SOFT" --neutral_npy "$NEU" \
    --test_prop 0.2