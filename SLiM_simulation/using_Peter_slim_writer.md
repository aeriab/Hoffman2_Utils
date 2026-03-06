SLiM to NumPy Pipeline: User Guide

Garud Lab - Brendan Aeria

This pipeline automates the generation of SLiM simulations and their immediate conversion into 2-channel, frequency-sorted NumPy arrays (.npy) for use in CNN training.1. Setup: Moving Files to your Working DirectoryTo avoid modifying the master scripts in the Utils folder, you should copy them to your specific scratch or project directory before running.From anywhere on Hoffman2, run these commands:Bash# 1. Create your working directory (change 'my_sims' to your preferred name)
mkdir -p ~/project-ngarud/my_sims
cd ~/project-ngarud/my_sims

# 2. Copy the simulation engine and the job generator
cp /u/project/ngarud/Garud_lab/Brendan/Utils/SLiM_simulation/Peter_slim_writer.py .
cp /u/project/ngarud/Garud_lab/Brendan/Utils/SLiM_simulation/000001_dann_slim.txt .
2. Execution WorkflowThe process follows a "Generator -> Scheduler -> Worker" flow.Step 1: Generate the Shell ScriptRun the Python "Writer" script. This will create a file named dann_slimulations.sh in your current folder.Bashpython Peter_slim_writer.py
Step 2: Submit to Hoffman2Submit the generated array job to the cluster. This will launch 1,500 tasks (500 replicates for Neutral, Hard, and Soft sweeps).Bashqsub dann_slimulations.sh
Step 3: Monitor JobsCheck the status of your tasks:Bashqstat -u $USER
3. Key Parameters for CustomizationIf you want to change the biological or computational scale of the simulations, look for these variables in the scripts:In Peter_slim_writer.py (The Manager)ParameterDefaultImpactall_rep_binrange(1, 501)Determines the total number of independent "seeds."NREP101How many simulations to run per task. (Total sims = Reps * NREP).SAMPLE200The number of haplotypes (rows) in your final .npy image.h_data8GMemory allocation. Increase if CHR_LENGTH is > 200kb.In 000001_dann_slim.txt (The Engine)ParameterLogicNoteCHR_LENGTH100000The physical length of the simulated DNA segment.N10000Population size. Affects genetic drift and time to fixation.S_DEL0.001The selection coefficient for deleterious mutations (m2).SWEEP_ENDVariableThe frequency the beneficial mutation must reach to trigger a "success" save.4. Output Data StructureOnce the jobs finish, your data will be organized as follows:Plaintextdann_slimulations_[JOB_ID]/
├── neutral/
│   ├── rep_1.1.npy  <-- (Samples, SNPs, 2)
│   └── ...
├── hard/
└── soft/
Channel Encoding:Channel 1: -1 (Major Allele), 1 (Minor Allele).Channel 2: -1 (Synonymous), 1 (Non-synonymous), 0 (Major/Missing).