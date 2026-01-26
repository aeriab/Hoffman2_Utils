import numpy as np
import subprocess
import glob
import os
import sys

# Usage: python processSLiMsims_ManagerScript.py <out.npy> <n_samps> <window> <in_dir> <channels>

# 1. Parse Arguments
OUTPUT_FILE = sys.argv[1] 
NUM_SAMPS   = int(sys.argv[2]) 
WINDOW_SIZE = int(sys.argv[3]) 
INPUT_DIR   = sys.argv[4]
CHANNELS    = int(sys.argv[5]) # 1 for (N, S, W), 2 for (N, S, W, 2)

# 2. Identify input files
input_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.txt")))
NUM_SIMS = len(input_files)

# 3. Define the shape based on requested channels
if CHANNELS == 2:
    SIM_SHAPE = (NUM_SAMPS, WINDOW_SIZE, 2)
else:
    SIM_SHAPE = (NUM_SAMPS, WINDOW_SIZE)

# 4. Preallocate the huge file on disk
# This creates the 'skeleton' that the workers will fill in
big_array = np.lib.format.open_memmap(
    OUTPUT_FILE, dtype=np.float32, mode="w+", shape=(NUM_SIMS,) + SIM_SHAPE
)
del big_array # Close it so workers can open it in r+ mode

# 5. Loop over files and launch the worker for each
for i, infile in enumerate(input_files):
    subprocess.run(
        [
            "python",
            "/u/project/ngarud/Garud_lab/Brendan/Utils/processSLiMsims.py",
            infile,
            OUTPUT_FILE,
            str(NUM_SAMPS),
            str(WINDOW_SIZE),
            str(i)
        ],
        check=True,
    )

print(f"Done! Processed {NUM_SIMS} files into {OUTPUT_FILE}")