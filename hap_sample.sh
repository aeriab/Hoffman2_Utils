#!/bin/bash

# Load Conda into the shell session
. /u/local/Apps/anaconda3/etc/profile.d/conda.sh
conda activate hap-env

INPUT_NPY=$1
PREFIX=$2

# Check if input file exists
if [ ! -f "$INPUT_NPY" ]; then
    echo "Error: File $INPUT_NPY not found."
    exit 1
fi

FINAL_BUNDLE="${PREFIX}_selected_10.npy"

# Perform selection and stacking in a single Python execution
python -c "
import numpy as np
import sys

input_path = sys.argv[1]
output_path = sys.argv[2]

# Load using mmap_mode to handle large files without filling RAM
data = np.load(input_path, mmap_mode='r')
total_images = data.shape[0]

# Calculate 10 evenly spaced indices
indices = np.linspace(0, total_images - 1, 10, dtype=int)
print(f'Selected indices: {indices}')

# Extract slices and stack. 
# We remove the trailing channel dimension (..., 0) if it exists
# to get the requested (10, samples, sites) shape.
selected_data = data[indices, :, :, 0]

# Save the final bundle
np.save(output_path, selected_data)
print(f'Successfully saved array of shape {selected_data.shape} to {output_path}')
" "$INPUT_NPY" "$FINAL_BUNDLE"

echo "Done."