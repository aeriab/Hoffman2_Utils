#!/bin/bash

# Load Conda into the shell session
. /u/local/Apps/anaconda3/etc/profile.d/conda.sh
conda activate hap-env

INPUT_NPY=$1
PREFIX=$2
CHANNEL=${3:-0} # Default to channel 0 (binary) if not provided

# Check if input file exists
if [ ! -f "$INPUT_NPY" ]; then
    echo "Error: File $INPUT_NPY not found."
    exit 1
fi

# Define suffix based on channel
if [ "$CHANNEL" -eq "0" ]; then
    SUFFIX="binary"
else
    SUFFIX="color"
fi

FINAL_BUNDLE="${PREFIX}_${SUFFIX}_selected_10.npy"

# Perform selection and stacking in a single Python execution
python -c "
import numpy as np
import sys

input_path = sys.argv[1]
output_path = sys.argv[2]
chan_idx = int(sys.argv[3])

# Load using mmap_mode for memory efficiency
data = np.load(input_path, mmap_mode='r')
total_images = data.shape[0]

# Calculate 10 evenly spaced indices
indices = np.linspace(0, total_images - 1, 10, dtype=int)
print(f'Selected indices: {indices}')

# Extracting based on the flexible channel argument
# Final shape: (10, samples, sites)
selected_data = data[indices, :, :, chan_idx]

# Save the final bundle
np.save(output_path, selected_data)
print(f'Successfully saved channel {chan_idx} ({sys.argv[4]}) data.')
print(f'Array shape: {selected_data.shape} -> {output_path}')
" "$INPUT_NPY" "$FINAL_BUNDLE" "$CHANNEL" "$SUFFIX"

echo "Done."