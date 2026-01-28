#!/bin/bash

# Get the directory where THIS script is located
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
HAP_PLOT_PATH="$SCRIPT_DIR/hap_plot.R"

INPUT_NPY=$1
PREFIX=$2

# Check if file exists
if [ ! -f "$INPUT_NPY" ]; then
    echo "Error: File $INPUT_NPY not found."
    exit 1
fi

# Load R module inside script to be safe
unset R_HOME
module load R/4.2.2

TOTAL_IMAGES=$(python -c "import numpy as np; print(np.load('$INPUT_NPY', mmap_mode='r').shape[0])")
indices=$(python -c "import numpy as np; print(' '.join(map(str, np.linspace(0, $TOTAL_IMAGES-1, 10, dtype=int))))")

for i in $indices; do
    TEMP_NPY="temp_idx_${i}.npy"
    OUT_PNG="${PREFIX}_idx_${i}.png"

    # Extract slice
    python -c "import numpy as np, sys; \
               arr = np.load(sys.argv[1], mmap_mode='r'); \
               np.save(sys.argv[2], arr[int(sys.argv[3]), :, :, 0][np.newaxis, :, :])" \
               "$INPUT_NPY" "$TEMP_NPY" "$i"

    # Call R using the absolute path we found earlier
    Rscript "$HAP_PLOT_PATH" -f "$TEMP_NPY" -i 1 -o "$OUT_PNG" --sort_method 'frequency'

    rm "$TEMP_NPY"
done