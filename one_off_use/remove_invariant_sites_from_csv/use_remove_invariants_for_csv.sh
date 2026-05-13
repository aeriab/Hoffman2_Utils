#!/bin/bash

# --- CONFIGURATION ---
# Path to the directory containing your CSV data
DATA_DIR="/path/to/your/data/folder"

# Filenames
INPUT_NAME="your_input_file.csv"
OUTPUT_NAME="filtered_variant_sites.csv"

# Path to your existing utility script
PY_SCRIPT="/u/project/ngarud/Garud_lab/Brendan/Utils/one_off_use/remove_invariant_sites_from_csv/remove_invariants_for_csv.py"
# ---------------------

# Change to the data directory
cd "$DATA_DIR" || { echo "Directory not found"; exit 1; }

echo "Starting invariant site removal..."
echo "Input: $DATA_DIR/$INPUT_NAME"

# Execute the Python script using the arguments defined above
python3 "$PY_SCRIPT" "$INPUT_NAME" "$OUTPUT_NAME"

echo "Done. Cleaned file saved as: $OUTPUT_NAME"