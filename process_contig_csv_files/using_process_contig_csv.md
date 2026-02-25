# A. finegoldii Contig Processor

This script filters and renames 86 contig CSV files based on a "good samples" list and sorts them by the number of sites (rows).

## Prerequisites
- **Python 3.x**
- **pandas** library (`pip install pandas`)
- **Files needed in directory**:
  - `Alistipes_finegoldii_56071_good_samples.txt`
  - Input CSVs (`AENZ01000001_haplotypes.csv` to `AENZ01000086_haplotypes.csv`)

## How to Use
1. Place the `.py` script in the folder containing your CSV files and `.txt` sample list.
2. Run the script:
   ```bash
   python /u/project/ngarud/Garud_lab/Brendan/Utils/process_contig_csv_files/process_contig_csv.py