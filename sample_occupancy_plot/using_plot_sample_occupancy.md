# A. finegoldii Occupancy Plotter

This script analyzes the processed `A_finegoldii_XX.csv` files to visualize how many genomic sites are shared across your samples, specifically filtering for high-quality data.

## Prerequisites
- **Python 3.x**
- **pandas** and **matplotlib** libraries (`pip install pandas matplotlib`)
- **Files needed in directory**:
  - `Alistipes_finegoldii_56071_good_samples.txt`
  - Processed CSVs (`A_finegoldii_01.csv` through `A_finegoldii_86.csv`)

## How to Use
1. Place the `.py` script in the folder containing your processed CSV files and `.txt` sample list.
2. Run the script:
   ```bash
   python /u/project/ngarud/Garud_lab/Brendan/Utils/sample_occupancy_plot/plot_sample_occupancy.py