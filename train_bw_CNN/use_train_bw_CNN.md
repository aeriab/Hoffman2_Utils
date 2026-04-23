## How to train the Color CNN (3-class: Neutral, Hard Sweep, Soft Sweep)

### Prerequisites
You need three `.npy` files (one per class) already processed with the new color encoding scheme. See the [HMP csv to numpy README](../hmp_to_color_npy/HMP_csv_to_numpy_README.md) or equivalent SLiM pipeline for how to generate these.

Each `.npy` should have shape `(num_sims, num_samples, window_h, 1)`.

---

### 1. Copy the job script to your working directory
```bash
cp /u/project/ngarud/Garud_lab/Brendan/Utils/train_bw_CNN/train_bw_CNN.sh .
```

## MAKE SURE TO MODIFY THE HARD CODED VARIABLES IN train_bw_CNN

### 2. Submit the job
```bash
qsub train_bw_CNN.sh
```
