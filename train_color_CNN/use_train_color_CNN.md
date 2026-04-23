## How to train the Color CNN (3-class: Neutral, Hard Sweep, Soft Sweep)

Each `.npy` should have shape `(num_sims, num_samples, window_h, 2)`.

---

### 1. Copy the job script to your working directory
```bash
cp /u/project/ngarud/Garud_lab/Brendan/Utils/train_color_CNN/train_color_CNN.sh .
```

## MAKE SURE TO MODIFY THE HARD CODED VARIABLES IN train_color_CNN

### 2. Submit the job
```bash
qsub train_color_CNN.sh
```
