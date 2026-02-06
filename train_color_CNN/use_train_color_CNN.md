## How to train the Color CNN (3-class: Neutral, Hard Sweep, Soft Sweep)

### Prerequisites
You need three `.npy` files (one per class) already processed with the new color encoding scheme. See the [HMP csv to numpy README](../hmp_to_color_npy/HMP_csv_to_numpy_README.md) or equivalent SLiM pipeline for how to generate these.

Each `.npy` should have shape `(num_sims, num_samples, window_h, 2)`.

---

### 1. Copy the job script to your working directory
```bash
cp /u/project/ngarud/Garud_lab/Brendan/Utils/train_color_CNN/train_color_CNN.sh .
```

### 2. Submit the job
```bash
qsub train_color_CNN.sh <hard_npy> <soft_npy> <neutral_npy> <output_model_name> [batch_size] [train_prop]
```

| Argument             | Required | Default | Description                                      |
|----------------------|----------|---------|--------------------------------------------------|
| `hard_npy`           | Yes      | —       | Path to Hard sweep `.npy` file                   |
| `soft_npy`           | Yes      | —       | Path to Soft sweep `.npy` file                   |
| `neutral_npy`        | Yes      | —       | Path to Neutral `.npy` file                      |
| `output_model_name`  | Yes      | —       | Base name for output model files                 |
| `batch_size`         | No       | 32      | Training batch size                              |
| `train_prop`         | No       | 1.0     | Proportion of data to use (e.g., 0.8 for 80%)   |

### Example
```bash
qsub train_color_CNN.sh \
    /path/to/hard_sorted_color.npy \
    /path/to/soft_sorted_color.npy \
    /path/to/neutral_sorted_color.npy \
    my_color_model
```

### 3. Output files
After training completes (30 epochs), you will find:

| File                                  | Description                              |
|---------------------------------------|------------------------------------------|
| `my_color_model_model.json`           | Model architecture                       |
| `my_color_model.01.weights.h5` – `.30.weights.h5` | Per-epoch weight checkpoints |
| `my_color_model.weights.h5`           | Final epoch weights                      |
| `training_multiclass_results.txt`     | Per-epoch loss and accuracy values       |
| `classifier_accuracy_hist.png`        | Training/validation accuracy plot        |
| `classifier_loss_hist.png`            | Training/validation loss plot            |
| `joblog.$JOB_ID`                      | SGE job log                              |

### 4. Choosing the best epoch
Open `training_multiclass_results.txt` or inspect the accuracy/loss plots. Pick the epoch with the best validation accuracy (or lowest validation loss), then load that checkpoint:
```python
from CNN_multiclass_data_mergeSims_A100 import *

model = create_model(np.load("neutral_sorted_color.npy", mmap_mode='r'))
model.load_weights("my_color_model.15.weights.h5")  # example: epoch 15
```

### Internal details
The job script calls the training Python script located at:
```
/u/project/ngarud/Garud_lab/Brendan/Utils/train_color_CNN/color_CNN_train.py
```
which imports utilities from:
```
/u/project/ngarud/Garud_lab/DANN/Utils/CNN_multiclass_data_mergeSims_A100.py
```