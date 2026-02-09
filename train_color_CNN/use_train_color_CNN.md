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
qsub train_color_CNN.sh <hard_npy> <soft_npy> <neutral_npy> <output_model_name> [batch_size] [train_prop] [test_prop]
```

| Argument             | Required | Default | Description                                                        |
|----------------------|----------|---------|--------------------------------------------------------------------|
| `hard_npy`           | Yes      | —       | Path to Hard sweep `.npy` file                                     |
| `soft_npy`           | Yes      | —       | Path to Soft sweep `.npy` file                                     |
| `neutral_npy`        | Yes      | —       | Path to Neutral `.npy` file                                        |
| `output_model_name`  | Yes      | —       | Base name for output model files                                   |
| `batch_size`         | No       | 32      | Training batch size                                                |
| `train_prop`         | No       | 1.0     | Proportion of data to use (e.g., 0.8 for 80%)                     |
| `test_prop`          | No       | 0.2     | Proportion of used data held out for testing (0.0 to skip testing) |

### Example
```bash
# Default: uses all data, holds out 20% for test evaluation
qsub train_color_CNN.sh \
    /path/to/hard_sorted_color.npy \
    /path/to/soft_sorted_color.npy \
    /path/to/neutral_sorted_color.npy \
    my_color_model

# Custom: use 50% of data, batch size 64, 25% test holdout
qsub train_color_CNN.sh \
    /path/to/hard_sorted_color.npy \
    /path/to/soft_sorted_color.npy \
    /path/to/neutral_sorted_color.npy \
    my_color_model 64 0.5 0.25

# Skip test evaluation entirely
qsub train_color_CNN.sh \
    /path/to/hard_sorted_color.npy \
    /path/to/soft_sorted_color.npy \
    /path/to/neutral_sorted_color.npy \
    my_color_model 32 1.0 0.0
```

### 3. Data splitting
The data is split as follows:
```
Full .npy array
  └─ train_prop slice (front portion of full array)
       ├─ Training+Validation (front 1-test_prop portion)
       │    ├─ Training data (front 90%)
       │    └─ Validation data (back 10%, used during training for early stopping)
       └─ Test holdout (back test_prop portion — never seen during training)
```

With defaults (`train_prop=1.0`, `test_prop=0.2`), 80% of each class is used for training+validation and 20% is held out for final evaluation.

### 4. Output files
After training completes (30 epochs), you will find:

**Training outputs:**

| File                                  | Description                              |
|---------------------------------------|------------------------------------------|
| `my_color_model_model.json`           | Model architecture                       |
| `my_color_model.01.weights.h5` – `.30.weights.h5` | Per-epoch weight checkpoints |
| `my_color_model.weights.h5`           | Final epoch weights                      |
| `training_multiclass_results.txt`     | Per-epoch loss and accuracy values       |
| `classifier_accuracy_hist.png`        | Training/validation accuracy plot        |
| `classifier_loss_hist.png`            | Training/validation loss plot            |
| `joblog.$JOB_ID`                      | SGE job log                              |

**Test evaluation outputs** (when `test_prop > 0`):

| File                                       | Description                                                                 |
|--------------------------------------------|-----------------------------------------------------------------------------|
| `my_color_model_test_auprc.png`            | Precision-recall curves for all 3 classes + combined sweeps                 |
| `my_color_model_test_confusion_matrix.png` | Side-by-side raw count and normalized confusion matrices                    |
| `my_color_model_test_results.txt`          | Test loss, accuracy, per-class accuracy, and AUPRC values                   |

The AUPRC plot includes four curves: Neutral, Hard Sweep, Soft Sweep, and a combined Sweeps class (hard OR soft), matching the evaluation format used in the DANN pipeline.

### 5. Choosing the best epoch
Open `training_multiclass_results.txt` or inspect the accuracy/loss plots. Pick the epoch with the best validation accuracy (or lowest validation loss), then load that checkpoint:
```python
from CNN_multiclass_data_mergeSims_A100 import *

model = create_model(np.load("neutral_sorted_color.npy", mmap_mode='r'))
model.load_weights("my_color_model.15.weights.h5")  # example: epoch 15
```

Note: the test evaluation at the end of training uses the final epoch (epoch 30) weights. If a different epoch performs better on validation, you may want to re-run evaluation manually with those weights.

### Internal details
The job script calls the training Python script located at:
```
/u/project/ngarud/Garud_lab/Brendan/Utils/train_color_CNN/color_CNN_train.py
```
which imports utilities from:
```
/u/project/ngarud/Garud_lab/DANN/Utils/CNN_multiclass_data_mergeSims_A100.py
```