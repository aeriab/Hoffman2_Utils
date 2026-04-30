### --------- load modules -------------------#
import sys
import os
import argparse
import numpy as np

# Utils_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '/u/project/ngarud/Garud_lab/DANN/Utils/'))
Utils_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '/u/project/ngarud/Garud_lab/Brendan/Utils/DANN_Utils/'))

sys.path.append(Utils_path)

from CNN_zero_as_missing_A100 import *
import gc

parser = argparse.ArgumentParser()
parser.add_argument('hard_npy', type=str, help="Path to Hard sweep npy file")
parser.add_argument('soft_npy', type=str, help="Path to Soft sweep npy file")
parser.add_argument('neutral_npy', type=str, help="Path to Neutral sweep npy file")
parser.add_argument('output_name', type=str, help="Name of the output model")
parser.add_argument('--batch_size', type=int, default=32)
parser.add_argument('--train_prop', type=float, default=1.0, help="Proportion of total available data to use")
parser.add_argument('--test_prop', type=float, default=0.2, help="Proportion of used data to hold out for testing (0.0 to skip)")
args = parser.parse_args()

### LOAD DATA
model_name = args.output_name
batch_size = args.batch_size

# Load using mmap_mode='r' to avoid filling up Hoffman2 RAM
mmap_HS_full = np.load(args.hard_npy, mmap_mode='r')
mmap_SS_full = np.load(args.soft_npy, mmap_mode='r')
mmap_neutral_full = np.load(args.neutral_npy, mmap_mode='r')

# Slice based on train_prop
def get_subset(mmap_array, prop):
    num_samples = int(mmap_array.shape[0] * prop)
    return mmap_array[:num_samples, :, :, :]

mmap_HS_used = get_subset(mmap_HS_full, args.train_prop)
mmap_SS_used = get_subset(mmap_SS_full, args.train_prop)
mmap_neutral_used = get_subset(mmap_neutral_full, args.train_prop)

# Split into train and test sets
def train_test_split(mmap_array, test_prop):
    """Split array into train (front) and test (back) portions."""
    if test_prop <= 0.0:
        return mmap_array, None
    n = mmap_array.shape[0]
    split_idx = int(n * (1 - test_prop))
    return mmap_array[:split_idx, :, :, :], mmap_array[split_idx:, :, :, :]

mmap_HS, test_HS = train_test_split(mmap_HS_used, args.test_prop)
mmap_SS, test_SS = train_test_split(mmap_SS_used, args.test_prop)
mmap_neutral, test_neutral = train_test_split(mmap_neutral_used, args.test_prop)

print(f"Using {args.train_prop*100:.0f}% of data, holding out {args.test_prop*100:.0f}% for testing:")
print(f"  Hard    — train: {mmap_HS.shape[0]}, test: {test_HS.shape[0] if test_HS is not None else 0}")
print(f"  Soft    — train: {mmap_SS.shape[0]}, test: {test_SS.shape[0] if test_SS is not None else 0}")
print(f"  Neutral — train: {mmap_neutral.shape[0]}, test: {test_neutral.shape[0] if test_neutral is not None else 0}")

### --------- Parameters -------------------#
val_split = 0.1  # validation carved from the training portion during training

### --------- Build and train model -------------------#

print("BUILDING MODEL")
model = create_model(mmap_neutral)
print(model.summary())

print("TRAINING MODEL")
model, score = train_model(
    model,
    mmap_neutral,
    mmap_HS,
    mmap_SS,
    val_split=val_split,
    batch_size=batch_size,
    path=model_name
)

### --------- Evaluate on held-out test set -------------------#
if args.test_prop > 0.0 and test_HS is not None:
    print("\nEVALUATING ON HELD-OUT TEST SET")
    import matplotlib
    matplotlib.use('Agg')  # non-interactive backend for cluster
    import matplotlib.pyplot as plt
    from sklearn.metrics import precision_recall_curve, auc, confusion_matrix
    from tensorflow.keras.utils import to_categorical

    # Load test data into memory from mmap slices
    X_test = np.concatenate([
        np.array(test_neutral),
        np.array(test_HS),
        np.array(test_SS)
    ], axis=0)

    # One-hot labels: 0=neutral, 1=hard, 2=soft
    y_test_onehot = np.concatenate([
        np.tile([1., 0., 0.], (test_neutral.shape[0], 1)),
        np.tile([0., 1., 0.], (test_HS.shape[0], 1)),
        np.tile([0., 0., 1.], (test_SS.shape[0], 1))
    ], axis=0)

    # Shuffle test set
    shuffle_idx = np.random.permutation(X_test.shape[0])
    X_test = X_test[shuffle_idx]
    y_test_onehot = y_test_onehot[shuffle_idx]
    y_test_int = np.argmax(y_test_onehot, axis=1)

    # --- Basic accuracy ---
    test_loss, test_acc = model.evaluate(X_test, y_test_onehot, batch_size=batch_size, verbose=1)
    print(f"Test Loss: {test_loss:.4f}")
    print(f"Test Accuracy: {test_acc:.4f}")

    y_pred_probs = model.predict(X_test, batch_size=batch_size)
    y_pred_classes = np.argmax(y_pred_probs, axis=1)

    # Per-class accuracy
    class_names = ['Neutral', 'Hard Sweep', 'Soft Sweep']
    for i, name in enumerate(class_names):
        mask = y_test_int == i
        if mask.sum() > 0:
            class_acc = (y_pred_classes[mask] == i).mean()
            print(f"  {name}: {class_acc:.4f} ({mask.sum()} samples)")

    # --- Precision-Recall with combined sweep class (from DANN eval script) ---
    n_classes = 4  # neutral, hard, soft, combined sweeps

    # Build combined sweep ground truth and predicted probability
    combine_sweeps_true = np.logical_or(y_test_onehot[:, 1], y_test_onehot[:, 2]).astype(float).reshape(-1, 1)
    y_true_extended = np.hstack([y_test_onehot, combine_sweeps_true])

    combine_sweeps_pred = np.maximum(y_pred_probs[:, 1], y_pred_probs[:, 2]).reshape(-1, 1)
    y_pred_extended = np.hstack([y_pred_probs, combine_sweeps_pred])

    # Compute PR curves
    precision = dict()
    recall = dict()
    pr_auc = dict()
    for i in range(n_classes):
        precision[i], recall[i], _ = precision_recall_curve(y_true_extended[:, i], y_pred_extended[:, i])
        pr_auc[i] = auc(recall[i], precision[i])

    # --- Plot AUPRC ---
    plot_labels = ["Neutral", "Hard Sweep", "Soft Sweep", "Sweeps (combined)"]
    colors = ["#BAB0AC", "#FF9D9A", "#A0CBE8", "#D4A6C8"]

    plt.figure(figsize=(8, 6))
    for i in range(n_classes):
        plt.plot(recall[i], precision[i], color=colors[i], linestyle='-', lw=2,
                 label=f'{plot_labels[i]} AUPRC = {pr_auc[i]:.3f}')
    plt.xlabel('Recall', fontsize=12)
    plt.ylabel('Precision', fontsize=12)
    plt.title('Precision-Recall Curves (Held-Out Test Set)', fontsize=14)
    plt.legend(fontsize=10)
    plt.tight_layout()

    auprc_plot_path = model_name + "_test_auprc.png"
    plt.savefig(auprc_plot_path, dpi=150)
    plt.close()
    print(f"AUPRC plot saved to: {auprc_plot_path}")

    # --- Confusion matrix plot ---
    cm = confusion_matrix(y_test_int, y_pred_classes)
    cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, data, title, fmt in [
        (axes[0], cm, 'Confusion Matrix (Counts)', 'd'),
        (axes[1], cm_norm, 'Confusion Matrix (Normalized)', '.3f')
    ]:
        im = ax.imshow(data, cmap='Blues', vmin=0)
        ax.set_xticks(range(3))
        ax.set_yticks(range(3))
        ax.set_xticklabels(class_names, fontsize=10)
        ax.set_yticklabels(class_names, fontsize=10)
        ax.set_xlabel('Predicted', fontsize=12)
        ax.set_ylabel('True', fontsize=12)
        ax.set_title(title, fontsize=13)
        for row in range(3):
            for col in range(3):
                ax.text(col, row, format(data[row, col], fmt),
                        ha='center', va='center',
                        color='white' if data[row, col] > data.max() / 2 else 'black')
        fig.colorbar(im, ax=ax, fraction=0.046)
    plt.tight_layout()

    cm_plot_path = model_name + "_test_confusion_matrix.png"
    plt.savefig(cm_plot_path, dpi=150)
    plt.close()
    print(f"Confusion matrix saved to: {cm_plot_path}")

    # --- Save text results ---
    results_path = model_name + "_test_results.txt"
    with open(results_path, 'w') as f:
        f.write(f"Test Loss: {test_loss:.4f}\n")
        f.write(f"Test Accuracy: {test_acc:.4f}\n\n")
        f.write("Per-class accuracy:\n")
        for i, name in enumerate(class_names):
            mask = y_test_int == i
            if mask.sum() > 0:
                class_acc = (y_pred_classes[mask] == i).mean()
                f.write(f"  {name}: {class_acc:.4f} ({mask.sum()} samples)\n")
        f.write("\nAUPRC:\n")
        for i, name in enumerate(plot_labels):
            f.write(f"  {name}: {pr_auc[i]:.4f}\n")
    print(f"Test results saved to: {results_path}")