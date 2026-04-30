import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, auc

# --- Path Setup ---
Utils_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '/u/project/ngarud/Garud_lab/Brendan/Utils/DANN_Utils/'))
sys.path.append(Utils_path)

from CNN_zero_as_missing_A100 import *

def get_args():
    parser = argparse.ArgumentParser(description="Compare Color and B&W CNN AUPRC curves.")
    parser.add_argument('--color_model', type=str, required=True, help="Prefix for Color model")
    parser.add_argument('--bw_model', type=str, required=True, help="Prefix for B&W model")
    parser.add_argument('--hard_npy', type=str, required=True)
    parser.add_argument('--soft_npy', type=str, required=True)
    parser.add_argument('--neutral_npy', type=str, required=True)
    parser.add_argument('--train_prop', type=float, default=1.0)
    parser.add_argument('--test_prop', type=float, default=0.2)
    parser.add_argument('--batch_size', type=int, default=32)
    return parser.parse_args()

def prepare_test_data(args):
    """Loads the full test set (all channels)."""
    def get_subset(mmap_array, prop):
        num_samples = int(mmap_array.shape[0] * prop)
        return mmap_array[:num_samples, :, :, :] # Load full channels here

    def get_test_portion(mmap_array, test_prop):
        n = mmap_array.shape[0]
        split_idx = int(n * (1 - test_prop))
        return mmap_array[split_idx:, :, :, :]

    print("Loading and splitting data...")
    m_HS = np.load(args.hard_npy, mmap_mode='r')
    m_SS = np.load(args.soft_npy, mmap_mode='r')
    m_Neu = np.load(args.neutral_npy, mmap_mode='r')

    test_HS = get_test_portion(get_subset(m_HS, args.train_prop), args.test_prop)
    test_SS = get_test_portion(get_subset(m_SS, args.train_prop), args.test_prop)
    test_Neu = get_test_portion(get_subset(m_Neu, args.train_prop), args.test_prop)

    X_test = np.concatenate([np.array(test_Neu), np.array(test_HS), np.array(test_SS)], axis=0)
    
    y_test_onehot = np.concatenate([
        np.tile([1., 0., 0.], (test_Neu.shape[0], 1)),
        np.tile([0., 1., 0.], (test_HS.shape[0], 1)),
        np.tile([0., 0., 1.], (test_SS.shape[0], 1))
    ], axis=0)

    return X_test, y_test_onehot

def get_pr_data(model, X_test, y_test_onehot, batch_size):
    """Computes precision/recall for Neutral and Combined Sweeps."""
    y_probs = model.predict(X_test, batch_size=batch_size)
    
    sweep_true = np.logical_or(y_test_onehot[:, 1], y_test_onehot[:, 2]).astype(float)
    sweep_pred = np.maximum(y_probs[:, 1], y_probs[:, 2])

    p_neu, r_neu, _ = precision_recall_curve(y_test_onehot[:, 0], y_probs[:, 0])
    auc_neu = auc(r_neu, p_neu)

    p_swp, r_swp, _ = precision_recall_curve(sweep_true, sweep_pred)
    auc_swp = auc(r_swp, p_swp)

    return (r_neu, p_neu, auc_neu), (r_swp, p_swp, auc_swp)

def main():
    args = get_args()
    X_test, y_test_onehot = prepare_test_data(args)

    # Load Models
    print(f"Loading Color Model: {args.color_model}")
    color_model = load_cnn_model(args.color_model)
    
    print(f"Loading B&W Model: {args.bw_model}")
    bw_model = load_cnn_model(args.bw_model)

    # Evaluate Color Model (uses full X_test)
    print("Evaluating Color model...")
    c_neu, c_swp = get_pr_data(color_model, X_test, y_test_onehot, args.batch_size)

    # Evaluate B&W Model (uses only genotype channel: [:, :, :, 0:1])
    print("Evaluating B&W model...")
    # This slice is critical to match the B&W training logic
    X_test_bw = X_test[:, :, :, 0:1] 
    b_neu, b_swp = get_pr_data(bw_model, X_test_bw, y_test_onehot, args.batch_size)

    # --- Plotting ---
    plt.figure(figsize=(9, 7))
    
    # Color CNN (Solid)
    plt.plot(c_neu[0], c_neu[1], color='#76b7b2', linestyle='-', lw=2.5, 
             label=f'Color Neutral (AUC={c_neu[2]:.3f})')
    plt.plot(c_swp[0], c_swp[1], color='#e15759', linestyle='-', lw=2.5, 
             label=f'Color Sweep (AUC={c_swp[2]:.3f})')

    # B&W CNN (Dashed)
    plt.plot(b_neu[0], b_neu[1], color='#76b7b2', linestyle='--', lw=2, 
             label=f'B&W Neutral (AUC={b_neu[2]:.3f})')
    plt.plot(b_swp[0], b_swp[1], color='#e15759', linestyle='--', lw=2, 
             label=f'B&W Sweep (AUC={b_swp[2]:.3f})')

    plt.xlabel('Recall', fontsize=12)
    plt.ylabel('Precision', fontsize=12)
    plt.title('AUPRC Comparison: Color vs Black & White CNN', fontsize=14)
    plt.legend(loc='lower left', fontsize=10)
    plt.grid(alpha=0.3)
    plt.tight_layout()

    save_path = "CNN_Comparison_AUPRC.png"
    plt.savefig(save_path, dpi=300)
    print(f"Comparison plot saved to: {save_path}")

if __name__ == "__main__":
    main()