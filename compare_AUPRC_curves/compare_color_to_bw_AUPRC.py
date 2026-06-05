import sys
import os
import argparse
import numpy as np
import tensorflow as tf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, auc
from tensorflow.keras.models import model_from_json
from tensorflow.keras.optimizers import Adam

# --- Path Setup ---
Utils_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '/u/project/ngarud/Garud_lab/Brendan/Utils/DANN_Utils/'))
sys.path.append(Utils_path)

# Import custom functions from your utility script
from CNN_zero_as_missing_A100 import (
    GradReverse, custom_bce, custom_categorical_ce, 
    custom_categorical_accuracy
)

def get_args():
    parser = argparse.ArgumentParser(description="Compare Color and B&W CNN AUPRC with explicit file paths.")
    
    # Color Model Explicit Files
    parser.add_argument('--color_json', type=str, required=True, help="Path to Color model JSON")
    parser.add_argument('--color_weights', type=str, required=True, help="Path to Color model weights (.h5)")
    
    # B&W Model Explicit Files
    parser.add_argument('--bw_json', type=str, required=True, help="Path to B&W model JSON")
    parser.add_argument('--bw_weights', type=str, required=True, help="Path to B&W model weights (.h5)")
    
    # Data Paths
    parser.add_argument('--hard_npy', type=str, required=True)
    parser.add_argument('--soft_npy', type=str, required=True)
    parser.add_argument('--neutral_npy', type=str, required=True)
    
    # Parameters
    parser.add_argument('--train_prop', type=float, default=1.0)
    parser.add_argument('--test_prop', type=float, default=0.2)
    parser.add_argument('--batch_size', type=int, default=32)
    
    # Zoom Parameter (Optional minimum limit for top-right focus)
    parser.add_argument('--zoom', type=float, default=None, 
                        help="Lower limit for Precision and Recall axes to zoom into top-right (e.g., 0.7). Leaves full range if omitted.")
    
    return parser.parse_args()

def load_model_explicit(json_path, weights_path):
    """Manually loads and compiles model using explicit file paths."""
    tf.keras.utils.get_custom_objects()['GradReverse'] = GradReverse
    tf.keras.utils.get_custom_objects()['custom_bce'] = custom_bce
    tf.keras.utils.get_custom_objects()['custom_categorical_ce'] = custom_categorical_ce
    tf.keras.utils.get_custom_objects()['custom_categorical_accuracy'] = custom_categorical_accuracy

    if not os.path.exists(json_path) or not os.path.exists(weights_path):
        print(f"ERROR: Missing model files: {json_path} or {weights_path}")
        sys.exit(1)

    with open(json_path, "r") as f:
        model = model_from_json(f.read())
    model.load_weights(weights_path)

    optimizer = Adam(learning_rate=1e-5)
    model.compile(optimizer=optimizer,
                  loss={'classifier': custom_categorical_ce},
                  metrics={'classifier': custom_categorical_accuracy})
    return model

def prepare_test_data(args):
    """Replicates the subsetting and splitting logic for the test set."""
    def get_subset(mmap_array, prop):
        num_samples = int(mmap_array.shape[0] * prop)
        return mmap_array[:num_samples, :, :, :]

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

    # Load Models Explicitly
    print(f"\n--- Loading Color Model ---\nJSON: {args.color_json}\nWeights: {args.color_weights}")
    color_model = load_model_explicit(args.color_json, args.color_weights)
    
    print(f"\n--- Loading B&W Model ---\nJSON: {args.bw_json}\nWeights: {args.bw_weights}")
    bw_model = load_model_explicit(args.bw_json, args.bw_weights)

    # Run Inference
    print("\nEvaluating Color model...")
    c_neu, c_swp = get_pr_data(color_model, X_test, y_test_onehot, args.batch_size)

    print("Evaluating B&W model...")
    X_test_bw = X_test[:, :, :, 0:1] # Slice for single-channel genotype
    b_neu, b_swp = get_pr_data(bw_model, X_test_bw, y_test_onehot, args.batch_size)

    # --- Plotting ---
    plt.figure(figsize=(9, 7))
    plt.plot(c_neu[0], c_neu[1], color='#76b7b2', linestyle='-', lw=2.5, label=f'Color Neutral (AUC={c_neu[2]:.3f})')
    plt.plot(c_swp[0], c_swp[1], color='#e15759', linestyle='-', lw=2.5, label=f'Color Sweep (AUC={c_swp[2]:.3f})')
    plt.plot(b_neu[0], b_neu[1], color='#76b7b2', linestyle='--', lw=2, label=f'B&W Neutral (AUC={b_neu[2]:.3f})')
    plt.plot(b_swp[0], b_swp[1], color='#e15759', linestyle='--', lw=2, label=f'B&W Sweep (AUC={b_swp[2]:.3f})')

    # Apply Dynamic Zoom if passed via CLI arguments
    if args.zoom is not None:
        if not (0.0 <= args.zoom < 1.0):
            print(f"WARNING: Zoom value {args.zoom} out of valid bounds (0.0 to 1.0). Ignoring zoom.")
        else:
            plt.xlim(args.zoom, 1.01)
            plt.ylim(args.zoom, 1.02)
            print(f"Applying zoom focused on top-right: [{args.zoom}, 1.0]")

    plt.xlabel('Recall', fontsize=12); plt.ylabel('Precision', fontsize=12)
    plt.title('AUPRC Comparison: Color vs Black & White CNN', fontsize=14)
    
    # Moves the legend dynamically based on zoom to avoid covering data
    legend_loc = 'lower left' if args.zoom is None else 'lower center'
    plt.legend(loc=legend_loc, fontsize=10)
    
    plt.grid(alpha=0.3); plt.tight_layout()

    save_path = "CNN_Comparison_AUPRC.png"
    plt.savefig(save_path, dpi=300)
    print(f"\nSuccess! Comparison plot saved to: {save_path}")

if __name__ == "__main__":
    main()