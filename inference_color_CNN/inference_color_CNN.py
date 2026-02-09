import sys
import os
import argparse
import numpy as np
import tensorflow as tf

# --- Set up path to load custom functions ---
Utils_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '/u/project/ngarud/Garud_lab/DANN/Utils/'))
sys.path.append(Utils_path)
from CNN_multiclass_data_mergeSims_A100 import load_cnn_model

parser = argparse.ArgumentParser(description="Run color CNN inference on colored HMP haplotype images")
parser.add_argument('model_name', type=str, help="Base name of trained model (expects <name>_model.json and <name>.weights.h5)")
parser.add_argument('image_npy', type=str, help="Path to colored .npy file with shape (num_images, num_samples, window_h, 2)")
parser.add_argument('--output', type=str, default='predictions.txt', help="Output filename (default: predictions.txt)")
parser.add_argument('--batch_size', type=int, default=32, help="Prediction batch size (default: 32)")
args = parser.parse_args()

# --- Load model ---
print(f"Loading model: {args.model_name}...")
model = load_cnn_model(args.model_name)

# --- Load colored HMP data ---
print(f"Loading image data: {args.image_npy}...")
data = np.load(args.image_npy, mmap_mode='r')
print(f"Data shape: {data.shape}")

# Validate that the data has the expected 2-channel color encoding
if data.ndim != 4 or data.shape[-1] != 2:
    print(f"ERROR: Expected shape (num_images, height, width, 2) but got {data.shape}")
    print("Make sure you are using colored .npy files (2-channel encoding), not single-channel data.")
    sys.exit(1)

print(f"Loaded {data.shape[0]} images, each {data.shape[1]}x{data.shape[2]} with {data.shape[3]} channels")

# --- Run predictions ---
print(f"Running predictions on {data.shape[0]} images...")
prediction_probabilities = model.predict(data, batch_size=args.batch_size, verbose=1)

# --- Interpret and save results ---
class_labels = ["Neutral", "Hard_sweep", "Soft_sweep"]

predicted_indices = np.argmax(prediction_probabilities, axis=1)
confidences = np.max(prediction_probabilities, axis=1)

print(f"Saving {len(predicted_indices)} predictions to {args.output}...")

with open(args.output, 'w') as f:
    f.write("Image_Index,Predicted_Label,Confidence,P_Neutral,P_Hard,P_Soft\n")
    for i in range(len(predicted_indices)):
        f.write(f"{i},{class_labels[predicted_indices[i]]},{confidences[i]:.6f},"
                f"{prediction_probabilities[i, 0]:.6f},"
                f"{prediction_probabilities[i, 1]:.6f},"
                f"{prediction_probabilities[i, 2]:.6f}\n")

# --- Print summary ---
print("\n--- Prediction Summary ---")
for i, label in enumerate(class_labels):
    count = np.sum(predicted_indices == i)
    pct = count / len(predicted_indices) * 100
    print(f"  {label}: {count} ({pct:.1f}%)")
print(f"\nResults saved to {args.output}")