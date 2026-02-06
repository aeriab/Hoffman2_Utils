### --------- load modules -------------------#
import sys
import os
import argparse
import numpy as np

Utils_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '/u/project/ngarud/Garud_lab/DANN/Utils/'))
sys.path.append(Utils_path)

from CNN_multiclass_data_mergeSims_A100 import *
import gc

parser = argparse.ArgumentParser()
parser.add_argument('hard_npy', type=str, help="Path to Hard sweep npy file")
parser.add_argument('soft_npy', type=str, help="Path to Soft sweep npy file")
parser.add_argument('neutral_npy', type=str, help="Path to Neutral sweep npy file")
parser.add_argument('output_name', type=str, help="Name of the output model")
parser.add_argument('--batch_size', type=int, default=32)
# Added parameter for training proportion (e.g., 0.8 uses 80% of the total data)
parser.add_argument('--train_prop', type=float, default=1.0, help="Proportion of total available data to use")
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

mmap_HS = get_subset(mmap_HS_full, args.train_prop)
mmap_SS = get_subset(mmap_SS_full, args.train_prop)
mmap_neutral = get_subset(mmap_neutral_full, args.train_prop)

print(f"Using {args.train_prop*100}% of data:")
print(f"Hard shape: {mmap_HS.shape}")
print(f"Soft shape: {mmap_SS.shape}")
print(f"Neutral shape: {mmap_neutral.shape}")

### --------- Parameters -------------------#
val_split = 0.1 # This will take the last 10% of the sliced data for validation

### --------- Build and train model -------------------#

print("BUILDING MODEL")
# Note: create_model uses the shape of mmap_neutral to define the Input layer
model = create_model(mmap_neutral)

print(model.summary())

# train model
# The util script uses 'val_split' to carve out the validation set from these inputs
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