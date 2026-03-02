## How to run Color CNN inference on colored HMP haplotype images

### Prerequisites
- A trained Color CNN model (a `_model.json` and `.weights.h5` file pair). See the [training README](../train_color_CNN/train_color_CNN_README.md) for how to train one.
- A colored `.npy` file with shape `(num_images, num_samples, window_h, 2)` using the dual-channel encoding scheme. See the [HMP csv to numpy README](../hmp_to_color_npy/HMP_csv_to_numpy_README.md) for how to generate these.

### Encoding Scheme
The inference script expects data encoded with two channels:

**Channel 0 (Genotype):**
| Value | Meaning       |
|-------|---------------|
| -1    | Major allele  |
|  0    | Missing data  |
|  1    | Minor allele  |

**Channel 1 (Mutation Type):**
| Value | Meaning                |
|-------|------------------------|
| -1    | Synonymous mutation    |
|  0    | Major allele / missing |
|  1    | Non-synonymous mutation|

If you pass in old single-channel data, the script will exit with an error telling you to use colored `.npy` files.

---


###  Run inference
``` bash
python /u/project/ngarud/Garud_lab/Brendan/Utils/inference_color_CNN/color_CNN_inference.py downsampled_color_CNN/downsampled_color_CNN_model.json downsampled_color_CNN/downsampled_color_CNN.17.weights.h5 r_bromii_120_rows_freq.npy --output r_bromii_predictions.txt
```





| Argument       | Required | Default          | Description                                                                 |
|----------------|----------|------------------|-----------------------------------------------------------------------------|
| `model_name`   | Yes      | —                | Base name of trained model (expects `<name>_model.json` and `<name>.weights.h5`) |
| `image_npy`    | Yes      | —                | Path to colored `.npy` file with shape `(num_images, height, width, 2)`     |
| `--output`     | No       | `predictions.txt`| Output filename                                                             |
| `--batch_size` | No       | 32               | Prediction batch size                                                       |

### Examples
```bash
# Basic usage
python color_CNN_inference.py CNN_color_multiclass_sims_trained GHIST_colored.npy

# Custom output file and batch size
python color_CNN_inference.py \
    CNN_color_multiclass_sims_trained \
    GHIST_colored.npy \
    --output GHIST_results.txt \
    --batch_size 64

# Using a specific epoch checkpoint instead of final weights
python color_CNN_inference.py \
    my_color_model.15 \
    species_data_colored.npy \
    --output species_epoch15_predictions.txt
```

Note: when using a per-epoch checkpoint, the model name should match the checkpoint prefix so that `<name>.weights.h5` resolves correctly. Check that both `<name>_model.json` and `<name>.weights.h5` exist.

### 3. Output format
The output is a CSV file with one row per image:

| Column           | Description                                       |
|------------------|---------------------------------------------------|
| `Image_Index`    | Zero-based index of the image in the `.npy` file  |
| `Predicted_Label`| Winning class: `Neutral`, `Hard_sweep`, or `Soft_sweep` |
| `Confidence`     | Probability of the winning class                  |
| `P_Neutral`      | Predicted probability for Neutral                 |
| `P_Hard`         | Predicted probability for Hard sweep              |
| `P_Soft`         | Predicted probability for Soft sweep              |

Example output:
```
Image_Index,Predicted_Label,Confidence,P_Neutral,P_Hard,P_Soft
0,Neutral,0.981234,0.981234,0.012345,0.006421
1,Hard_sweep,0.874561,0.063210,0.874561,0.062229
2,Soft_sweep,0.712345,0.143210,0.144445,0.712345
```

The script also prints a summary to stdout showing the count and percentage of images classified into each category.

### 4. Tips
- **Batch size:** 32 is fine for most cases. If you run into memory issues with very large `.npy` files, reduce it to 16.
- **Choosing the right model weights:** The default training pipeline saves per-epoch checkpoints. If your test evaluation or validation curves suggest a mid-training epoch performed best, use that checkpoint rather than the final weights.
- **Downstream analysis:** Since all three class probabilities are included in the output, you can threshold or filter predictions by confidence in downstream scripts rather than relying solely on the argmax label.

### Internal details
The inference script imports the model loading function from:
```
/u/project/ngarud/Garud_lab/DANN/Utils/CNN_multiclass_data_mergeSims_A100.py
```