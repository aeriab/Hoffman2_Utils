# Minimal — no genomic positions, plot iLDS_pval as raw value
python /u/project/ngarud/Garud_lab/Brendan/Utils/plot_inference_results/compare_results_plot.py predictions.csv ricky_results.csv

# With site_indices and treating iLDS_pval as a p-value (gives -log10 scale)
python /u/project/ngarud/Garud_lab/Brendan/Utils/plot_inference_results/compare_results_plot.py predictions.csv ricky_results.csv \
    --site_indices windows_site_indices.npy \
    --ricky_pval \
    --ricky_col iLDS_pval \
    --ricky_contig FP929051 \
    --title "R. bromii sweep scan — FP929051" \
    --output FP929051_comparison.png

# Try rNrS_pval instead
python /u/project/ngarud/Garud_lab/Brendan/Utils/plot_inference_results/compare_results_plot.py predictions.csv ricky_results.csv \
    --site_indices windows_site_indices.npy \
    --ricky_pval --ricky_col rNrS_pval