
# Sort the entire image by row frequency:

python /u/project/ngarud/Garud_lab/Brendan/Utils/one_off_use/numpy_postprocessing/sorting/row_frequency_sort.py hard_unsorted.npy hard_sorted_fixed.npy



# Sort according to row frequency of the middle X SNPs:

python /u/project/ngarud/Garud_lab/Brendan/Utils/one_off_use/numpy_postprocessing/sorting/middle_sort_row_frequency.py unsorted.npy  middle_sorted.npy --window 200
