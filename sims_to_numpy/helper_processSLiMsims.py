#!/usr/bin/env python3
"""
Helper script to process a single SLiM simulation file.
Called by SLiMsims_to_numpy.py manager script.

REVISED ORDER OF OPERATIONS:
1. Parse SLiM output
2. Create full haplotype matrix
3. Downsample individuals (NUM_SAMPS)
4. Filter out invariant sites (relative only to those sampled individuals)
5. Standardize Major/Minor alleles
6. Crop to fixed SNP-window size
7. Sort haplotypes
"""

import numpy as np
import sys


def parse_slim_sims(path_sims):
    """
    Parses a SLiM simulation output file to extract mutation types, positions, and genome sequences.
    """
    mut_type_dict = {}
    mut_pos_dict = {}
    genomes_dict = {}
    section = None
    
    with open(path_sims, 'r') as file:
        for line in file:
            line = line.strip()
            if line == "Mutations:":
                section = "mutations"
                continue
            elif line == "Individuals:":
                section = None
                continue
            elif line == "Genomes:":
                section = "genomes"
                continue

            if section == "mutations":
                muts_line = line.split()
                if len(muts_line) >= 5:
                    mut_id = muts_line[0]
                    mut_type = muts_line[4]  # selection coefficient
                    mut_pos = muts_line[3]
                    mut_type_dict[mut_id] = mut_type
                    mut_pos_dict[mut_id] = mut_pos
            elif section == "genomes":
                genome_line = line.split()
                if len(genome_line) >= 2:
                    sample_id = genome_line[0]
                    genome_sequence = genome_line[2:]
                    genomes_dict[sample_id] = genome_sequence

    sorted_positions = sorted(set(mut_pos_dict.values()), key=lambda x: int(x))
    pos_to_idx_dict = {pos: idx for idx, pos in enumerate(sorted_positions)}
                    
    return mut_type_dict, mut_pos_dict, pos_to_idx_dict, genomes_dict    


def create_haplotype_matrix(mut_type_dict, mut_pos_dict, pos_to_idx_dict, genomes_dict, 
                            missing_positions=None):
    """Creates a haplotype matrix from genome sequences and mutation data."""
    n_samples = len(genomes_dict)
    n_sites = len(pos_to_idx_dict)
    
    mut_array = np.zeros((n_samples, n_sites, 2), dtype=np.int8)
    mut_array[:, :, 0] = -1  # Default to major allele
    
    missing_set = set(missing_positions) if missing_positions is not None else set()
    
    for pos_idx in missing_set:
        mut_array[:, pos_idx, 0] = 0
        mut_array[:, pos_idx, 1] = 0
    
    for i, sample_id in enumerate(genomes_dict):
        mutation_ids = genomes_dict[sample_id]
        for mut_id in mutation_ids:
            pos_idx = pos_to_idx_dict[mut_pos_dict[mut_id]]
            if pos_idx in missing_set:
                continue
            
            mut_array[i, pos_idx, 0] = 1
            selection_coef = float(mut_type_dict[mut_id])
            if selection_coef != 0:
                mut_array[i, pos_idx, 1] = 1   # Non-synonymous
            else:
                mut_array[i, pos_idx, 1] = -1  # Synonymous

    return mut_array


def sample_haplotypes(mut_array, n_samples):
    """Randomly sample haplotypes from the matrix."""
    total_samples = mut_array.shape[0]
    if n_samples > total_samples:
        raise ValueError(f"Requested {n_samples} samples but only {total_samples} available")
    
    random_indices = np.random.choice(total_samples, n_samples, replace=False)
    random_indices.sort()
    return mut_array[random_indices, :, :]


def crop(mut_array, window_size, sparse=False):
    """
    Revised crop logic: Grabs the first window_size variant sites. 
    This avoids throwing away images where the 'middle' is thin.
    """
    current_size = mut_array.shape[1]
    
    if current_size < window_size:
        raise ValueError(f"Window size {window_size} > available SNPs {current_size}")
    elif current_size == window_size:
        return mut_array
    else:
        if sparse:
            random_sites = np.random.choice(current_size, window_size, replace=False)
            random_sites.sort()
            return mut_array[:, random_sites, :]
        else:
            # Replaced center-crop with leading-edge crop to ensure we get a full window
            return mut_array[:, :window_size, :]


def major_minor(mut_array):
    """Ensure minor allele is 1 based on current sample frequencies."""
    ch1 = mut_array[:, :, 0].copy()
    for site_idx in range(ch1.shape[1]):
        site_data = ch1[:, site_idx]
        non_missing_mask = site_data != 0
        non_missing = site_data[non_missing_mask]
        
        if len(non_missing) > 0:
            minor_freq = np.sum(non_missing == 1) / len(non_missing)
            if minor_freq > 0.5:
                mut_array[non_missing_mask, site_idx, 0] *= -1
    return mut_array


def clusterHaps(numSamples, samples_dict):
    """Cluster identical haplotypes together."""
    haps = {}
    for j in range(numSamples):
        hap_str = ','.join(map(str, samples_dict[j]))
        haps.setdefault(hap_str, [])
        haps[hap_str].append(j)
    
    haps_clumped = {}
    haps_clumped_count = {}
    compared = {}
    
    for key1 in haps.keys():
        if key1 not in compared:
            compared[key1] = 1
            haps_clumped[key1] = list(haps[key1])
            haps_clumped_count[key1] = 1
            
            for key2 in haps.keys():
                if key2 not in compared and haps[key2][0] not in haps_clumped[key1]:
                    distance, s1 = hamming_distance_clump(key1, key2, 0.75)
                    
                    if distance == 0 and key1 != s1:
                        haps_clumped_count[s1] = haps_clumped_count[key1]
                        haps_clumped[s1] = haps_clumped[key1]
                        del haps_clumped_count[key1]
                        del haps_clumped[key1]
                        key1 = s1
                    
                    if distance <= 0:
                        haps_clumped[key1].extend(haps[key2])
                        haps_clumped_count[key1] += 1
                        compared[key2] = 1
    
    return haps_clumped, haps_clumped_count


def hamming_distance_clump(s1, s2, missing_thresh):
    """Hamming distance with missing data support."""
    list_s1 = s1.split(',')
    list_s2 = s2.split(',')
    numMissing_s1 = list_s1.count('0')
    numMissing_s2 = list_s2.count('0')
    
    if (numMissing_s1 >= int(len(list_s1) * missing_thresh) or 
        numMissing_s2 >= int(len(list_s2) * missing_thresh)):
        return len(list_s1), s1
    
    distance = 0
    list_s1_modified = list(list_s1)
    for i in range(len(list_s1)):
        if list_s1_modified[i] != list_s2[i]:
            if list_s1_modified[i] == '0':
                list_s1_modified[i] = list_s2[i]
            elif list_s2[i] == '0':
                pass
            else:
                distance += 1
                if distance > 0:
                    return distance, s1
    
    new_s1 = ','.join(list_s1_modified)
    return distance, new_s1


def sort_haplotypes(mut_array, ordering=None):
    """Sort haplotypes by frequency."""
    if ordering is None or ordering == "none":
        return mut_array
    
    numSamples = mut_array.shape[0]
    samples_dict = {j: row[:, 0] for j, row in enumerate(mut_array)}
    samples_dict_full = {j: row for j, row in enumerate(mut_array)}
    
    haps_clumped, haps_clumped_count = clusterHaps(numSamples, samples_dict)

    if ordering == 'rows_freq':
        haps_sorted = dict(sorted(haps_clumped_count.items(), 
                                   key=lambda x: x[1], reverse=True))
        sorted_data = []
        for hap in haps_sorted.keys():
            for idx in haps_clumped[hap]:
                sorted_data.append(samples_dict_full[idx])
        mut_array = np.array(sorted_data, dtype=np.int8)

    return mut_array


# -----------------------------------------------------------------------------
# Main Execution Block (Order Reordered per Mentor Suggestion)
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python helper_processSLiMsims.py <input> <output> <n_samps> <window> <index> <channels> <sort>")
        sys.exit(1)
    
    PATH_SIMS   = sys.argv[1]
    OUTPUT_PATH = sys.argv[2]
    NUM_SAMPS   = int(sys.argv[3])
    WINDOW_SIZE = int(sys.argv[4])
    INDEX       = int(sys.argv[5])
    CHANNELS    = int(sys.argv[6])
    SORT_ORDER  = sys.argv[7]

    # 1. Parse SLiM output
    mut_type_dict, mut_pos_dict, pos_to_idx_dict, genomes_dict = parse_slim_sims(PATH_SIMS)
    
    # 2. Create initial haplotype matrix (all 30,000 samples)
    mut_array = create_haplotype_matrix(
        mut_type_dict, mut_pos_dict, pos_to_idx_dict, genomes_dict,
        missing_positions=None
    )
    
    # 3. Randomly downsample individuals (e.g., 30,000 -> 154)
    mut_array = sample_haplotypes(mut_array, n_samples=NUM_SAMPS)

    # 4. Filter out invariant sites (looking ONLY at the 154 samples)
    # This removes sites that became fixed/invariant due to the downsampling
    is_variant = np.var(mut_array[:, :, 0], axis=0) > 0
    mut_array = mut_array[:, is_variant, :]
    
    # 5. Correct Major/Minor alleles based on the new 154-sample frequencies
    mut_array = major_minor(mut_array)

    # 6. Crop to window size (using the first WINDOW_SIZE SNPs)
    # This will now succeed as long as the 154 samples have > WINDOW_SIZE variant sites total
    mut_array = crop(mut_array, window_size=WINDOW_SIZE, sparse=False)
    
    # 7. Sort haplotypes
    if SORT_ORDER != "none":
        mut_array = sort_haplotypes(mut_array, SORT_ORDER)

    # 8. Write to memmap
    big_array = np.lib.format.open_memmap(OUTPUT_PATH, dtype=np.int8, mode="r+")
    if CHANNELS == 2:
        big_array[INDEX] = mut_array
    else:
        big_array[INDEX] = mut_array[:, :, 0]

    del big_array