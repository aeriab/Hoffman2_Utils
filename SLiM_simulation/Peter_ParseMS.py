# Peter Laurin
# July 9, 2024
# updated to reflect multi-allelic sites
# adapted from Mariana Harris (ParseMS.py)

from optparse import OptionParser
import sys
import numpy as np
import pandas as pd

def parseMS(File, windowSize, fitnessEffects):
    if not fitnessEffects:
        MS = pd.read_csv(File, skiprows=2,header=None)  # ,sep="\s*",engine='python'
        variant_start = 1
        selection_cos = np.zeros(len(MS.iloc[1,0]))
    if fitnessEffects:

        # NEED TO LOOK HOW MULTI-ALLELIC SITES AFFECT SELECTION COEFFICIENTS!!
        MS = pd.read_csv(File, skiprows=3,header=None)
        variant_start = 2
        selection_cos = list(MS.iloc[1, :])[0].split()[1:]
        selection_cos = -np.abs(np.array(selection_cos, dtype=float)) #convert all to negative, for now
    pos = list(MS.iloc[0, :])[0].split()
    pos = np.round(np.array(pos[1:], dtype= float) * windowSize)
    uniq_pos, counts = np.unique(pos, return_counts=True)
    duplicated_pos = uniq_pos[counts > 1]
    for d_pos in duplicated_pos:
        indices = np.where(np.array(pos) == d_pos)[0]
        for indx in indices[1:]:
            next_pos = d_pos
            while next_pos in pos:
                next_pos += 1
            pos[indx] = next_pos

    MS = (MS.iloc[variant_start:, :])[0].apply(lambda x: pd.Series(list(x)))
    MS.columns = (np.rint(pos)).astype(int)
    MS = MS.transpose().astype(int)
    site_types = np.repeat("nonsyn",len(selection_cos))
    site_types[selection_cos == 0] = "syn"
    MS.insert(0, "site_type", site_types)
    #There are now some sites out of order
    MS.sort_index(inplace=True)
    #MS = MS.loc[MS.index.drop_duplicates(keep=False)]#drop multiple allele sites
    

    return MS
    
def main():
    parser = OptionParser()
    parser.add_option("-i", "--inFile", type="string", default=None)
    parser.add_option("-o", "--outFile", type="string")
    parser.add_option("-w", "--windowSize", type="int", default=29999)
    parser.add_option("--fitnessEffects", action="store_true", default=False)

    options, args = parser.parse_args()
    infile = options.inFile
    outfile = options.outFile
    window_size = options.windowSize
    fitness_effects = options.fitnessEffects

    site_df = parseMS(infile, window_size, fitness_effects)

    site_df.to_csv(outfile, index=True, index_label="site_pos")
    

if __name__ == "__main__":
    main()