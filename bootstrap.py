#!/usr/bin/python3
import numpy as np
from pathlib import Path
import shutil
import sys
import pandas as pd


input_folder = Path(sys.argv[1]).resolve()
config_name = Path(sys.argv[2]).stem
bin_df = pd.read_csv(input_folder / 'bin_info.txt', delimiter='\t')
bin_numbers = bin_df['bin']

df = pd.read_csv(input_folder / f'{config_name}::fit_results.txt', delimiter='\t', index_col=False)
df.sort_values(['Bin', 'likelihood', 'total_intensity_err'], ascending=[True, False, True], inplace=True)

def mask_first(x):
    result = np.zeros_like(x)
    result[0] = 1
    return result

mask = df.groupby(['Bin'])['Bin'].transform(mask_first).astype(bool)
df_filtered = df.loc[mask]

best_result_by_bin = [df_filtered[df_filtered['Bin'] == n].iloc[0] for n in bin_numbers]


bin_folders = [input_folder / f"{n}" for n in bin_numbers]
for n in bin_numbers:
    config_file = list(bin_folders[n].glob(f"{config_name}_*.cfg"))[0]
    config_dest = config_file.parent / (config_file.stem + "_bootstrap.cfg")
    with open(config_file, 'r') as config_old:
        config_old_lines = config_old.readlines()
    with open(config_dest, 'w') as config_bootstrap:
        for line in config_old_lines:
            if "ROOTDataReader" in line:
                line = line.replace("ROOTDataReader", "ROOTDataReaderBootstrap")
                line = line.replace("\n", " @seed\n")
            if line.startswith("initialize"):
                line = line.replace("polar", "cartesian") # use cartesian coordinates for bootstrap
                wave_name = line.split(" ")[1].split("::")[2]
                fields = line.split(" ")
                fields[3] = str(best_result_by_bin[n][wave_name + "_re"])
                fields[4] = str(best_result_by_bin[n][wave_name + "_im"])
                line = " ".join(fields)
                line += "\n"
            config_bootstrap.write(line)
            
