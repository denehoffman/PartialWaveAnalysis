#!/usr/bin/env python3
"""
bootstrap.py: This program sets up bootstrapping config files using results from an
existing fit.

This script is meant to run internally from run_amptools.py, and it takes two arguments
in order: the path to the input folder and the path to the config file.

Author: Nathaniel Dene Hoffman - Carnegie Mellon University - GlueX Collaboration
Creation Date: 2 August 2021
"""

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

input_folder = Path(sys.argv[1]).resolve()
config_name = Path(sys.argv[2]).stem
bin_df = pd.read_csv(input_folder / "bin_info.txt", delimiter="\t")
bin_numbers = bin_df["bin"]

df = pd.read_csv(input_folder / f"{config_name}::fit_results.txt",
                 delimiter="\t",
                 index_col=False)
df.sort_values(
    ["Bin", "likelihood", "total_intensity_err"],
    ascending=[True, False, True],
    inplace=True,
) # sort the results by bin number, likelihood, and total intensity error


def mask_first(x):
    """mask_first is a helper function that returns a 0-list with a 1 in the first
    position.
    
    :param x: The list from which to create a mask list

    """
    result = np.zeros_like(x)
    result[0] = 1
    return result

# this mask selects the best fit in each bin by the sorting done before
mask = df.groupby(["Bin"])["Bin"].transform(mask_first).astype(bool)
df_filtered = df.loc[mask]

best_result_by_bin = [
    df_filtered[df_filtered["Bin"] == n].iloc[0] for n in bin_numbers
] # a list with the best iteration in each bin

bin_folders = [input_folder / f"{n}" for n in bin_numbers]
for n in bin_numbers: # for each bin
    config_file = list(bin_folders[n].glob(f"{config_name}_*.cfg"))[0]
    config_dest = config_file.parent / (config_file.stem + "_bootstrap.cfg")
    with open(config_file, "r") as config_old:
        config_old_lines = config_old.readlines() # read in the original config
    with open(config_dest, "w") as config_bootstrap:
        for line in config_old_lines:
            if "ROOTDataReader" in line: # replace with bootstrap data reader
                line = line.replace("ROOTDataReader", "ROOTDataReaderBootstrap")
                line = line.replace("\n", " @seed\n")
            if line.startswith("initialize"): # set initialization to best value
                line = line.replace("polar", "cartesian")
                wave_name = line.split(" ")[1].split("::")[2]
                fields = line.split(" ")
                fields[3] = str(best_result_by_bin[n][wave_name + "_re"])
                fields[4] = str(best_result_by_bin[n][wave_name + "_im"])
                line = " ".join(fields)
                line += "\n"
            config_bootstrap.write(line)
