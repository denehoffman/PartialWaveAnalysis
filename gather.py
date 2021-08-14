#!/usr/bin/env python3
import argparse
import errno
import os
import subprocess
import sys
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from colorama import Fore
from tqdm import tqdm

"""
gather_fits.py: This program collects data from AmpTools fits and stores them in a single location.
    Run it without arguments for a more detailed description of its inputs.
    Author: Nathaniel Dene Hoffman - Carnegie Mellon University - GlueX Collaboration
    Creation Date: 13 July 2021
"""


def gather(output_dir, config_file, bootstrap, n_iterations):
    """Gathers cumulative results of fits, calculates intensities with get_fit_results, and outputs to a tab-separated file
    
    This method goes through all the bin directories, finds converged fits, and
    collects fit_results.txt files produced by run_fit.py

    :param output_dir: directory which contains bin subdirectories
    :type output_dir: Path
    :param config_file: template of AmpTools configuration file
    :type config_file: Path
    :param bootstrap: collect bootstrap files
    :type bootstrap: bool
    :param n_iterations: number of iterations to collect
    :type n_iterations: int
    :returns: None

    """
    print("Gathering Results")
    headers = []
    amplitudes = []
    with open(config_file, 'r') as config:
        for line in config.readlines():
            if line.startswith("amplitude"):
                amplitudes.append(line.split()[1].strip())
    for wave_name in amplitudes:
        wave_parts = wave_name.split("::")
        if wave_parts[1].endswith("Re"):
            headers.append(wave_parts[2] + "_INT")
            headers.append(wave_parts[2] + "_err_INT")
            headers.append(wave_parts[2] + "_AC_INT")
            headers.append(wave_parts[2] + "_err_AC_INT")
    for wave_name in amplitudes:
        wave_parts = wave_name.split("::")
        if wave_parts[1].endswith("Re"):
            headers.append(wave_parts[2] + "_re")
            headers.append(wave_parts[2] + "_im")
    amp_pairs = combinations(amplitudes, 2)
    for amp_pair in amp_pairs:
        wave_name1, wave_name2 = amp_pair
        wave_parts1 = wave_name1.split("::")
        wave_parts2 = wave_name2.split("::")
        if wave_parts1[1].endswith("Re") and wave_parts2[1].endswith("Re"):
            if wave_parts1[1] == wave_parts2[1]:
                print(wave_name1, wave_name2)

                headers.append(wave_parts1[2] + "_" + wave_parts2[2] + "_PHASE")
                headers.append(wave_parts1[2] + "_" + wave_parts2[2] +
                               "_err_PHASE")
    headers.append("total_intensity")
    headers.append("total_intensity_err")
    headers.append("total_intensity_AC")
    headers.append("total_intensity_err_AC")
    headers.append("likelihood")
    if bootstrap:
        output_file_name = config_file.stem + "::bootstrap.txt"
    else:
        output_file_name = config_file.stem + "::fit_results.txt"
    with open(output_dir / output_file_name, 'w') as out_file:
        header = "\t".join(headers)
        out_file.write(f"Bin\tIteration\tConvergence\t{header}\n")
        bin_dirs = [
            bindir for bindir in output_dir.glob("*") if bindir.is_dir()
        ]
        bin_converged_total = np.zeros_like(bin_dirs)
        bin_total_iterations = np.zeros_like(bin_dirs)
        for bin_dir in tqdm(bin_dirs):
            bin_num_string = bin_dir.name
            iter_dirs = [
                iterdir for iterdir in bin_dir.glob("*")
                if iterdir.is_dir() and int(iterdir.name) < n_iterations
            ]
            for iteration_dir in iter_dirs:
                iteration_num_string = iteration_dir.name
                if bootstrap:
                    fit_results = iteration_dir / (config_file.stem +
                                                   "::bootstrap.txt")
                else:
                    fit_results = iteration_dir / (config_file.stem +
                                                   "::fit_results.txt")
                bin_total_iterations[int(bin_num_string)] += 1
                converged = "U"
                if fit_results.exists():
                    with open(fit_results, "r") as fit_results_reader:
                        line = fit_results_reader.readline()
                        if line != "":
                            converged = line.split("\t")[2]
                    if converged == "C":
                        bin_converged_total[int(bin_num_string)] += 1
                    with open(fit_results) as fit_reader:
                        out_file.write(fit_reader.read())
        print("Convergence Results:")
        for i, bin_converged_num in enumerate(bin_converged_total):
            percent_converged = bin_converged_total[i] / bin_total_iterations[i]
            if percent_converged == 0:
                color = Fore.RED
            elif percent_converged <= 0.25:
                color = Fore.YELLOW
            elif percent_converged <= 0.80:
                color = Fore.BLUE
            else:
                color = Fore.GREEN
            print(
                f"{color}Bin {i}: {bin_converged_total[i]}/{bin_total_iterations[i]}\t{Fore.WHITE}",
                end='')
        print()


"""
Script begins here:
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Runs AmpTools fits on each mass bin")
    parser.add_argument("-d",
                        "--directory",
                        required=True,
                        help="the input directory (output of divide_data.py)")
    parser.add_argument("-c",
                        "--config",
                        required=True,
                        help="path to the AmpTools config template file")
    parser.add_argument(
        "--bootstrap",
        action='store_true',
        help=
        "use bootstrapping (must have run a fit already and run bootstrap.py)")
    parser.add_argument("-n",
                        "--niterations",
                        type=int,
                        required=True,
                        help="number of iterations")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    bin_directory = Path(args.directory).resolve()
    if bin_directory.is_dir():
        print(f"Input Directory: {bin_directory}")
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                args.directory)

    config_template = Path(args.config).resolve()
    if config_template.is_file():
        print(f"Config Template: {config_template}")
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                args.config)

    gather(bin_directory, config_template, args.bootstrap, args.niterations)
