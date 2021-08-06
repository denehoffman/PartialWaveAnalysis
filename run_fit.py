#!/usr/bin/python3

"""
run_amptools.py: This program performs AmpTools fits over mass-binned data created by divide_data.py.
    Run it without arguments for a more detailed description of its inputs.

Author: Nathaniel Dene Hoffman - Carnegie Mellon University - GlueX Collaboration
Creation Date: 19 July 2021
"""

import argparse
import errno
import sys
import os
import numpy as np
import shutil
from pathlib import Path
from multiprocessing import Pool
import subprocess
import time
from tqdm import tqdm

def resample_params(iteration, config_file):
    """Reads in an AmpTools configuration file and generates a new file with resampled initial fit parameters.

    :param iteration: fit iteration to use in the output file name
    :type iteration: int
    :param config_file: path to template configuration file
    :type config_file: str

    :rtype: str
    :return: path to new resampled configuration file
    """
    with open(config_file, 'r') as config:
        config_lines = config.readlines() # read in the lines from the config template file
    output_filename = config_file[:-4] + f"-{iteration}.cfg" # output file has the same stem as template but with iteration info
    with open(output_filename, 'w') as config:
        output_lines = []
        for line in config_lines:
            if "@seed" in line:
                line = line.replace("@seed", str(np.random.randint(1, high=100000)))
            if line.startswith('initialize'): # if a line starts with initialize...
                line_parts = line.split() # split it on spaces and set the 3rd and 4th fields to randoms
                if line_parts[2] == "cartesian":
                    if line_parts[3] == "@uniform":
                        line_parts[3] = str(np.random.uniform(low=-100.0, high=100.0))
                    if line_parts[4] == "@uniform":
                        line_parts[4] = str(np.random.uniform(low=-100.0, high=100.0))
                elif line_parts[2] == "polar":
                    if line_parts[3] == "@uniform":
                        line_parts[3] = str(np.random.uniform(low=0.0, high=100.0))
                    if line_parts[4] == "@uniform":
                        line_parts[4] = str(np.random.uniform(low=0.0, high=2 * np.pi))
                line = " ".join(line_parts)
                line += "\n"
            output_lines.append(line)
        config.writelines(output_lines) # write the randomized lines to the new output config
    return output_filename


def run_fit(bin_number, iteration, seed, reaction, log_dir, bootstrap, configstem):
    """Runs an iteration of an AmpTools fit on a given bin.
    :param bin_n: the bin number of the fit
    :param iteration: the iteration number of the fit
    :param seed: the fit seed
    :param reaction: the reaction name

    :rtype: None
    :return: None
    """
    start_time = time.time()
    log_dir = Path(log_dir).resolve()
    log_file = log_dir / f"bin_{bin_number}_iteration_{iteration}_seed_{seed}_reaction_{reaction}.log"
    err_file = log_dir / f"bin_{bin_number}_iteration_{iteration}_seed_{seed}_reaction_{reaction}.err"
    os.chdir(str(bin_number)) # cd into the bin directory
    Path(f"./{iteration}").mkdir(exist_ok=True) # create a directory for this iteration if it doesn't already exist
    root_files = Path(".").glob("*.root") # get all the ROOT files for this bin
    print("Copying files")
    for root_file in root_files:
        source = root_file.resolve()
        destination = source.parent / str(iteration) / source.name
        shutil.copy(str(source), str(destination)) # copy all the ROOT files into the iteration directory
    elapsed = time.time() - start_time
    start_time = time.time()
    print(f"Copied: {int(elapsed)} seconds")
    np.random.seed(seed) # set the seed
    if bootstrap:
        config_file = Path(f"{configstem}_{bin_number}_bootstrap.cfg") # locate the (bootstrap) config file
    else:
        config_file = Path(f"{configstem}_{bin_number}.cfg") # locate the config file
    iteration_config_file = resample_params(iteration, str(config_file)) # resample initialization and create a new file for this iteration
    source = Path(iteration_config_file).resolve()
    destination = source.parent / str(iteration) / source.name
    source.replace(destination) # move the iteration config file into its corresponding folder
    os.chdir(str(iteration)) # cd into the iteration folder
    # run the fit: use the iteration config file, send output to stdout, send errors to sterr
    print("Fitting")
    process = subprocess.run(['fit', '-c', iteration_config_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    elapsed = time.time() - start_time
    start_time = time.time()
    print(f"Fit: {int(elapsed)} seconds")
    with open(err_file, 'w') as err_writer:
        err_writer.write(process.stderr) # write any errors to an error file
    with open(log_file, 'w') as log_writer:
        log_writer.write(process.stdout) # write output to a log file
    fit_result = process.stdout
    # check if the fit converged, we'll add this to the filename later 
    if "STATUS=CONVERGED" in fit_result:
        status = "CONVERGED"
        convergence = "C"
    elif "STATUS=FAILED" in fit_result:
        status = "FAILED"
        convergence = "F"
    elif "STATUS=CALL LIMIT" in fit_result:
        status = "CALL_LIMIT"
        convergence = "L"
    else:
        status = "UNKNOWN"
        convergence = "U"
    fit_output = Path(reaction + ".fit").resolve() # locate the output .fit file generated by AmpTools
    if bootstrap:
        fit_output_destination = Path(configstem + "::" + fit_output.stem + f"::{status}.bootstrap").resolve()
    else:
        fit_output_destination = Path(configstem + "::" + fit_output.stem + f"::{status}.fit").resolve()
    fit_output.replace(fit_output_destination) # rename it to contain the fit status (convergence)
    root_files_in_iteration = Path(".").glob("*.root")
    for root_file in root_files_in_iteration:
        root_file.unlink() # remove all the ROOT files in the iteration subdirectory (no longer need them)

    print("Gathering results")
    # Get the fit results here
    commands = []
    with open(destination, 'r') as config:
        for line in config.readlines():
            if line.startswith("amplitude"): # find amplitude lines
                command = []
                wave_name = line.split()[1].strip() # get the parameter name (KsKs::PositiveIm::S0+, for example)
                wave_parts = wave_name.split("::")
                if wave_parts[1].endswith("Re"):
                    command.append("intensity")
                    for pol in ["_AMO", "_000", "_045", "_090", "_135"]:
                        command.append(wave_parts[0] + pol + "::" + wave_parts[1] + "::" + wave_parts[2])
                        command.append(wave_parts[0] + pol + "::" + wave_parts[1].replace("Re", "Im") + "::" + wave_parts[2])
                    commands.append(command)
    with open(destination, 'r') as config:
        command = ["realImag"]
        for line in config.readlines():
            if line.startswith("amplitude"):
                wave_name = line.split()[1].strip() # get the parameter name
                wave_parts = wave_name.split("::")
                if wave_parts[1].endswith("Re"): # we constrain the <reflectivity>Re and <reflectivity>Im amplitudes to be equal
                    command.append(wave_parts[0] + "_000" + "::" + wave_parts[1] + "::" + wave_parts[2]) # polarizations are constrained too
        commands.append(command)
    commands.append(["intensityTotal"])
    commands.append(["likelihood"])
    if bootstrap:
        output_file_name = configstem + "::bootstrap.txt"
    else:
        output_file_name = configstem + "::fit_results.txt"
    with open(output_file_name, 'w') as out_file:
        if convergence == 'C' or convergence == 'L':
            outputs = []
            for command in commands:
                print(command)
                process = subprocess.run(['get_fit_results', str(fit_output_destination), *command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                outputs.append(process.stdout.split("#")[1]) # AmpTools makes a silly warning sometimes
            output = "\t".join(outputs)
            out_file.write(f"{bin_number}\t{iteration}\t{convergence}\t{output}\n") # write fit results to output file (in no particular row order)
    elapsed = time.time() - start_time
    print(f"Finished gathering results: {int(elapsed)} seconds")
    os.chdir("../..")


def main():
    run_fit(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), str(sys.argv[4]), str(sys.argv[5]), str(sys.argv[6]) == "True", str(sys.argv[7]))


if __name__ == "__main__":
    main()
