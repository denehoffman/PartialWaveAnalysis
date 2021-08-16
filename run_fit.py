#!/usr/bin/env python3
"""
run_amptools.py: This program performs AmpTools fits over mass-binned data created by divide_data.py.
    Run it without arguments for a more detailed description of its inputs.

Author: Nathaniel Dene Hoffman - Carnegie Mellon University - GlueX Collaboration
Creation Date: 19 July 2021
"""

import argparse
import errno
import logging
import os
import shutil
import subprocess
import sys
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from tqdm import tqdm


def resample_params(iteration, config_file):
    """Reads in an AmpTools configuration file and generates a new file with resampled initial fit parameters.

    :param iteration: fit iteration to use in the output file name
    :type iteration: int
    :param config_file: path to template configuration file
    :type config_file: str
    :returns: path to new resampled configuration file

    """
    with open(config_file, 'r') as config:
        config_lines = config.readlines()
    output_filename = config_file[:-4] + f"-{iteration}.cfg"
    with open(output_filename, 'w') as config:
        output_lines = []
        for line in config_lines:
            if "@seed" in line:
                bootseed = str(np.random.randint(1, high=100000))
                line = line.replace("@seed", bootseed)
                logging.info(f"Bootstrap data sampling seed: {bootseed}")
            if line.startswith('initialize'):
                line_parts = line.split()
                if line_parts[2] == "cartesian":
                    if line_parts[3] == "@uniform":
                        line_parts[3] = str(
                            np.random.uniform(low=-100.0, high=100.0))
                    if line_parts[4] == "@uniform":
                        line_parts[4] = str(
                            np.random.uniform(low=-100.0, high=100.0))
                    logging.info(
                        f"Initializing {line_parts[1]} to {line_parts[3]} + {line_parts[4]}i"
                    )
                elif line_parts[2] == "polar":
                    if line_parts[3] == "@uniform":
                        line_parts[3] = str(
                            np.random.uniform(low=0.0, high=100.0))
                    if line_parts[4] == "@uniform":
                        line_parts[4] = str(
                            np.random.uniform(low=0.0, high=2 * np.pi))
                    logging.info(
                        f"Initializing {line_parts[1]} to {line_parts[3]} x Exp(i{line_parts[4]})"
                    )
                line = " ".join(line_parts)
                line += "\n"
            output_lines.append(line)
        config.writelines(output_lines)
    return output_filename


def run_fit(bin_number, iteration, seed, reaction, log_dir, bootstrap,
            configstem):
    """Runs an iteration of an AmpTools fit on a given bin.

    :param bin_n: the bin number of the fit
    :param iteration: the iteration number of the fit
    :param seed: the fit seed
    :param reaction: the reaction name
    :param bin_number: 
    :param log_dir: 
    :param bootstrap: 
    :param configstem: 
    :returns: None

    """
    log_dir = Path(log_dir).resolve()
    log_file = log_dir / f"{configstem}_{bin_number}_{iteration}.log"
    logging.basicConfig(filename=log_file,
                        filemode='a',
                        format="%(asctime)s - %(levelname)s:%(message)s",
                        datefmt="%d/%m/%y %H:%M:%S")
    os.chdir(str(bin_number)) # cd into the bin directory
    logger.info("---------- Start of Fit ----------\n\n")
    logging.info(
        f"Starting AmpTools fit for {configstem}\tBin = {bin_number}\tIteration = {iteration}\tSeed = {seed}"
    )
    if bootstrap:
        logging.info("Bootstrapping is enabled")
    Path(f"./{iteration}").mkdir(exist_ok=True)
    root_files = Path(".").glob("*.root")
    logging.info("Copying ROOT files")
    for root_file in root_files:
        source = root_file.resolve()
        destination = source.parent / str(iteration) / source.name
        shutil.copy(str(source), str(destination))
        logging.debug(f"{str(source)} -> {str(destination)}")
    np.random.seed(seed)
    if bootstrap:
        config_file = Path(f"{configstem}_{bin_number}_bootstrap.cfg")
    else:
        config_file = Path(f"{configstem}_{bin_number}.cfg")
    logging.info(f"Using configuration in {str(config_file)}")
    iteration_config_file = resample_params(iteration, str(config_file))
    source = Path(iteration_config_file).resolve()
    destination = source.parent / str(iteration) / source.name
    source.replace(destination)
    logging.debug(
        f"Renaming resampled config: {str(source)} -> {str(destination)}")
    os.chdir(str(iteration))

    logging.info("Running Fit")
    process = subprocess.run(['fit', '-c', iteration_config_file],
                             stdout=subprocess.PIPE,
                             universal_newlines=True)
    logging.info(process.stdout)
    fit_result = process.stdout

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
    logging.info(f"AmpTools has completed the fit with status = {status}")
    fit_output = Path(reaction + ".fit").resolve()
    if bootstrap:
        fit_output_destination = Path(configstem + "::" + fit_output.stem +
                                      f"::{status}.bootstrap").resolve()
    else:
        fit_output_destination = Path(configstem + "::" + fit_output.stem +
                                      f"::{status}.fit").resolve()
    fit_output.replace(fit_output_destination)
    logging.debug(
        f"Renamed fit file: {str(fit_output)} -> {str(fit_output_destination)}")
    root_files_in_iteration = Path(".").glob("*.root")
    for root_file in root_files_in_iteration:
        root_file.unlink()
    logging.info(f"ROOT files have been removed from the iteration directory")
    logging.info("Gathering Results")

    commands = []
    amplitudes = []
    polarizations = []
    with open(destination, 'r') as config:
        for line in config.readlines():
            if line.startswith("amplitude"):
                amplitudes.append(line.split()[1].strip())
            if line.startswith("define polAngle"):
                polarizations.append(line.split()[1].replace("polAngle", ""))
                logging.debug("Found Polarization: " +
                              line.split()[1].replace("polAngle", ""))
    for wave_name in amplitudes:
        command = []
        wave_parts = wave_name.split("::")
        if wave_parts[1].endswith("Re"):
            command.append("intensity")
            for pol in polarizations:
                command.append(wave_parts[0] + pol + "::" + wave_parts[1] +
                               "::" + wave_parts[2])
                command.append(wave_parts[0] + pol + "::" +
                               wave_parts[1].replace("Re", "Im") + "::" +
                               wave_parts[2])
            commands.append(command)
    command = ["realImag"]
    for wave_name in amplitudes:
        wave_parts = wave_name.split("::")
        if wave_parts[1].endswith("Re"):
            command.append(wave_parts[0] + polarizations[0] + "::" +
                           wave_parts[1] + "::" + wave_parts[2])
    commands.append(command)
    amp_pairs = combinations(amplitudes, 2)
    for amp_pair in amp_pairs:
        wave_name1, wave_name2 = amp_pair
        wave_parts1 = wave_name1.split("::")
        wave_parts2 = wave_name2.split("::")
        if wave_parts1[1].endswith("Re") and wave_parts2[1].endswith("Re"):
            if wave_parts1[1] == wave_parts2[1]:

                command = ["phaseDiff"]
                command.append(wave_parts1[0] + polarizations[0] + "::" +
                               wave_parts1[1] + "::" + wave_parts1[2])
                command.append(wave_parts2[0] + polarizations[0] + "::" +
                               wave_parts2[1] + "::" + wave_parts2[2])
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
                logging.debug(
                    f"Running command: get_fit_results {str(fit_output_destination)} {' '.join(command)}"
                )
                process = subprocess.run(
                    ['get_fit_results',
                     str(fit_output_destination), *command],
                    stdout=subprocess.PIPE,
                    universal_newlines=True)
                if "#" in process.stdout:
                    outputs.append(process.stdout.split("#")[1])
                else:
                    logging.error("An error occurred in get_fit_results!")
                    logging.error(process.stdout)
            output = "\t".join(outputs)
            logging.info(f"Writing fit output to {output_file_name}")
            out_file.write(f"{bin_number}\t{iteration}\t{convergence}\t{output}\n") # write fit results to output file (in no particular row order)
    os.chdir("../..")
    logging.info("----------- End of Fit -----------\n\n")


def main():
    """ """
    run_fit(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]),
            str(sys.argv[4]), str(sys.argv[5]),
            str(sys.argv[6]) == "True", str(sys.argv[7]))


if __name__ == "__main__":
    main()
