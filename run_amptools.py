#!/usr/bin/python3

"""
run_amptools.py: This program performs AmpTools fits over mass-binned data created by divide_data.py.
    Run it without arguments for a more detailed description of its inputs.

Author: Nathaniel Dene Hoffman - Carnegie Mellon University - GlueX Collaboration
Creation Date: 5 July 2021
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
from tqdm import tqdm
from colorama import Fore
import time

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


def run_fit(bin_iterations_seed_reaction_tuple):
    """Runs an iteration of an AmpTools fit on a given bin.
    :param bin_iterations_seed_reaction_tuple: tuple(int, int, int, str), contains the bin number, iteration number, seed, and reaction name

    :rtype: None
    :return: None
    """
    bin_number, iteration, seed, reaction = bin_iterations_seed_reaction_tuple # get info for this run
    log_file = log_dir / f"bin_{bin_number}_iteration_{iteration}_seed_{seed}_reaction_{reaction}.log"
    err_file = log_dir / f"bin_{bin_number}_iteration_{iteration}_seed_{seed}_reaction_{reaction}.err"
    os.chdir(str(bin_number)) # cd into the bin directory
    Path(f"./{iteration}").mkdir(exist_ok=True) # create a directory for this iteration if it doesn't already exist
    root_files = Path(".").glob("*.root") # get all the ROOT files for this bin
    for root_file in root_files:
        source = root_file.resolve()
        destination = source.parent / str(iteration) / source.name
        shutil.copy(str(source), str(destination)) # copy all the ROOT files into the iteration directory

    np.random.seed(seed) # set the seed
    config_file = [f for f in Path(".").glob(f"*_{bin_number}.cfg")][0] # locate the config file
    iteration_config_file = resample_params(iteration, str(config_file)) # resample initialization and create a new file for this iteration
    source = Path(iteration_config_file).resolve()
    destination = source.parent / str(iteration) / source.name
    source.replace(destination) # move the iteration config file into its corresponding folder
    os.chdir(str(iteration)) # cd into the iteration folder
    # run the fit: use the iteration config file, send output to stdout, send errors to sterr
    process = subprocess.run(['fit', '-c', iteration_config_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
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
    fit_output_destination = Path(fit_output.stem + f"::{status}.fit").resolve()
    fit_output.replace(fit_output_destination) # rename it to contain the fit status (convergence)
    root_files_in_iteration = Path(".").glob("*.root")
    for root_file in root_files_in_iteration:
        root_file.unlink() # remove all the ROOT files in the iteration subdirectory (no longer need them)

    # Get the fit results here
    with open(destination, 'r') as config:
        config_lines = config.readlines() # read in the lines from the config template file
        commands = []
        for line in config_lines:
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
        commands.append(["intensityTotal"])
        commands.append(["likelihood"])
    with open("fit_results.txt", 'w') as out_file:
        if convergence == 'C' or convergence == 'L':
            outputs = []
            for command in commands:
                process = subprocess.run(['get_fit_results', str(fit_output_destination), *command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                outputs.append(process.stdout.split("#")[1]) # AmpTools makes a silly warning sometimes
            output = "\t".join(outputs)
        out_file.write(f"{bin_number}\t{iteration}\t{convergence}\t{output}\n") # write fit results to output file (in no particular row order)
    os.chdir("../..")


def gather(output_dir, config_file):
    """Gathers cumulative results of fits, calculates intensities with get_fit_results, and outputs to a tab-separated file

    This method goes through all the bin directories, finds converged fits, and runs a C program that uses IUAmpTools/FitResults.h
    to collect the amplitude fit values and their errors, as well as the total likelihood. A particular bin and iteration
    converged fit represents one line in the final output file, which will be a tab-separated file with unordered rows.
    The first column is the bin number, the second is the iteration number, and the rest are paired columns of parameters
    and parameter errors. The input to the C program will be the fit file followed by a list of parameters.

    :param output_dir: directory which contains bin subdirectories
    :type output_dir: Path
    :param config_file: template of AmpTools configuration file
    :type config_file: Path

    :rtype: None
    :return: None
    """
    print("Gathering Results")
    with open(config_file, 'r') as config:
        config_lines = config.readlines() # read in the lines from the config template file
        headers = []
        for line in config_lines:
            if line.startswith("amplitude"): # find amplitude lines
                wave_name = line.split()[1].strip() # get the parameter name (KsKs::PositiveIm::S0+, for example)
                wave_parts = wave_name.split("::")
                if wave_parts[1].endswith("Re"):
                    headers.append(wave_parts[2])
                    headers.append(wave_parts[2] + "_err")
        headers.append("total_intensity")
        headers.append("total_intensity_err")
        headers.append("likelihood")
    with open(output_dir / "fit_results.txt", 'w') as out_file:
        header = "\t".join(headers)
        out_file.write(f"Bin\tIteration\tConvergence\t{header}\n") # print the header to the output file
        bin_dirs = [bindir for bindir in output_dir.glob("*") if bindir.is_dir()]
        bin_converged_total = np.zeros_like(bin_dirs)
        bin_total_iterations = np.zeros_like(bin_dirs)
        for bin_dir in tqdm(bin_dirs): # for each bin subdirectory
            bin_num_string = bin_dir.name
            iter_dirs = [iterdir for iterdir in bin_dir.glob("*") if iterdir.is_dir()] # for each iteration subdirectory
            for iteration_dir in iter_dirs:
                iteration_num_string = iteration_dir.name
                fit_files = [fit.resolve() for fit in iteration_dir.glob("*.fit")]
                latest_fit_file = max(fit_files, key=os.path.getctime)
                fit_results = iteration_dir / "fit_results.txt"
                bin_total_iterations[int(bin_num_string)] += 1
                if "CONVERGED" in latest_fit_file.name: # only collect converged fits
                    bin_converged_total[int(bin_num_string)] += 1
                with open(fit_results) as fit_reader:
                    out_file.write(fit_reader.read()) # write fit results to output file (in no particular row order)
        print("Convergence Results:")
        for i, bin_converged_num in enumerate(bin_converged_total):
            percent_converged = bin_converged_total[i] / bin_total_iterations[i]
            if percent_converged == 0:
                color = Fore.RED
            elif percent_converged <= 0.25:
                color = Fore.YELLOW
            elif percent_converged <=0.80:
                color = Fore.BLUE
            else:
                color = Fore.GREEN
            print(f"{color}Bin {i}: {bin_converged_total[i]}/{bin_total_iterations[i]}\t{Fore.WHITE}", end='')
        print()


"""
Script begins here:
"""
parser = argparse.ArgumentParser(description="Runs AmpTools fits on each mass bin")
parser.add_argument("-d", "--directory", required=True, help="the input directory (output of divide_data.py)")
parser.add_argument("-p", "--processes", default=60, type=int, help="number of processes to spawn")
parser.add_argument("-s", "--seed", default=1, type=int, help="set the seed (default = 1)")
parser.add_argument("-i", "--iterations", default=1, type=int, help="set the number of fits to perform for each bin")
parser.add_argument("-c", "--config", required=True, help="path to the AmpTools config template file")
parser.add_argument("--slurm", action='store_true', help="use SLURM's sbatch command instead of python3 multiprocessing")
queue_group = parser.add_mutually_exclusive_group()
queue_group.add_argument('-r', '--red', action='store_true', help="run on red queue (default)")
queue_group.add_argument('-g', '--green', action='store_true', help="run on green queue")
if len(sys.argv) == 1: # if the user doesn't supply any arguments, print the help string and exit
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
if args.green:
    cpu_memory = 1590
    threads = 4
    queue = "green"
else:
    cpu_memory = 1990
    threads = 4
    queue = "red"
memory = threads * cpu_memory

bin_directory = Path(args.directory).resolve()
if bin_directory.is_dir(): # check if directory with all the separated bins exists
    print(f"Input Directory: {bin_directory}")
else:
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.directory)

config_template = Path(args.config).resolve()
if config_template.is_file(): # check if config file exists
    print(f"Config Template: {config_template}")
else:
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.config)

reaction = "REACTION"
with open(config_template, 'r') as config_temp:
    lines = config_temp.readlines()
    for line in lines:
        if line.startswith("reaction"):
            reaction = line.split()[1] # get the reaction name from the config file
            # the line reads "reaction KsKs" or whatever you call your reaction
            # this is also what AmpTools calls the .fit file later (like KsKs.fit)
print(f"Reaction Name: {reaction}")

log_dir = Path("logs").resolve()
log_dir.mkdir(exist_ok=True) # create a log directory if it doesn't already exist

np.random.seed(args.seed) # seed the RNG
seeds = [np.random.randint(1, high=100000) for _ in range(args.iterations)] # create seeds for each iteration
n_bins = len([bin_dir for bin_dir in bin_directory.glob("*") if bin_dir.is_dir()])
bin_iterations_seed_reaction_tuple = [(i, j, seeds[j], reaction) for i in range(n_bins) for j in range(args.iterations)]
# tuple contains the bin number, the iteration number, a seed specific to the iteration (iteration "j" across all bins will have the same seed, so in theory
# it should also have the same starting parameters, in case that is important later on), and the reaction name

os.chdir(str(bin_directory)) # cd into the directory containing bin subdirectories
username = os.getlogin()
if args.slurm:
    for tup in tqdm(bin_iterations_seed_reaction_tuple):
        bin_number, iteration, seed, reaction = tup
        bin_n, it_n, _, rxn = tup
        log_out = log_dir / f"{rxn}_{bin_n}_{it_n}_SLURM.out"
        log_err = log_dir / f"{rxn}_{bin_n}_{it_n}_SLURM.err"
        os.system(f"sbatch --job-name={rxn}_{bin_n}_{it_n} --ntasks={threads} --partition={queue} --mem={memory} --time=30:00 --output={str(log_out)} --error={str(log_err)} --quiet ../run_fit_slurm.csh {bin_number} {iteration} {seed} {reaction} {str(log_dir)}")
        time.sleep(1)

    finished_running = False
    while not finished_running:
        n_jobs_queued = len(subprocess.run(['squeue', '-h', '-u', username], stdout=subprocess.PIPE).stdout.decode('utf-8').splitlines())
        n_jobs_running = len(subprocess.run(['squeue', '-h', '-u', username, '-t', 'running'], stdout=subprocess.PIPE).stdout.decode('utf-8').splitlines())
        len_queued = len(str(n_jobs_queued - n_jobs_running))
        len_running = len(str(n_jobs_running))
        print(f"\r{n_jobs_queued - n_jobs_running} job(s) currently queued, {n_jobs_running} job(s) running", end=" " * (len_queued + len_running + 4))
        if n_jobs_queued == 0:
            finished_running = True
        time.sleep(1)

else:
    with Pool(processes=args.processes) as pool: # create a multiprocessing pool
        res = list(tqdm(pool.imap(run_fit, bin_iterations_seed_reaction_tuple), total=args.iterations * n_bins)) # imap(x, y) spawns processes which run a method x(y)
        # imap vs map just spawns an iterator so tqdm can make a progress bar


gather(bin_directory, config_template)
