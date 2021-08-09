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
from simple_term_menu import TerminalMenu
from colorama import Fore
import time


def config_menu(root):
    config_menu_title = "Select a fit_results.txt file:"
    configs = [p for p in root.glob("*.cfg")]
    config_menu_items = ["Cancel"] + [f.name + (" " * (40 - len(f.name))) + ("(bootstrapped)" if (f.parent / (f.stem + "::bootstrap.txt")).exists() else ("(fit)" if (f.parent / (f.stem + "::fit_results.txt")).exists() else "*")) for f in configs] 
    config_menu_cursor = "> "
    config_menu_cursor_style = ("fg_red", "bold")
    config_menu_style = ("bg_black", "fg_green")
    config_menu = TerminalMenu(
        menu_entries=config_menu_items,
        title=config_menu_title,
        menu_cursor=config_menu_cursor,
        menu_cursor_style=config_menu_cursor_style,
        menu_highlight_style=config_menu_style
    )
    selection_index = config_menu.show()
    if selection_index == 0:
        print("No config file chosen")
        sys.exit(1)
    else:
        return root / configs[selection_index - 1]

def run_pool(tup):
    global log_dir
    os.system(f"python3 ../run_fit.py {int(tup[0])} {int(tup[1])} {int(tup[2])} {str(tup[3])} {str(log_dir)} {str(tup[4])} {str(tup[5])}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Runs AmpTools fits on each mass bin")
    parser.add_argument("-d", "--directory", required=True, help="the input directory (output of divide_data.py)")
    parser.add_argument("-p", "--processes", default=60, type=int, help="number of processes to spawn")
    parser.add_argument("-s", "--seed", default=1, type=int, help="set the seed (default = 1)")
    parser.add_argument("-i", "--iterations", default=1, type=int, help="set the number of fits to perform for each bin")
    parser.add_argument("--parallel", default="SLURM", help="method of parallelization, options are SLURM (default) and Pool (Python Multiprocessing)")
    parser.add_argument("--bootstrap", action='store_true', help="use bootstrapping (must have run a fit already and run bootstrap.py)")
    parser.add_argument("--rerun", action='store_true', help="rerun existing fit")
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

    config_template = config_menu(bin_directory)

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
    bin_iterations_seed_reaction_bootstrap_configstem_tuple = [(i, j, seeds[j], reaction, args.bootstrap, config_template.stem) for i in range(n_bins) for j in range(args.iterations)]
    # tuple contains the bin number, the iteration number, a seed specific to the iteration (iteration "j" across all bins will have the same seed, so in theory
    # it should also have the same starting parameters, in case that is important later on), and the reaction name

    if args.bootstrap:
        if (bin_directory / f"{config_template.stem}::fit_results.txt").exists():
            os.system(f"python3 bootstrap.py {str(bin_directory)} {config_template.stem}") # run bootstrap script (assuming a fit has already been performed)
        else:
            print("You need to run a fit with that configuration file before you can bootstrap it!")
            sys.exit(1)

    os.chdir(str(bin_directory)) # cd into the directory containing bin subdirectories
    if args.parallel == "SLURM":
        for tup in tqdm(bin_iterations_seed_reaction_bootstrap_configstem_tuple):
            bin_number, iteration, seed, reaction, bootstrap, configstem = tup
            fit_file = Path(".") / f"{bin_number}/{iteration}/{configstem}::fit_results.txt"
            bootstrap_file = Path(".") / f"{bin_number}/{iteration}/{configstem}::bootstrap.txt"
            log_out = log_dir / f"{reaction}_{bin_number}_{iteration}_SLURM.out"
            log_err = log_dir / f"{reaction}_{bin_number}_{iteration}_SLURM.err"
            if args.rerun or (not fit_file.exists()) and (not (bootstrap and bootstrap.exists())):
                os.system(f"sbatch --job-name={reaction}_{bin_number}_{iteration} --ntasks={threads} --partition={queue} --mem={memory} --time=30:00 --output={str(log_out)} --error={str(log_err)} --quiet ../run_fit_slurm.csh {bin_number} {iteration} {seed} {reaction} {str(log_dir)} {str(bootstrap)} {str(configstem)}")
                time.sleep(1)

        # Wait for all jobs to finish before gathering results
        finished_running = False
        username = os.getlogin() # get username for checking SLURM
        while not finished_running:
            n_jobs_queued = len(subprocess.run(['squeue', '-h', '-u', username], stdout=subprocess.PIPE).stdout.decode('utf-8').splitlines())
            n_jobs_running = len(subprocess.run(['squeue', '-h', '-u', username, '-t', 'running'], stdout=subprocess.PIPE).stdout.decode('utf-8').splitlines())
            len_queued = len(str(n_jobs_queued - n_jobs_running))
            len_running = len(str(n_jobs_running))
            print(f"\r{n_jobs_queued - n_jobs_running} job(s) currently queued, {n_jobs_running} job(s) running", end=" " * (len_queued + len_running + 4))
            if n_jobs_queued == 0:
                finished_running = True
            time.sleep(1)
        os.system(f"python3 ../gather.py -d {bin_directory} -c {config_template} {'--bootstrap' if args.bootstrap else ''}")

    elif args.parallel == "Pool":
        with Pool(processes=args.processes) as pool: # create a multiprocessing pool
            res = list(tqdm(pool.imap(run_pool, bin_iterations_seed_reaction_bootstrap_configstem_tuple), total=args.iterations * n_bins)) # imap(x, y) spawns processes which run a method x(y)
            # imap vs map just spawns an iterator so tqdm can make a progress bar
        os.system(f"python3 ../gather.py -d {bin_directory} -c {config_template} {'--bootstrap' if args.bootstrap else ''}")
    else:
        print("Please select a supported parallelization method!")
