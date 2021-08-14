#!/usr/bin/env python3
"""
run_amptools.py: This program performs AmpTools fits over mass-binned data created by divide_data.py.
    Run it without arguments for a more detailed description of its inputs.

Author: Nathaniel Dene Hoffman - Carnegie Mellon University - GlueX Collaboration
Creation Date: 5 July 2021
"""

import argparse
import errno
import logging
import os
import shutil
import subprocess
import sys
import threading
import time
from multiprocessing import Pool
from pathlib import Path

import enlighten
import numpy as np
from colorama import Fore
from simple_term_menu import TerminalMenu
from tqdm import tqdm


def config_menu(root):
    """

    :param root: 

    """
    config_menu_title = "Select a fit_results.txt file:"
    configs = [p for p in root.glob("*.cfg")]
    config_menu_items = ["Cancel"] + [
        f.name + (" " * (40 - len(f.name))) +
        ("(bootstrapped)" if
         (f.parent / (f.stem + "::bootstrap.txt")).exists() else
         ("(fit)" if
          (f.parent / (f.stem + "::fit_results.txt")).exists() else "*"))
        for f in configs
    ]
    config_menu_cursor = "> "
    config_menu_cursor_style = ("fg_red", "bold")
    config_menu_style = ("bg_black", "fg_green")
    config_menu = TerminalMenu(menu_entries=config_menu_items,
                               title=config_menu_title,
                               menu_cursor=config_menu_cursor,
                               menu_cursor_style=config_menu_cursor_style,
                               menu_highlight_style=config_menu_style)
    selection_index = config_menu.show()
    if selection_index == 0:
        print("No config file chosen")
        sys.exit(1)
    else:
        return root / configs[selection_index - 1]


def run_pool(tup):
    """

    :param tup: 

    """
    global log_dir
    subprocess.run(
        f"python3 ../run_fit.py {int(tup[0])} {int(tup[1])} {int(tup[2])} {str(tup[3])} {str(log_dir)} {str(tup[4])} {str(tup[5])}"
        .split(),
        stdout=subprocess.PIPE)


def update_counters(active_jobs, bin_directory, configstem, total_counter,
                    total_call_limit, total_failed, counters,
                    counters_call_limit, counters_failed):
    """

    :param active_jobs: 
    :param bin_directory: 
    :param configstem: 
    :param total_counter: 
    :param total_call_limit: 
    :param total_failed: 
    :param counters: 
    :param counters_call_limit: 
    :param counters_failed: 

    """
    global args
    active_indices = np.argwhere(active_jobs != 0)
    while (len(active_indices) != 0):
        for index in active_indices:
            bin_str = str(index[0])
            iteration_str = str(index[1])
            fit_directory = bin_directory / bin_str / iteration_str
            if args.bootstrap:
                potential_fit_files = list(
                    fit_directory.glob(f"{configstem}*.bootstrap"))
                potential_results_file = fit_directory / f"{configstem}::bootstrap.txt"
            else:
                potential_fit_files = list(
                    fit_directory.glob(f"{configstem}*.fit"))
                potential_results_file = fit_directory / f"{configstem}::fit_results.txt"
            if len(potential_fit_files) != 0 and potential_results_file.exists(
            ):
                fit_file = potential_fit_files[0]
                active_jobs[index[0], index[1]] = 0
                if "CONVERGED" in fit_file.name:
                    counters[index[0]].update()
                    total_counter.update()
                elif "CALL_LIMIT" in fit_file.name:
                    counters_call_limit[index[0]].update()
                    total_call_limit.update()
                elif "FAILED" in fit_file.name:
                    counters_failed[index[0]].update()
                    total_failed.update()
        active_indices = np.argwhere(active_jobs != 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Runs AmpTools fits on each mass bin")
    parser.add_argument("-d",
                        "--directory",
                        required=True,
                        help="the input directory (output of divide_data.py)")
    parser.add_argument("-p",
                        "--processes",
                        default=60,
                        type=int,
                        help="number of processes to spawn")
    parser.add_argument("-s",
                        "--seed",
                        default=1,
                        type=int,
                        help="set the seed (default = 1)")
    parser.add_argument("-i",
                        "--iterations",
                        default=1,
                        type=int,
                        help="set the number of fits to perform for each bin")
    parser.add_argument(
        "--parallel",
        default="SLURM",
        help=
        "method of parallelization, options are SLURM (default) and Pool (Python Multiprocessing)"
    )
    parser.add_argument(
        "--bootstrap",
        action='store_true',
        help=
        "use bootstrapping (must have run a fit already and run bootstrap.py)")
    parser.add_argument("--rerun",
                        action='store_true',
                        help="rerun existing fit")
    queue_group = parser.add_mutually_exclusive_group()
    queue_group.add_argument('-r',
                             '--red',
                             action='store_true',
                             help="run on red queue (default)")
    queue_group.add_argument('-g',
                             '--green',
                             action='store_true',
                             help="run on green queue")
    parser.add_argument("-v",
                        "--verbose",
                        action='store_true',
                        help="verbose output")
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

    config_template = config_menu(bin_directory)

    reaction = "REACTION"
    with open(config_template, 'r') as config_temp:
        lines = config_temp.readlines()
        for line in lines:
            if line.startswith("reaction"):
                reaction = line.split()[1]
    print(f"Reaction Name: {reaction}")

    logging.basicConfig(filename=f"{reaction}.log",
                        filemode='a',
                        format="%(asctime)s - %(levelname)s:%(message)s",
                        datefmt="%d/%m/%y %H:%M:%S")

    logging.info("Beginning AmpTools Fitting Procedure")

    log_dir = Path("logs").resolve()
    log_dir.mkdir(exist_ok=True)

    logging.info(f"Setting seed to {args.seed}")
    np.random.seed(args.seed)
    seeds = [np.random.randint(1, high=100000) for _ in range(args.iterations)]
    n_bins = len(
        [bin_dir for bin_dir in bin_directory.glob("*") if bin_dir.is_dir()])
    logging.info(
        f"Using {n_bins} bins and {args.iterations} iterations per bin ({n_bins * args.iterations} total fits)"
    )
    run_tuple = [(i, j, seeds[j], reaction, args.bootstrap,
                  config_template.stem)
                 for i in range(n_bins)
                 for j in range(args.iterations)]

    if args.bootstrap:
        if (bin_directory /
                f"{config_template.stem}::fit_results.txt").exists():
            logging.info(
                "Bootstrapping enabled, generating new configuration file")
            os.system(
                f"python3 bootstrap.py {str(bin_directory)} {config_template.stem}"
            )
        else:
            print(
                "You need to run a fit with that configuration file before you can bootstrap it!"
            )
            logging.error(
                "No configuration file found at {str(bin_directory)}/{config_template.stem}::fit_results.txt"
            )
            sys.exit(1)

    os.chdir(str(bin_directory))
    if args.parallel == "SLURM":
        logging.debug("Using SLURM to distribute fit jobs")
        if args.green:
            cpu_memory = 1590
            threads = 4
            queue = "green"
        else:
            cpu_memory = 1990
            threads = 4
            queue = "red"
        memory = threads * cpu_memory
        np.random.shuffle(run_tuple)
        if args.rerun:
            logging.warning(
                "Rerun flag has been specified, existing fit files will be deleted!"
            )
            for tup in run_tuple:
                bin_number, iteration, seed, reaction, bootstrap, configstem = tup
                fit_file = (
                    Path(".") /
                    f"{bin_number}/{iteration}/{configstem}::fit_results.txt"
                ).resolve()
                bootstrap_file = (
                    Path(".") /
                    f"{bin_number}/{iteration}/{configstem}::bootstrap.txt"
                ).resolve()
                if bootstrap:
                    if bootstrap_file.exists():
                        bootstrap_file.unlink()
                else:
                    if fit_file.exists():
                        fit_file.unlink()
                    if bootstrap_file.exists():
                        bootstrap_file.unlink()
        logging.info("Submitting Jobs to Cluster")
        if args.verbose:
            manager = enlighten.get_manager()
            active_jobs = np.array(
                [[1 for _ in range(args.iterations)] for _ in range(n_bins)])
            total_counter = manager.counter(total=(args.iterations * n_bins),
                                            color='white',
                                            desc="Total")
            total_failed = total_counter.add_subcounter(color='red')
            total_call_limit = total_counter.add_subcounter(color='yellow')
            counters = [
                manager.counter(
                    total=args.iterations,
                    desc=(f"Bin {bin_num}" + " " * (4 - len(str(bin_num)))),
                    color='green',
                    bar_format=
                    '{desc}{desc_pad}{percentage:3.0f}%|{bar}| {count:{len_total}d}/{total:d}'
                ) for bin_num in range(n_bins)
            ]
            counters_failed = [
                counter.add_subcounter(color='red') for counter in counters
            ]
            counters_call_limit = [
                counter.add_subcounter(color='yellow') for counter in counters
            ]
            status_bar = manager.status_bar('Submitting SLURM Jobs',
                                            color='white_on_black',
                                            justify=enlighten.Justify.CENTER)
            bars = threading.Thread(
                target=update_counters,
                args=(active_jobs, bin_directory, config_template.stem,
                      total_counter, total_call_limit, total_failed, counters,
                      counters_call_limit, counters_failed))
            bars.start()

            for i, tup in enumerate(run_tuple):
                bin_number, iteration, seed, reaction, bootstrap, configstem = tup
                fit_file = (
                    Path(".") /
                    f"{bin_number}/{iteration}/{configstem}::fit_results.txt"
                ).resolve()
                bootstrap_file = (
                    Path(".") /
                    f"{bin_number}/{iteration}/{configstem}::bootstrap.txt"
                ).resolve()
                log_out = log_dir / f"{reaction}_{bin_number}_{iteration}_SLURM.out"
                log_err = log_dir / f"{reaction}_{bin_number}_{iteration}_SLURM.err"
                if not fit_file.exists() and not (bootstrap and
                                                  bootstrap_file.exists()):
                    status_bar.update(
                        f"Submitting SLURM Jobs: {i}/{len(run_tuple)}")
                    os.system(
                        f"sbatch --job-name={reaction}_{bin_number}_{iteration} --ntasks={threads} --partition={queue} --mem={memory} --time=30:00 --output={str(log_out)} --error={str(log_err)} --quiet ../run_fit_slurm.csh {bin_number} {iteration} {seed} {reaction} {str(log_dir)} {str(bootstrap)} {str(configstem)}"
                    )
                    time.sleep(1)

            finished_running = False
            username = os.getlogin()
            while not finished_running:
                n_jobs_queued = len(
                    subprocess.run(['squeue', '-h', '-u', username],
                                   stdout=subprocess.PIPE).stdout.decode(
                                       'utf-8').splitlines())
                n_jobs_running = len(
                    subprocess.run(
                        ['squeue', '-h', '-u', username, '-t', 'running'],
                        stdout=subprocess.PIPE).stdout.decode(
                            'utf-8').splitlines())
                len_queued = len(str(n_jobs_queued - n_jobs_running))
                len_running = len(str(n_jobs_running))
                status_bar.update(
                    f"{n_jobs_queued - n_jobs_running} job(s) currently queued, {n_jobs_running} job(s) running"
                )
                if n_jobs_queued == 0:
                    finished_running = True
                time.sleep(1)
            manager.stop()
        else:
            manager = enlighten.get_manager()
            total_counter = manager.counter(total=(args.iterations * n_bins),
                                            color='white',
                                            desc="Total")
            status_bar = manager.status_bar('Submitting SLURM Jobs',
                                            color='white_on_black',
                                            justify=enlighten.Justify.CENTER)
            username = os.getlogin()
            for tup in run_tuple:
                bin_number, iteration, seed, reaction, bootstrap, configstem = tup
                fit_file = (
                    Path(".") /
                    f"{bin_number}/{iteration}/{configstem}::fit_results.txt"
                ).resolve()
                bootstrap_file = (
                    Path(".") /
                    f"{bin_number}/{iteration}/{configstem}::bootstrap.txt"
                ).resolve()
                log_out = log_dir / f"{reaction}_{bin_number}_{iteration}_SLURM.out"
                log_err = log_dir / f"{reaction}_{bin_number}_{iteration}_SLURM.err"
                total_counter.update()
                if not fit_file.exists() and not (bootstrap and
                                                  bootstrap_file.exists()):
                    os.system(
                        f"sbatch --job-name={reaction}_{bin_number}_{iteration} --ntasks={threads} --partition={queue} --mem={memory} --time=30:00 --output={str(log_out)} --error={str(log_err)} --quiet ../run_fit_slurm.csh {bin_number} {iteration} {seed} {reaction} {str(log_dir)} {str(bootstrap)} {str(configstem)}"
                    )

                    n_jobs_queued = len(
                        subprocess.run(['squeue', '-h', '-u', username],
                                       stdout=subprocess.PIPE).stdout.decode(
                                           'utf-8').splitlines())
                    n_jobs_running = len(
                        subprocess.run(
                            ['squeue', '-h', '-u', username, '-t', 'running'],
                            stdout=subprocess.PIPE).stdout.decode(
                                'utf-8').splitlines())
                    len_queued = len(str(n_jobs_queued - n_jobs_running))
                    len_running = len(str(n_jobs_running))
                    status_bar.update(
                        f"{n_jobs_queued - n_jobs_running} job(s) currently queued, {n_jobs_running} job(s) running"
                    )
                    time.sleep(1)
            total_counter.close()

            finished_running = False
            while not finished_running:
                n_jobs_queued = len(
                    subprocess.run(['squeue', '-h', '-u', username],
                                   stdout=subprocess.PIPE).stdout.decode(
                                       'utf-8').splitlines())
                n_jobs_running = len(
                    subprocess.run(
                        ['squeue', '-h', '-u', username, '-t', 'running'],
                        stdout=subprocess.PIPE).stdout.decode(
                            'utf-8').splitlines())
                len_queued = len(str(n_jobs_queued - n_jobs_running))
                len_running = len(str(n_jobs_running))
                status_bar.update(
                    f"{n_jobs_queued - n_jobs_running} job(s) currently queued, {n_jobs_running} job(s) running"
                )
                if n_jobs_queued == 0:
                    finished_running = True
                time.sleep(1)
            manager.stop()
        logger.info("All jobs have been processed")
        logger.info("Gathering results")
        os.system(
            f"python3 ../gather.py -d {bin_directory} -c {config_template} {'--bootstrap' if args.bootstrap else ''} -n {args.iterations}"
        )

    elif args.parallel == "Pool":
        logger.info(
            f"Running jobs in multiprocessing pool with {args.processes} simultaneous processes"
        )
        np.random.shuffle(run_tuple)
        if args.rerun:
            logging.warning(
                "Rerun flag has been specified, existing fit files will be deleted!"
            )
            for tup in run_tuple:
                bin_number, iteration, seed, reaction, bootstrap, configstem = tup
                fit_file = (
                    Path(".") /
                    f"{bin_number}/{iteration}/{configstem}::fit_results.txt"
                ).resolve()
                bootstrap_file = (
                    Path(".") /
                    f"{bin_number}/{iteration}/{configstem}::bootstrap.txt"
                ).resolve()
                if bootstrap:
                    if bootstrap_file.exists():
                        bootstrap_file.unlink()
                else:
                    if fit_file.exists():
                        fit_file.unlink()
                    if bootstrap_file.exists():
                        bootstrap_file.unlink()
        if args.verbose:
            manager = enlighten.get_manager()
            active_jobs = np.array(
                [[1 for _ in range(args.iterations)] for _ in range(n_bins)])
            total_counter = manager.counter(total=(args.iterations * n_bins),
                                            color='white',
                                            desc="Total")
            total_failed = total_counter.add_subcounter(color='red')
            total_call_limit = total_counter.add_subcounter(color='yellow')
            counters = [
                manager.counter(
                    total=args.iterations,
                    desc=(f"Bin {bin_num}" + " " * (4 - len(str(bin_num)))),
                    color='green',
                    bar_format=
                    '{desc}{desc_pad}{percentage:3.0f}%|{bar}| {count:{len_total}d}/{total:d}'
                ) for bin_num in range(n_bins)
            ]
            counters_failed = [
                counter.add_subcounter(color='red') for counter in counters
            ]
            counters_call_limit = [
                counter.add_subcounter(color='yellow') for counter in counters
            ]
            status_bar = manager.status_bar('Running Pool Jobs',
                                            color='white_on_black',
                                            justify=enlighten.Justify.CENTER)
            bars = threading.Thread(
                target=update_counters,
                args=(active_jobs, bin_directory, config_template.stem,
                      total_counter, total_call_limit, total_failed, counters,
                      counters_call_limit, counters_failed))
            bars.start()
            with Pool(processes=args.processes) as pool:
                res = pool.map(run_pool, run_tuple)
            manager.stop()
        else:
            with Pool(processes=args.processes) as pool:
                res = list(
                    tqdm(pool.imap(run_pool, run_tuple),
                         total=args.iterations * n_bins))
        logger.info("All jobs have been processed")
        logger.info("Gathering results")
        os.system(
            f"python3 ../gather.py -d {bin_directory} -c {config_template} {'--bootstrap' if args.bootstrap else ''} -n {args.iterations}"
        )
    else:
        logger.error(
            "Parallelization method {args.parallel} is not currently supported!"
        )
        print("Please select a supported parallelization method!")
        sys.exit(1)
