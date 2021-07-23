#!/usr/bin/python3

"""
divide_data.py: Runs the split_mass program from halld_sim over a set of AmpTools ROOT files to divide it into mass bins.
    Run it without arguments for a more detailed description of its inputs.

Author: Nathaniel Dene Hoffman - Carnegie Mellon University - GlueX Collaboration
Creation Date: 5 July 2021
"""

import argparse
import errno
import sys
import os
import numpy as np
from pathlib import Path
parser = argparse.ArgumentParser(description="Splits up an AmpTools tree into mass bins (invariant mass of all final state particles, ignoring beam and recoil proton)")
parser.add_argument("--low", type=float, help="low bin in GeV")
parser.add_argument("--high", type=float, help="high bin in GeV")
parser.add_argument("-n", default=0, type=int, help="number of bins (run without this argument for ~40 MeV bins)")
parser.add_argument("-g", "--generated", help="path to the generated Monte Carlo AmpTools ROOT file")
parser.add_argument("-a", "--accepted", help="path to the accepted Monte Carlo  AmpTools ROOT file")
parser.add_argument("-d", "--data", help="path to the data AmpTools ROOT file")
parser.add_argument("-b", "--background", help="path to the background AmpTools ROOT file")
parser.add_argument("-t", "--tree", help="optionally specify the input tree name", default='kin')
parser.add_argument("-o", "--output", help="path to the output folder")
parser.add_argument("-c", "--config", help="path to a template amptools config file (keywords are @DATAFILE, @GENFILE, @ACCFILE, and @NIFILE)", required=True)

# Handle argument errors
if len(sys.argv) == 1: # If no arguments are provided, print the help string and exit
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
if args.high <= args.low: # high bin must be >= low bin
    print("The high bin value must be greater than the low bin value!")
    sys.exit(1)
if args.n == 0: # If no bin number is specified (or if 0 is chosen), set it to create ~40 MeV bins (as close as possible)
    span = args.high - args.low
    args.n = int(span / 0.040)

# Check to see if all the files actually exist:
config_path = Path(args.config).resolve()
if config_path.is_file():
    print(f"AmpTools Config Template: {config_path}")
else:
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.config)

if args.generated != None:
    generated_path = Path(args.generated).resolve()
    if generated_path.is_file():
        print(f"Generated file: {generated_path}")
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.generated)
if args.accepted != None:
    accepted_path = Path(args.accepted).resolve()
    if accepted_path.is_file():
        print(f"Accepted file: {accepted_path}")
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.accepted)
if args.data != None:
    data_path = Path(args.data).resolve()
    if data_path.is_file():
        print(f"Data file: {data_path}")
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.data)
if args.background != None:
    background_path = Path(args.background).resolve()
    if background_path.is_file():
        print(f"Background file: {background_path}")
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.background)
if args.output != None:
    output_dir = Path(args.output).resolve()
    if output_dir.is_dir():
        print(f"Output_Folder: {output_dir}")
    else:
        output_dir.mkdir(parents=True)
        print(f"Output_Folder: {output_dir}")


# removed for now, it seems the split_mass program sets a value by default and you can't set this
# AND set an optional tree name
#max_events = 1e9 # Maximum number of events to process, this may need to be modified in the future

output_dirs = [output_dir / f"{bin_num}" for bin_num in np.arange(args.n)]
for bin_dir in output_dirs:
    bin_dir.mkdir(exist_ok=True) # create subdirectories for each bin

bin_info_file = output_dir / "bin_info.txt"
with open(bin_info_file, 'w') as bin_info_writer:
    bin_content = np.linspace(args.low, args.high, args.n)
    lines_to_write = ["bin\tmass\n"]
    for i, bin_c in enumerate(bin_content):
        lines_to_write.append(f"{i}\t{bin_c}\n")
    bin_info_writer.writelines(lines_to_write) # write a tab-separated file containing the masses of each bin

# Process Generated
if args.generated != None:
    print("Splitting Generated Monte Carlo")
    os.system(f"split_mass {generated_path} {generated_path.stem + '_GEN_'} {args.low} {args.high} {args.n} -T {args.tree}:kin")

# Process Accepted
if args.accepted != None:
    print("Splitting Accepted Monte Carlo")
    os.system(f"split_mass {accepted_path} {accepted_path.stem + '_ACC_'} {args.low} {args.high} {args.n} -T {args.tree}:kin")

# Process Data
if args.data != None:
    print("Splitting Data")
    os.system(f"split_mass {data_path} {data_path.stem + '_DAT_'} {args.low} {args.high} {args.n} -T {args.tree}:kin")


# Process Background
if args.background != None:
    print("Splitting Background")
    os.system(f"split_mass {background_path} {background_path.stem + '_BKG_'} {args.low} {args.high} {args.n} -T {args.tree}:kin")

# Copy in a template config file to each bin and change the @TAG file paths inside to point to the proper files
for bin_num in np.arange(args.n):
    bin_files = Path('.').glob(f"*_{bin_num}.root")
    for bin_file in bin_files:
        destination = output_dir / str(bin_num) / bin_file.name
        bin_file.replace(destination)
        with open(config_path) as config:
            config_text = config.read()
            if args.generated != None:
                config_text = config_text.replace("@GENFILE", generated_path.stem + f"_{bin_num}.root")
            if args.accepted != None:
                config_text = config_text.replace("@ACCFILE", accepted_path.stem + f"_{bin_num}.root")
            if args.data != None:
                config_text = config_text.replace("@DATAFILE", data_path.stem + f"_{bin_num}.root")
            if args.background != None:
                config_text = config_text.replace("@BKGFILE", background_path.stem + f"_{bin_num}.root")
            config_text = config_text.replace("@NIFILE", f"{bin_num}_ni.txt")
        config_bin_name = config_path.stem + f"_{bin_num}.cfg"
        with open(output_dir / str(bin_num) / config_bin_name, 'w') as config_bin:
            config_bin.write(config_text)
