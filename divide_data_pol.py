#!/usr/bin/env python3
"""
divide_data_pol.py: Runs the split_mass program from halld_sim over a set of AmpTools
ROOT files to divide it into mass bins. Run it without arguments for a more detailed
description of its inputs.

Author: Nathaniel Dene Hoffman - Carnegie Mellon University - GlueX Collaboration
Creation Date: 5 July 2021
"""

import argparse
import errno
import os
import shutil
import sys
from pathlib import Path

import numpy as np


def gen_config(config_template, working_dir):
    """gen_config creates configuration files from a template by replacing tags which
    begin with '@' and describe paths to various files.

    :param config_template: path to configuration template
    :param working_dir: path to bin directory

    """
    main_config_destination = working_dir / config_template.name
    shutil.copy(config_template, main_config_destination)
    for bin_dir in [
            bin_dir for bin_dir in working_dir.glob("*") if bin_dir.is_dir()
    ]:
        bin_num = int(bin_dir.name)
        with open(config_template) as config:
            config_text = config.read()
            if "@GENFILE" in config_text:
                generated_path_AMO = [
                    p for p in bin_dir.glob(f"*_GEN__{bin_num}.root")
                    if "AMO" in p.name
                ][0]
                generated_path_000 = [
                    p for p in bin_dir.glob(f"*_GEN__{bin_num}.root")
                    if "PARA_0" in p.name
                ][0]
                generated_path_045 = [
                    p for p in bin_dir.glob(f"*_GEN__{bin_num}.root")
                    if "PERP_45" in p.name
                ][0]
                generated_path_090 = [
                    p for p in bin_dir.glob(f"*_GEN__{bin_num}.root")
                    if "PERP_90" in p.name
                ][0]
                generated_path_135 = [
                    p for p in bin_dir.glob(f"*_GEN__{bin_num}.root")
                    if "PARA_135" in p.name
                ][0]
                config_text = config_text.replace("@GENFILE_AMO",
                                                  generated_path_AMO.name)
                config_text = config_text.replace("@GENFILE_000",
                                                  generated_path_000.name)
                config_text = config_text.replace("@GENFILE_045",
                                                  generated_path_045.name)
                config_text = config_text.replace("@GENFILE_090",
                                                  generated_path_090.name)
                config_text = config_text.replace("@GENFILE_135",
                                                  generated_path_135.name)
            if "@ACCFILE" in config_text:
                accepted_path_AMO = [
                    p for p in bin_dir.glob(f"*_ACCEPT__{bin_num}.root")
                    if "AMO" in p.name
                ][0]
                accepted_path_000 = [
                    p for p in bin_dir.glob(f"*_ACCEPT__{bin_num}.root")
                    if "PARA_0" in p.name
                ][0]
                accepted_path_045 = [
                    p for p in bin_dir.glob(f"*_ACCEPT__{bin_num}.root")
                    if "PERP_45" in p.name
                ][0]
                accepted_path_090 = [
                    p for p in bin_dir.glob(f"*_ACCEPT__{bin_num}.root")
                    if "PERP_90" in p.name
                ][0]
                accepted_path_135 = [
                    p for p in bin_dir.glob(f"*_ACCEPT__{bin_num}.root")
                    if "PARA_135" in p.name
                ][0]
                config_text = config_text.replace("@ACCFILE_AMO",
                                                  accepted_path_AMO.name)
                config_text = config_text.replace("@ACCFILE_000",
                                                  accepted_path_000.name)
                config_text = config_text.replace("@ACCFILE_045",
                                                  accepted_path_045.name)
                config_text = config_text.replace("@ACCFILE_090",
                                                  accepted_path_090.name)
                config_text = config_text.replace("@ACCFILE_135",
                                                  accepted_path_135.name)
            if "@DATAFILE" in config_text:
                data_path_AMO = [
                    p for p in bin_dir.glob(f"*_DATA__{bin_num}.root")
                    if "AMO" in p.name
                ][0]
                data_path_000 = [
                    p for p in bin_dir.glob(f"*_DATA__{bin_num}.root")
                    if "PARA_0" in p.name
                ][0]
                data_path_045 = [
                    p for p in bin_dir.glob(f"*_DATA__{bin_num}.root")
                    if "PERP_45" in p.name
                ][0]
                data_path_090 = [
                    p for p in bin_dir.glob(f"*_DATA__{bin_num}.root")
                    if "PERP_90" in p.name
                ][0]
                data_path_135 = [
                    p for p in bin_dir.glob(f"*_DATA__{bin_num}.root")
                    if "PARA_135" in p.name
                ][0]
                config_text = config_text.replace("@DATAFILE_AMO",
                                                  data_path_AMO.name)
                config_text = config_text.replace("@DATAFILE_000",
                                                  data_path_000.name)
                config_text = config_text.replace("@DATAFILE_045",
                                                  data_path_045.name)
                config_text = config_text.replace("@DATAFILE_090",
                                                  data_path_090.name)
                config_text = config_text.replace("@DATAFILE_135",
                                                  data_path_135.name)
            if "@BKGFILE" in config_text:
                background_path_AMO = [
                    p for p in bin_dir.glob(f"*_BKG__{bin_num}.root")
                    if "AMO" in p.name
                ][0]
                background_path_000 = [
                    p for p in bin_dir.glob(f"*_BKG__{bin_num}.root")
                    if "PARA_0" in p.name
                ][0]
                background_path_045 = [
                    p for p in bin_dir.glob(f"*_BKG__{bin_num}.root")
                    if "PERP_45" in p.name
                ][0]
                background_path_090 = [
                    p for p in bin_dir.glob(f"*_BKG__{bin_num}.root")
                    if "PERP_90" in p.name
                ][0]
                background_path_135 = [
                    p for p in bin_dir.glob(f"*_BKG__{bin_num}.root")
                    if "PARA_135" in p.name
                ][0]
                config_text = config_text.replace("@BKGFILE_AMO",
                                                  background_path_AMO.name)
                config_text = config_text.replace("@BKGFILE_000",
                                                  background_path_000.name)
                config_text = config_text.replace("@BKGFILE_045",
                                                  background_path_045.name)
                config_text = config_text.replace("@BKGFILE_090",
                                                  background_path_090.name)
                config_text = config_text.replace("@BKGFILE_135",
                                                  background_path_135.name)
            config_text = config_text.replace(
                "@NIFILE_AMO", f"{config_template.stem}_{bin_num}_ni_AMO.txt")
            config_text = config_text.replace(
                "@NIFILE_000", f"{config_template.stem}_{bin_num}_ni_000.txt")
            config_text = config_text.replace(
                "@NIFILE_045", f"{config_template.stem}_{bin_num}_ni_045.txt")
            config_text = config_text.replace(
                "@NIFILE_090", f"{config_template.stem}_{bin_num}_ni_090.txt")
            config_text = config_text.replace(
                "@NIFILE_135", f"{config_template.stem}_{bin_num}_ni_135.txt")
        config_bin_name = config_template.stem + f"_{bin_num}.cfg"
        with open(output_dir / str(bin_num) / config_bin_name,
                  'w') as config_bin:
            config_bin.write(config_text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=("Splits up an AmpTools tree into mass bins (invariant mass of all "
                     "final state particles, ignoring beam and recoil proton).\n\nRun "
                     "with only -o <output dir> -c <config> to add a config file to an "
                     "existing output directory"))
    parser.add_argument("--low", type=float, help="low bin in GeV")
    parser.add_argument("--high", type=float, help="high bin in GeV")
    parser.add_argument(
        "-n",
        default=0,
        type=int,
        help="number of bins (run without this argument for ~40 MeV bins)")
    parser.add_argument(
        "-g",
        "--generated",
        help=
        "path to the folder containing generated Monte Carlo AmpTools ROOT files"
    )
    parser.add_argument("-tG",
                        "--treeGenerated",
                        help="optionally specify the input Generated tree name",
                        default='kin')
    parser.add_argument(
        "-a",
        "--accepted",
        help=
        "path to the folder containing accepted Monte Carlo  AmpTools ROOT files"
    )
    parser.add_argument("-tA",
                        "--treeAccepted",
                        help="optionally specify the input Accepted tree name",
                        default='kin')
    parser.add_argument(
        "-d",
        "--data",
        help="path to the folder containing data AmpTools ROOT files")
    parser.add_argument("-tD",
                        "--treeData",
                        help="optionally specify the input Data tree name",
                        default='kin')
    parser.add_argument(
        "-b",
        "--background",
        help="path to the folder containing background AmpTools ROOT files")
    parser.add_argument(
        "-tB",
        "--treeBackground",
        help="optionally specify the input Bakcground tree name",
        default='kin')
    parser.add_argument("-t",
                        "--tree",
                        help="optionally specify the input tree name",
                        default='kin')
    parser.add_argument("-o",
                        "--output",
                        help="path to the output folder",
                        required=True)
    parser.add_argument("-c",
                        "--config",
                        help="path to a template AmpTools config file",
                        required=True)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    if len(sys.argv) == 5 and (args.output != None) and (args.config != None):
        # if only config and output are provided, add the new config instead of
        # re-splitting all the data files again
        output_dir = Path(args.output).resolve()
        config_path = Path(args.config).resolve()
        print((f"Adding config file {str(config_path)} to existing binned analysis "
               f"directory {str(output_dir)}"))
        assert config_path.is_file(
        ), "Config file {str(config_path)} either does not exist or is not a file!"
        assert output_dir.is_dir(
        ), f"{output_dir} does not exist or is not a directory!"
        print(f"AmpTools Config Template: {config_path}")
        print(f"Output_Folder: {output_dir}")
        gen_config(config_path, output_dir)
        sys.exit(1)
    if args.high <= args.low:
        print("The high bin value must be greater than the low bin value!")
        sys.exit(1)
    if args.n == 0:
        span = args.high - args.low
        args.n = int(span / 0.040)

    config_path = Path(args.config).resolve()
    if config_path.is_file():
        print(f"AmpTools Config Template: {config_path}")
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                args.config)
    # Collect files to divide:
    if args.generated != None:
        generated_dir = Path(args.generated).resolve()
        if generated_dir.is_dir():
            generated_glob = generated_dir.glob("*.root")
            generated_print = "\n".join([f"\t{f.name}" for f in generated_glob])
            print(f"Generated files:\n---\n{generated_print}\n---\n")
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    args.generated)
    if args.accepted != None:
        accepted_dir = Path(args.accepted).resolve()
        if accepted_dir.is_dir():
            accepted_glob = accepted_dir.glob("*.root")
            accepted_print = "\n".join([f"\t{f.name}" for f in accepted_glob])
            print(f"Accepted files:\n---\n{accepted_print}\n---\n")
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    args.accepted)
    if args.data != None:
        data_dir = Path(args.data).resolve()
        if data_dir.is_dir():
            data_glob = data_dir.glob("*.root")
            data_print = "\n".join([f"\t{f.name}" for f in data_glob])
            print(f"Data files:\n---\n{data_print}\n---\n")
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    args.data)
    if args.background != None:
        background_dir = Path(args.background).resolve()
        if background_dir.is_dir():
            background_glob = background_dir.glob("*.root")
            background_print = "\n".join(
                [f"\t{f.name}" for f in background_glob])
            print(f"Background files:\n---\n{background_print}\n---\n")
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    args.background)
    if args.output != None:
        output_dir = Path(args.output).resolve()
        if output_dir.is_dir():
            print(f"Output_Folder: {output_dir}")
        else:
            output_dir.mkdir(parents=True)
            print(f"Output_Folder: {output_dir}")

    output_dirs = [output_dir / f"{bin_num}" for bin_num in np.arange(args.n)]
    for bin_dir in output_dirs:
        bin_dir.mkdir(exist_ok=True)

    bin_info_file = output_dir / "bin_info.txt"
    with open(bin_info_file, 'w') as bin_info_writer:
        bin_content = np.linspace(args.low, args.high, args.n)
        lines_to_write = ["bin\tmass\n"]
        for i, bin_c in enumerate(bin_content):
            lines_to_write.append(f"{i}\t{bin_c}\n")
        bin_info_writer.writelines(lines_to_write)

    tmp = output_dir / "tmp" # do all the split_mass commands in a temporary directory
    tmp.mkdir(exist_ok=True)
    os.chdir(str(tmp))

    if args.generated != None:
        print("Splitting Generated Monte Carlo")
        generated_glob = generated_dir.glob("*.root")
        for generated_path in generated_glob:
            os.system((f"split_mass {generated_path} {generated_path.stem + '_GEN_'} "
                       f"{args.low} {args.high} {args.n} -T {args.treeGenerated}:kin"))

    if args.accepted != None:
        print("Splitting Accepted Monte Carlo")
        accepted_glob = accepted_dir.glob("*.root")
        for accepted_path in accepted_glob:
            os.system((f"split_mass {accepted_path} {accepted_path.stem + '_ACCEPT_'} "
                       f"{args.low} {args.high} {args.n} -T {args.treeAccepted}:kin"))

    if args.data != None:
        print("Splitting Data")
        data_glob = data_dir.glob("*.root")
        for data_path in data_glob:
            os.system((f"split_mass {data_path} {data_path.stem + '_DATA_'} "
                       f"{args.low} {args.high} {args.n} -T {args.treeData}:kin"))

    if args.background != None:
        print("Splitting Background")
        background_glob = background_dir.glob("*.root")
        for background_path in background_glob:
            os.system((f"split_mass {background_path} {background_path.stem + '_BKG_'} "
                       f"{args.low} {args.high} {args.n} -T {args.treeBackground}:kin"))

    os.chdir("..")

    for bin_num in np.arange(args.n):
        bin_files = tmp.glob(f"*_{bin_num}.root")
        for bin_file in bin_files:
            destination = output_dir / str(bin_num) / bin_file.name
            bin_file.replace(destination) # move ROOT files to their proper bins

    tmp.rmdir() # get rid of tmp
    gen_config(config_path, output_dir) # generate configurations for each bin
