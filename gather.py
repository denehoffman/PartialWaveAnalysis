#!/usr/bin/python3

"""
gather_fits.py: This program collects data from AmpTools fits and stores them in a single location.
    This is useful if you need to collect additional fit data or the program run_amptools.py exited
    before it got to this step for some reason. Note that this code is run at the end of run_amptools.py.
    Run it without arguments for a more detailed description of its inputs.
    Author: Nathaniel Dene Hoffman - Carnegie Mellon University - GlueX Collaboration
    Creation Date: 13 July 2021
"""

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
    with open(config_file, 'r') as config:
        config_lines = config.readlines() # read in the lines from the config template file
        parameters = []
        func_names = []
        for line in config_lines:
            if line.startswith("amplitude"): # find amplitude lines
                param_name = line.split()[1] # get the parameter name (KsKs::PositiveIm::S0+, for example)
                func_name = line.split()[2] # get function type (Zlm or TwoPSAngles)
                parameters.append(param_name) 
                func_names.append(func_name)
        param_header_list = [] # make a tab-separated header for each parameter
        for i, parameter in enumerate(parameters):
            if func_names[i] == "Zlm":
                if "Re" in parameter: # we use both the Re and Im parameters to calculate one intensity, so ignore the "Im" ones in the header
                    param_header_list.append(parameter.replace("Re", ""))
                    param_header_list.append(parameter.replace("Re", "") +"_err")
            elif func_names[i] == "TwoPSAngles":
                param_header_list.append(parameter)
                param_header_list.append(parameter + "_err")
            else:
                print(f"Function not found: {func_names[i]}")
        param_header_list.append("total_intensity")
        param_header_list.append("total_intensity_err")
        param_header_list.append("likelihood")
        param_header = "\t".join(param_header_list)
    with open(output_dir / "fit_results.txt", 'w') as out_file:
        out_file.write(f"Bin\tIteration\t{param_header}\n") # print the header to the output file
        for bin_dir in [bindir for bindir in output_dir.glob("*") if bindir.is_dir()]: # for each bin subdirectory
            bin_num_string = bin_dir.name
            for iteration_dir in [iterdir for iterdir in bin_dir.glob("*") if iterdir.is_dir()]: # for each iteration subdirectory
                iteration_num_string = iteration_dir.name
                fit_file = [fit for fit in iteration_dir.glob("*.fit")][0].resolve() # should be only one .fit file in this directory
                if "CONVERGED" in fit_file.name: # only collect converged fits
                    process = subprocess.run(['get_fit_results', str(fit_file), func_names[0], *parameters], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                    out_file.write(f"{bin_num_string}\t{iteration_num_string}\t{process.stdout}\n") # write fit results to output file (in no particular row order)


"""
Script begins here:
"""
parser = argparse.ArgumentParser(description="Runs AmpTools fits on each mass bin")
parser.add_argument("-d", "--directory", required=True, help="the input directory (output of divide_data.py)")
parser.add_argument("-c", "--config", required=True, help="path to the AmpTools config template file")
if len(sys.argv) == 1: # if the user doesn't supply any arguments, print the help string and exit
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

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

gather(bin_directory, config_template)
