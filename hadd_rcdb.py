#!/usr/bin/env python3

##########################################################################################################################

import glob
import os
import os.path
import re
import subprocess
import sys
from optparse import OptionParser

#################################################### RCDB ENVIRONMENT ####################################################
#os.environ["RCDB_HOME"] = "/group/halld/Software/builds/Linux_CentOS7.7-x86_64-gcc4.8.5/rcdb/rcdb_0.06.00/"
#sys.path.append("/group/halld/Software/builds/Linux_CentOS7.7-x86_64-gcc4.8.5/rcdb/rcdb_0.06.00/python")
import rcdb

db = rcdb.RCDBProvider("mysql://rcdb@hallddb.jlab.org/rcdb")

#################################################### GLOBAL VARIABLES ####################################################

VERBOSE = True

####################################################### FIND FILES #######################################################


def find_files(SIGNATURE):
    """

    :param SIGNATURE: 

    """

    # SEARCH FOR THE FILES
    file_signature = SIGNATURE
    file_list = glob.glob(file_signature)
    if (VERBOSE == True):
        print("size of file_list is " + str(len(file_list)))

    return file_list


def main(argv):
    """

    :param argv: 

    """
    parser_usage = "hadd.py minrun maxrun directory [basename]"
    parser = OptionParser(usage=parser_usage)
    (options, args) = parser.parse_args(argv)

    if (len(args) < 3):
        parser.print_help()
        return

    BASENAME = "tree"  #Default is to look for all tree*.root files

    # GET ARGUMENTS
    MINRUN = int(args[0])
    MAXRUN = int(args[1])
    DIR = args[2]
    if len(args) == 4:
        BASENAME = args[3]

    print("Getting lists of good runs to choose from")
    # GET LISTS OF GOOD RUNS TO CHOOSE FROM
    RCDB_RUNS = db.select_runs("@is_production and @status_approved",
                               run_min=MINRUN,
                               run_max=MAXRUN)
    if MINRUN > 40000 and MINRUN < 60000:
        RCDB_RUNS = db.select_runs("@is_2018production and @status_approved",
                                   run_min=MINRUN,
                                   run_max=MAXRUN)
    if MINRUN > 70000:
        RCDB_RUNS = db.select_runs("@is_dirc_production and @status_approved",
                                   run_min=MINRUN,
                                   run_max=MAXRUN)
    RCDB_RUN_NUMBERS = [run.number for run in RCDB_RUNS]

    # GET LIST OF FILES ON DISK
    os.chdir(DIR)
    print(f"Directory: {DIR}")
    FOUND_HIST_FILES = find_files("%s*.root" % BASENAME)

    HADD_RUNS_TOT = "hadd tree_sum_%d_%d.root tree_sum_*" % (MINRUN, MAXRUN)
    HADD_RUNS_AMO = "hadd tree_sum_AMO_%d_%d.root" % (MINRUN, MAXRUN)
    HADD_RUNS_PARA = "hadd tree_sum_PARA_0_%d_%d.root" % (MINRUN, MAXRUN)
    HADD_RUNS_PERP = "hadd tree_sum_PERP_90_%d_%d.root" % (MINRUN, MAXRUN)
    HADD_RUNS_PARA_45_135 = "hadd tree_sum_PARA_135_%d_%d.root" % (MINRUN,
                                                                   MAXRUN)
    HADD_RUNS_PERP_45_135 = "hadd tree_sum_PERP_45_%d_%d.root" % (MINRUN,
                                                                  MAXRUN)

    AMO = 0
    PARA = 0
    PERP = 0
    PARA_45_135 = 0
    PERP_45_135 = 0

    MISSING_RUNS = []

    # FIND/ADD JOBS
    for RCDB_RUN in RCDB_RUNS:

        RUN = RCDB_RUN.number

        # Check RCDB status for each run number
        if RUN not in RCDB_RUN_NUMBERS and RUN != 10000:
            continue

        # Format run and file numbers
        FORMATTED_RUN = "%d" % RUN
        #print FORMATTED_RUN

        # Check if expected file is on disk
        if not any(FORMATTED_RUN in x for x in FOUND_HIST_FILES):
            #print FORMATTED_RUN
            MISSING_RUNS.append("0%s" % FORMATTED_RUN)
            continue  # skip files that are missing

        # Histogram file name to create hadd
        HIST_FILE = next(filter(lambda a: FORMATTED_RUN in a, FOUND_HIST_FILES))

        # AMO, PARA, and PERP
        conditions_by_name = RCDB_RUN.get_conditions_by_name()
        conditions = RCDB_RUN.get_conditions_by_name().keys()
        if 'RL' in str(RCDB_RUN.get_condition('radiator_type')) or 'Al' in str(
                RCDB_RUN.get_condition('radiator_type')):
            AMO += 1
            HADD_RUNS_AMO += " %s" % HIST_FILE

        # Spring 2017 when polarization_angle is defined
        if RCDB_RUN.get_condition('polarization_angle'):
            if RCDB_RUN.get_condition('polarization_angle').value == 0:
                PARA += 1
                HADD_RUNS_PARA += " %s" % HIST_FILE
            elif RCDB_RUN.get_condition('polarization_angle').value == 90:
                PERP += 1
                HADD_RUNS_PERP += " %s" % HIST_FILE
            elif RCDB_RUN.get_condition('polarization_angle').value == 135:
                PARA_45_135 += 1
                HADD_RUNS_PARA_45_135 += " %s" % HIST_FILE
            elif RCDB_RUN.get_condition('polarization_angle').value == 45:
                PERP_45_135 += 1
                HADD_RUNS_PERP_45_135 += " %s" % HIST_FILE
        else:  # Spring 2016 when only polarization_direction was defined
            if RCDB_RUN.get_condition('polarization_direction').value == "PARA":
                PARA += 1
                HADD_RUNS_PARA += " %s" % HIST_FILE
            elif RCDB_RUN.get_condition(
                    'polarization_direction').value == "PERP":
                PERP += 1
                HADD_RUNS_PERP += " %s" % HIST_FILE

    # HADD SUMMARY RUNS
    if MINRUN != 10000 or MAXRUN != 10000:
        if AMO > 0:
            subprocess.call(HADD_RUNS_AMO, shell=True)
        if PARA > 0:
            subprocess.call(HADD_RUNS_PARA, shell=True)
        if PERP > 0:
            subprocess.call(HADD_RUNS_PERP, shell=True)
        if PARA_45_135 > 0:
            subprocess.call(HADD_RUNS_PARA_45_135, shell=True)
        if PERP_45_135 > 0:
            subprocess.call(HADD_RUNS_PERP_45_135, shell=True)
        subprocess.call(HADD_RUNS_TOT, shell=True)

    # notify about missing runs
    if len(MISSING_RUNS) > 0:
        print()
        print("Missing files for these runs")
        for RUN in MISSING_RUNS:
            print(RUN)


if __name__ == "__main__":
    main(sys.argv[1:])
