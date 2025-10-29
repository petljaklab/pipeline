"""
validate_files.py [last_run]

Validates the presence of petljaklab analyses' output files, where those analyses are marked "complete".
Prints a table of files, containing their timestamps and names, and whether they are a raw or summary file.

If last_run file is provided, checks for:
    1. new files that did not exist before
    2. changed timestamps

The differences are exported to STDERR for logging/email notifications
"""
import os
import datetime
import argparse
import sys

import petljakapi
import petljakapi.select
import petljakapi.translate

from lib.input_functions import gateway


parser = argparse.ArgumentParser(
    prog="validate_files.py",
    description="Checks timestamps of all analyses marked as complete in the Petljakdb. Optionally checks if any are new or modified"
)

parser.add_argument("-i", "--input", help = "previous_run", required = False)

args = parser.parse_args()

analyses = petljakapi.select.multi_select("petljakdb", "analyses", filters = {"analysis_complete":"True"}, headers = True)

anames = analyses[0]
analyses = analyses[1:]
scratch_dir = "/gpfs/data/petljaklab/data/"
prod_dir = scratch_dir

analyses_dict = {}
old_analyses_dict = {}
if args.input:
    with open(args.input, "r") as f:
        old_table = [line.strip().split("\t") for line in f]
    [old_analyses_dict.update({line[0]:line[1]}) for line in old_table]

for ana in analyses:
    critical_files = gateway(ana[anames.index("analysis_type")], petljakapi.translate.idtostring(ana[anames.index("samples_id")], "MPS"), scratch_dir, prod_dir, "petljakdb", bench = False, ref = ana[anames.index("reference_genome")], dryrun = True, quiet = True)
    ## ignore this annoying file path
    critical_files = [crit for crit in critical_files if not crit.endswith("line_of_origin.txt")]
    critical_file_times = [datetime.datetime.fromtimestamp(os.path.getmtime(f)).strftime("%Y-%m-%d %H:%M:%S") for f in critical_files]
    [analyses_dict.update({file:stamp}) for file,stamp in zip(critical_files, critical_file_times)]
    [print(f"{file}\t{stamp}") for file, stamp in zip(critical_files, critical_file_times)]

if args.input:
    for file, stamp in analyses_dict.items():
        if old_analyses_dict.get(file):
            if stamp != old_analyses_dict[file]:
                print(f"Warning! File {file} with old timestamp {old_analyses_dict[file]} has been modified since yesterday, at {stamp}", file = sys.stderr)

