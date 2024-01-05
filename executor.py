import os
import subprocess
import argparse
import snakemake
import yaml
import sys


import petljakapi
from petljakapi import connection
import petljakapi.select
import petljakapi.inserts
import petljakapi.translate
from modules.db_deps import db_deps, module_inputs, module_outputs
import lib.input_functions

## Hardcoded for now
db = "petljakdb_devel"

parser = argparse.ArgumentParser(
    prog="executor.py",
    description="Executes Snakemake and handles Database & API calls as needed"
)
group = parser.add_mutually_exclusive_group(required=True)

group.add_argument('--id', help = "String IDs", nargs = "+", action = "extend")
group.add_argument('--idfile', help = "List of IDs")
parser.add_argument("--pipeline", help = "Analysis requested")
parser.add_argument("--cell_id", required=False)

script_args, unknown = parser.parse_known_args()

print(script_args.id)

## Parse IDs
if script_args.id is not None:
    id = script_args.id
elif script_args.idfile is not None:
    with open(script_args.idfile, "r") as f:
        id = f.readlines()
    id = [line.strip() for line in id]


available_pipelines = list(module_outputs.keys())
available_analyses = list(db_deps.keys())

prefixes = ["MPP", "MPS", "MPR", "MPA"]

## Config
configfile = os.path.join(os.path.dirname("__file__"), "config.yaml")
with open(configfile) as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

## Check inputs
if script_args.pipeline not in available_analyses:
    raise ValueError(f"Requested pipeline not in {available_analyses}")

if not unknown:
    unknown = []
## Delete argparse arguments?
sys.argv = [sys.argv[0]]
sys.argv.extend(unknown)
#print(sys.argv)


## Sanity check the IDs
for i in id:
    if not i.startswith(tuple(prefixes)):
        raise ValueError(f"Supplied ID does not start with one of {', '.join(prefixes)}")
    analysis_path = lib.input_functions.gateway(analysis_name=script_args.pipeline, given_id=i, scratch_dir = config["SCRATCH_DIR"], prod_dir = config["PROD_DIR"], db = config["db"])
    unknown.extend(analysis_path)




## This is going to be the TESTS tp53 test for now
# First get the analysis ID for what we're trying to make so we can construct it

#analysis_path = lib.input_functions.gateway(analysis_name=script_args.pipeline, given_id=id[0], scratch_dir = config["SCRATCH_DIR"], prod_dir = config["PROD_DIR"], db = config["db"])

print(unknown)

snakemake.main(argv = unknown)
#snakemake.snakemake(snakefile="./Snakefile", dryrun = True, targets = analysis_path)
