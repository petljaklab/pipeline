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
parser.add_argument('--id', help = "String IDs", nargs = "+", action = "extend", required = True)
parser.add_argument("--pipeline", help = "Analysis requested")

script_args, unknown = parser.parse_known_args()

available_pipelines = ["TESTS", "SRA", "GATK_BAM"]
available_analyses = ["FASTQ", "WGS_MERGE_BAM"]

prefixes = ["MPP", "MPS", "MPR", "MPA"]

## Config
configfile = os.path.join(os.path.dirname("__file__"), "config.yaml")
with open(configfile) as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

## Check inputs
if script_args.pipeline not in available_pipelines:
    raise ValueError(f"Requested pipeline not in {available_pipelines}")

id = script_args.id

## Sanity check the IDs
for i in id:
    if not i.startswith(tuple(prefixes)):
        raise ValueError(f"Supplied ID does not start with one of {', '.join(prefixes)}")

id_dict = petljakapi.select.parent_ids(in_id = id[0], db = db)


## Check that the entries exist
## Remember - here we are not inserting into the PSR tables, you need to do that before you execute this
#study = petljakapi.select.simple_select("petljakdb_devel", "studies", "id", study)

#def pipeline_init(terminal_id, analysis_name, id_type = None, ):
    ## Takes a given ID and analysis name and does a few things
    ## First thing we need to do is get an analysis ID, either from existing matching analysis or make a new db entry with pending status
    ## Return the relevant path with AIDs filled in
    ## e.g. file path: "{arbitrary path}/studies/MPP000001/samples/MPS000001/runs/MPS000001/analyses/TESTS/MPA000001/flatfile.txt"
#    return(True)

## Run the initial input function
## This is going to be the TESTS tp53 test for now
# First get the analysis ID for what we're trying to make so we can construct it
analysis_path = lib.input_functions.gateway(analysis_name=script_args.pipeline, id_dict = id_dict, scratch_dir = config["SCRATCH_DIR"], prod_dir = config["PROD_DIR"])

## Delete argparse arguments?
sys.argv = [sys.argv[0]]
sys.argv.extend(unknown)
print(sys.argv)
if not unknown:
    unknown = []
unknown.extend(analysis_path)
snakemake.main(argv = unknown)
#snakemake.snakemake(snakefile="./Snakefile", dryrun = True, targets = analysis_path)
