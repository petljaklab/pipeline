import os
import subprocess
import argparse
import snakemake


import petljakapi
from petljakapi import connection
import petljakapi.select
import petljakapi.inserts
import petljakapi.translate
from modules.db_deps import db_deps, module_inputs, module_outputs

## Hardcoded for now
db = "petljakdb_devel"

parser = argparse.ArgumentParser(
    prog="executor.py",
    description="Executes Snakemake and handles Database & API calls as needed"
)
parser.add_argument("-i", '--id', help = "String IDs", nargs = "+", action = "extend", required = True)
parser.add_argument("-p", "--pipeline", help = "Pipeline name to be run. Dependency pipelines will be run as needed.")
parser.add_argument("-n", '--dryrun', action=argparse.BooleanOptionalAction)

args = parser.parse_args()
available_pipelines = ["TESTS", "SRA", "GATK_BAM"]

prefixes = ["MPP", "MPS", "MPR", "MPA"]

## Check inputs
if args.pipeline not in available_pipelines:
    raise ValueError(f"Requested pipeline not in {available_pipelines}")

id = args.id

## Sanity check the IDs
for i in id:
    if not i.startswith(tuple(prefixes)):
        raise ValueError(f"Supplied ID does not start with one of {', '.join(prefixes)}")
    
id_dict = petljakapi.select.parent_ids(id, db = db)


## Create db entry of the incoming analysis
## Try doing this as a function I guess? Maybe can repurpose it for the input function (and then put it into a shared library for the pipeline to import)
def gateway(analysis_name, id_dict):
    ## Determine what the analysis needs to find in the dict
    deps = db_deps[analysis_name]
    if not all(elem in id_dict.keys() for elem in deps):
        raise ValueError(f"Requested analysis {analysis_name} requires db entries {deps} but only provided {id_dict.keys()}")
    ## Determine what the required input pipeline is
    ## First determine the input type we need e.g. FASTQ, WGS_MERGE_BAM
    req_inputs = module_inputs[analysis_name]
    ## Identify the database level we need for the input to this pipeline
    ## Should be the final one (ie. should ensure that the db_deps is always ordered hierarchically)
    terminal_dep = deps[-1]
    ## Handle fastq inputs
    if req_inputs[0] == "FASTQ":
        ## Get the entry of the "run" table corresponding to the given ID
        runs_entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = id_dict[terminal_dep])
        ## index 4 is the "source" column
        source = runs_entry[0][4]
        if source == "SRA":
            INPIPE = "SRA"
        elif source == "synthetic_test":
            INPIPE = "TESTS"
        else:
            raise ValueError(f"Handling of fastq source {source} not yet implemented")
    else:
        raise ValueError(f"Handling of {req_inputs[0]} not yet supported")
    
    ## Now we construct the input file to make
    path = f"studies/{id_dict['studies']}/"
    if "samples" in id_dict.keys():
        path = path + f"samples/{id_dict['samples']}/"
    if "runs" in id_dict.keys():
        path = path + f"runs/{id_dict['runs']}/"
    # Now we insert the analysis entry into the database, and use that ID to construct the path
    path = path + "analyses/"
        ## You were here on Dec 29th at about 3pm
    return(True)




## Check that the entries exist
## Remember - here we are not inserting into the PSR tables, you need to do that before you execute this
#study = petljakapi.select.simple_select("petljakdb_devel", "studies", "id", study)

def pipeline_init(terminal_id, analysis_name, id_type = None, ):
    ## Takes a given ID and analysis name and does a few things
    ## First thing we need to do is get an analysis ID, either from existing matching analysis or make a new db entry with pending status
    ## Return the relevant path with AIDs filled in
    ## e.g. file path: "{arbitrary path}/studies/MPP000001/samples/MPS000001/runs/MPS000001/analyses/TESTS/MPA000001/flatfile.txt"
    return(True)

snakemake.snakemake(snakefile="./Snakefile", dryrun=True)