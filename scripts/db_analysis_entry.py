import petljakapi
import petljakapi.select
import petljakapi.inserts

import argparse

parser = argparse.ArgumentParser(
    prog="db_analysis_entry.py",
    description="Database wrapper script for inserting an analysis entry into the petljakdb"
)
parser.add_argument('-t', '--type', type = str, help = "What classification of input - study, sample, or run?")
parser.add_argument('-p', '--pipeline', type = str, help = "Pipeline name")
parser.add_argument('-e', '--version', type = str, help = "Version of the pipeline")
parser.add_argument('-d', '--dir', type = str, help = "analysis directory")
parser.add_argument('-c', '--container', type = str)
parser.add_argument('-p', '--path', type = str)
parser.add_argument('-o', '--outfile', type = str)

args = parser.parse_args()

## First check if we added this entry to studies. The insert should run basically once ever.
proj_id = petljakapi.select.simple_select(args.db, "studies", "rname", args.study)
if proj_id:
    proj_id = proj_id[0][0]
else:
    proj_id = petljakapi.inserts.generic_insert({"rname":args.study}, "studies")[0][0]
## Now check samples
sample_id = petljakapi.select.simple_select(args.db, "samples", "rname", args.sample)
if sample_id:
    sample_id = sample_id[0][0]
else:
    sample_id = petljakapi.inserts.generic_insert({"rname":args.sample, "study_id":proj_id}, "samples", args.db)[0][0]
run_id = petljakapi.select.simple_select(args.db, "runs", "rname", args.run)
if run_id:
    run_id = run_id[0][0]
else:
    run_id = petljakapi.inserts.generic_insert({"rname":args.run, "study_id":proj_id, "sample_id":sample_id, "source":"local", "local_path":args.path}, "runs", args.db)[0][0]

with open(args.outfile, "w") as f:
    pass
