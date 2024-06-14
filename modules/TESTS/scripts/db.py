import petljakapi
import petljakapi.select
import petljakapi.inserts

import argparse

parser = argparse.ArgumentParser(
    prog="db.py",
    description="Database wrapper script for pipeljak.tests pipeline module"
)
parser.add_argument('-s', '--study', type = str)
parser.add_argument('-m', '--sample', type = str)
parser.add_argument('-r', '--run', type = str)
parser.add_argument('-d', '--db', type = str)
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
