import petljakapi
import petljakapi.select
import petljakapi.inserts
import petljakapi.translate

import argparse

parser = argparse.ArgumentParser(
    prog="add_element.py",
    description="Database wrapper script for inserting an entry into the petljakdb"
)
parser.add_argument('-t', '--type', type = str, help = "What classification to input - study, sample, or run?")
parser.add_argument('-r', '--rname', type = str, help = "rname of what we're inserting")
parser.add_argument('-s', '--source', type = str)
parser.add_argument('-d', '--db', type = str)
parser.add_argument('--study', type = str, help = "study id", required= False)
parser.add_argument('--sample', type = str, help = "sample id", required = False)

args = parser.parse_args()

if args.type == "study":
    petljakapi.inserts.generic_insert({"rname":args.rname}, "studies", args.db)
elif args.type == "sample":
    petljakapi.inserts.generic_insert({"rname":args.rname, "study_id":petljakapi.translate.stringtoid(args.study)}, "samples", args.db)
elif args.type == "run":
    petljakapi.inserts.generic_insert({"rname":args.rname, "study_id":petljakapi.translate.stringtoid(args.study), "sample_id":petljakapi.translate.stringtoid(args.sample), "source":args.source}, "runs", args.db)
