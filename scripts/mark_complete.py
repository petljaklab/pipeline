import argparse
import petljakapi
from petljakapi import connection
import petljakapi.select
import petljakapi.inserts
import petljakapi.translate
import petljakapi.update
from modules.db_deps import db_deps, module_inputs, module_outputs
import os

parser = argparse.ArgumentParser(
    prog="mark_complete.py",
    description="Marks an analysis as complete, updates the output files entry and updates the timestamp in the petljakdb"
)

parser.add_argument("-i", "--id", help = "analysis ID")
parser.add_argument("-d", "--db", help = "database name")
parser.add_argument("outputs", nargs = "+", help = "List of output files")

args = parser.parse_args()

analysis_id = petljakapi.translate.stringtoid(args.id)

from datetime import datetime
now = datetime.now()
timestamp = now.strftime("%d/%m/%Y %H:%M:%S")
petljakapi.update.update(args.db, "analyses", {"id":analysis_id}, "analysis_time", timestamp)
petljakapi.update.update(args.db, "analyses", {"id":analysis_id}, "analysis_complete", "True")
petljakapi.update.update(args.db, "analyses", {"id":analysis_id}, "analysis_dir", os.path.dirname(args.outputs[0]))


