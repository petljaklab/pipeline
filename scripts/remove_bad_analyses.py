import petljakapi
import petljakapi.select
import petljakapi.inserts
import petljakapi.translate

import os
import glob

analysis_paths = glob.glob("/gpfs/data/petljaklab/data/studies/*/samples/*/analyses/*/M*")

for p in analysis_paths:
    terminal = petljakapi.translate.stringtoid(os.path.basename(p))
    ret = petljakapi.select.simple_select("petljakdb", "analyses", "id", terminal)
    if len(ret) == 0:
        print(p)