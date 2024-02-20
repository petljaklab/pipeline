"""dump_flatfiles.py
Dumps database flatfiles to petljak lab directories
The idea is that whenever someone is faced with a wall of MP[A-Z][0-9]+ IDs, they can cat the flatfile to check what everything is
"""
from petljakapi import connection, q
import petljakapi.select
import petljakapi.translate
import petljakapi.dbs
import pandas as pd
import os

db = "petljakdb_devel"

## Start with studies
studies_dump = petljakapi.select.multi_select(db, "studies", filters = None)
studies_dump = pd.DataFrame(studies_dump, columns = ["id", "rname", "study_pmid", "ncbi_bioproject_id"])
studies_dump["study_id"] = studies_dump["id"].apply(petljakapi.translate.idtostring, prefix = "MPP", by_row = "compat")
with open("/gpfs/data/petljaklab/data/studies/flatfile.txt", "w") as f:
    studies_dump.to_csv(f, sep = "\t", na_rep = "NA", index = False)


## Iterate over the studies now
for study_index, study_row in studies_dump.iterrows():
    study_id = study_row["id"]
    study_str_id = study_row["study_id"]
    path = os.path.join(f"/gpfs/data/petljaklab/data/studies/", study_str_id)
    #### Future: Check for analyses
    # Study-level analyses flatfile code here. For now there are no study-level analyses, so we move on
    ####
    samp_tbl = pd.DataFrame(petljakapi.select.select_join_2(db=db, tbl1 = "samples", tbl1_fkey="cell_id", tbl2 = "cells", tbl2_fkey = "id", tbl2_cols = "rname", filters = {"study_id":study_id}, headers = True))
    samp_tbl = samp_tbl.rename(columns=samp_tbl.iloc[0]).drop(samp_tbl.index[0]).reset_index(drop=True)
    samp_tbl["sample_id"] = samp_tbl["id"].apply(petljakapi.translate.idtostring, prefix = "MPS", by_row = "compat")
    with open(os.path.join(path, "samples", "flatfile.txt"), "w") as f:
        samp_tbl.to_csv(f, sep = "\t", na_rep = "NA", index = False)
    ## Now recurse through the samples/sample analyses
    for sample_index, sample_row in samp_tbl.iterrows():
        sample_str_id = sample_row["sample_id"]
        sample_id = sample_row["id"]
        #print(sample_row)
        samp_path = os.path.join(path, "samples", sample_str_id)
        ## Analyses
        analyses_tbl = pd.DataFrame(petljakapi.select.multi_select(db, "analyses", {"samples_id":sample_id, "input_table":"samples"}, headers = True))
        analyses_tbl = analyses_tbl.rename(columns=analyses_tbl.iloc[0]).drop(analyses_tbl.index[0]).reset_index(drop=True)
        analyses_tbl["analysis_id"] = analyses_tbl["id"].apply(petljakapi.translate.idtostring, prefix = "MPA", by_row = "compat")
        if not analyses_tbl.empty:
            os.makedirs(os.path.join(samp_path, "analyses"), exist_ok = True)
            with open(os.path.join(samp_path, "analyses", "flatfile.txt"), "w") as f:
                analyses_tbl.to_csv(f, sep = "\t", na_rep = "NA", index = False)
        ## Now get runs table
        run_tbl = pd.DataFrame(petljakapi.select.multi_select(db=db, table = "runs", filters = {"sample_id":sample_id}, headers = True))
        run_tbl = run_tbl.rename(columns=run_tbl.iloc[0]).drop(run_tbl.index[0]).reset_index(drop=True)
        run_tbl["run_id"] = run_tbl["id"].apply(petljakapi.translate.idtostring, prefix = "MPR", by_row = "compat") 
        os.makedirs(os.path.join(samp_path, "runs"), exist_ok=True)
        with open(os.path.join(samp_path, "runs", "flatfile.txt"), "w") as f:
            run_tbl.to_csv(f, sep = "\t", na_rep = "NA", index = False)
        for run_index, run_row in run_tbl.iterrows():
            run_str_id = run_row["run_id"]
            run_id = run_row["id"]
            run_path = os.path.join(samp_path, "runs", run_str_id)
            analyses_tbl = pd.DataFrame(petljakapi.select.multi_select(db, "analyses", {"runs_id":run_id, "input_table":"runs"}, headers = True))
            analyses_tbl = analyses_tbl.rename(columns=analyses_tbl.iloc[0]).drop(analyses_tbl.index[0]).reset_index(drop=True)
            analyses_tbl["analysis_id"] = analyses_tbl["id"].apply(petljakapi.translate.idtostring, prefix = "MPA", by_row = "compat")
            if not analyses_tbl.empty:
                os.makedirs(os.path.join(run_path, "analyses"), exist_ok = True)
                with open(os.path.join(run_path, "analyses", "flatfile.txt"), "w") as f:
                    analyses_tbl.to_csv(f, sep = "\t", na_rep = "NA", index = False)
            


    
