import petljakapi
import petljakapi.select
import petljakapi.inserts
import petljakapi.translate
from modules.db_deps import db_deps, module_inputs, module_outputs
## Create db entry of the incoming analysis
## Try doing this as a function I guess? Maybe can repurpose it for the input function (and then put it into a shared library for the pipeline to import)
def gateway(analysis_name, id_dict, scratch_dir, prod_dir, db = "petljakdb_devel"):
    ## Determine what the analysis needs to find in the dict
    deps = db_deps[analysis_name]
    if not all(elem in id_dict.keys() for elem in deps):
        raise ValueError(f"Requested analysis {analysis_name} requires db entries {deps} but only provided {id_dict.keys()}")
    ## Determine what the required input pipeline is
    ## First determine the input type we need e.g. FASTQ, WGS_MERGE_BAM
    req_inputs = module_outputs[analysis_name]
    print(req_inputs)
    ## Identify the database level we need for the input to this pipeline
    ## Should be the final one (ie. should ensure that the db_deps is always ordered hierarchically)
    terminal_dep = deps[-1]
    ## Handle fastq inputs
    if req_inputs == "FASTQ":
        ## Get the entry of the "run" table corresponding to the given ID
        entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]))
        ## index 4 is the "source" column
        source = entry[0][4]
        table_id = entry[0][0]
        end_path = ["fq/reads1.fq", "fq/reads2.fq"]
        path_prefix = scratch_dir
        if source == "SRA":
            INPIPE = "SRA"
        elif source == "synthetic_test":
            INPIPE = "TESTS"
        else:
            raise ValueError(f"Handling of fastq source {source} not yet implemented")
    else:
        raise ValueError(f"Handling of {analysis_name} not yet supported")
    ## Check for an entry already
    petljakapi.select.multi_select(db = db, table = terminal_dep, filters = {"id":source})
    ## Now we construct the input file to make
    path = f"{path_prefix}studies/{id_dict['studies']}/"
    if "samples" in id_dict.keys():
        path = path + f"samples/{id_dict['samples']}/"
    if "runs" in id_dict.keys():
        path = path + f"runs/{id_dict['runs']}/"
    # Now we insert the analysis entry into the database, and use that ID to construct the path
    path = path + "analyses/" + INPIPE + "/"
    terminal_dep_string = terminal_dep + "_id"
    analysis_entry = {
        "pipeline_name":INPIPE,
        "pipeline_version":"1.0.0",
        "analysis_dir":path,
        "input_table":terminal_dep,
        terminal_dep_string:table_id
    }
    ## Insert/get the relevant ID
    analysis_id = petljakapi.inserts.analysis_insert(analysis_entry, "analyses")[0][0]
    ## Construct final IDs
    out_paths = [path + petljakapi.translate.idtostring(analysis_id, "MPA") + "/" + p for p in end_path]
    return(out_paths)
