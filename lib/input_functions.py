import petljakapi
import petljakapi.select
import petljakapi.inserts
import petljakapi.translate
from modules.db_deps import db_deps, module_inputs, module_outputs
def gateway(analysis_name, given_id, scratch_dir, prod_dir, db = "petljakdb") -> list:
    """
    Creates a database entry for the requested analysis and returns the directory that will be used for that entry. 

    Parameters:
    analysis_name (str): Name of the analysis, e.g. "FASTQ" or "WGS_MERGE_BAM".
    given_id: Human-readable lowest-level ID. ie. run ID formatted as MPRXXXXXX, where Xs = Numbers. If no run, provide MPS ID. If no MPS, provide MPP.
    scratch_dir: Directory for scratch space for the pipeline. Should be specified in the config for pipeline usage

    Returns:
    list: List of paths associated with the output of the pipeline to be run. E.g. running with the FASTQ analysis_name returns a path to the .fq files to be created by the pipeline
    """
    ## initialize empty dict for analysis entry
    ## We do this because we add to the dict according to what we're making
    analysis_entry = {}
    #print(given_id)
    
    id_dict = petljakapi.select.parent_ids(in_id = given_id, db = db)
    #print(id_dict)
    ## Determine what the analysis needs to find in the dict
    deps = db_deps[analysis_name]
    if not all(elem in id_dict.keys() for elem in deps):
        raise ValueError(f"Requested analysis {analysis_name} requires db entries {deps} but only provided {id_dict.keys()}")
    ## Identify the database level we need for the input to this pipeline
    ## Should be the final one (ie. should ensure that the db_deps is always ordered hierarchically)
    terminal_dep = deps[-1]
    ## Handle fastq inputs
    ## Maybe to be replaced with a function for better readability
    if analysis_name == "FASTQ":
        ## Get the entry of the "run" table corresponding to the given ID
        entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]))
        ## index 4 is the "source" column
        source = entry[0][4]
        table_id = entry[0][0]
        end_path = ["fq/reads1.fq", "fq/reads2.fq"]
        path_prefix = scratch_dir
        if source == "SRA":
            INPIPE = "SRA"
        elif source == "EGA":
            INPIPE = "EGA"
        elif source == "synthetic_test":
            INPIPE = "TESTS"
        elif source == "external_bam":
            INPIPE = "EXTERNAL_BAM"
        else:
            raise ValueError(f"Handling of fastq source {source} not yet implemented")
    ## Handle Merge inputs
    elif analysis_name == "WGS_MERGE_BAM":
        ## Select the entry of the relevant sample table
        entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]))
        ## Index 0 is the id
        INPIPE = "GATK_BAM"
        table_id = entry[0][0]
        end_path = ["merge/hg19/merged.cram", "merge/hg19/merged.cram.crai", "merge/hg19/merged.done"]
        path_prefix = prod_dir
    elif analysis_name == "MUTECT":
        ## Select the entry of the relevant sample table
        entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]))
        table_id = entry[0][0]
        cell_id = entry[0][6]
        parent_id = entry[0][5]
        if cell_id is not None and cell_id != "NULL":
            INPIPE = "MUTECT_CELLLINE"
        else:
            raise ValueError(f"Handling of non-cell line mutect pipeline runs not yet supported")
        if parent_id is not None and parent_id != "NULL":
            end_path = ["mutect2/hg19/germ/filtered.vcf",
                        "mutect2/hg19/std/filtered.vcf"]
        else:
            end_path = ["mutect2/hg19/parental/table_raw.txt"]
        #print(end_path)
        path_prefix = prod_dir

    elif analysis_name == "LOAD_EXTERNAL_BAM":
        entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]))
        table_id = entry[0][0]
        end_path = ["db_load_runs/runs.loaded"]
        path_prefix = scratch_dir
        INPIPE = "LOAD_EXT_BAM"
    elif analysis_name == "EXTERNAL_BAM":
        entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]))
        table_id = entry[0][0]
        end_path = ["splitfq/split.done"]
        path_prefix = scratch_dir
        INPIPE = "EXTERNAL_BAM"
    else:
        raise ValueError(f"Handling of {analysis_name} not yet supported")
    ## Now we construct the input file to make
    path = f"{path_prefix}studies/{id_dict['studies']}/"
    analysis_entry.update({"studies_id":petljakapi.translate.stringtoid(id_dict['studies'])})
    if id_dict["samples"]:
        path = path + f"samples/{id_dict['samples']}/"
        analysis_entry.update({"samples_id":petljakapi.translate.stringtoid(id_dict['samples'])})
    if id_dict["runs"]:
        path = path + f"runs/{id_dict['runs']}/"
        analysis_entry.update({"runs_id":petljakapi.translate.stringtoid(id_dict['runs'])})
    # Now we insert the analysis entry into the database, and use that ID to construct the path
    path = path + "analyses/" + INPIPE + "/"
    terminal_dep_string = terminal_dep + "_id"
    analysis_entry.update({
        "pipeline_name":INPIPE,
        "pipeline_version":"1.0.0",
        "analysis_type":analysis_name,
        "analysis_dir":path,
        "input_table":terminal_dep,
        terminal_dep_string:table_id
    })
    ## Insert/get the relevant ID
    analysis_id = petljakapi.inserts.analysis_insert(analysis_entry, "analyses", db = db)[0][0]
    ## Construct final IDs
    out_paths = [path + petljakapi.translate.idtostring(analysis_id, "MPA") + "/" + p for p in end_path]
    return(out_paths)
