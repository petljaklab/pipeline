import glob
import petljakapi
import petljakapi.select
import petljakapi.inserts
import petljakapi.translate
from modules.db_deps import db_deps, module_inputs, module_outputs
#def get_seq_strategies(sample_id, db = "petljakdb") -> dict:
    

def gateway(analysis_name, given_id, scratch_dir, prod_dir, db = "petljakdb", bench = False, ref = "hg19", dryrun = False, quiet = False) -> list:
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
    direct = False
    ## Default ref value because some analyses don't have a ref
    analysis_ref = None
    ## Maybe to be replaced with a function for better readability
    if analysis_name == "SOMATIC":
        outlist = []
        for pipe in ["SBS", "INDEL"]:
            outlist.extend(gateway(analysis_name=pipe, given_id=given_id, scratch_dir=scratch_dir, prod_dir=prod_dir, db=db, ref=ref))
        return(outlist)
    ## Handle fastq inputs
    if analysis_name == "FASTQ":
        ## Get the entry of the "run" table corresponding to the given ID
        if terminal_dep:
            entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]), bench = bench)
        else:
            entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict), bench = bench)
        ## index 4 is the "source" column
        source = entry[0][4]
        table_id = entry[0][0]
        samp = entry[0][2]
        ## Check if we already have a processed CRAM to rip the reads from
        existing_done = petljakapi.select.multi_select(db = db, table = "analyses", filters = {"pipeline_name":"GATK_BAM", "analysis_complete":"True", "samples_id":samp, "reference_genome":ref})
        existing = petljakapi.select.multi_select(db = db, table = "analyses", filters = {"pipeline_name":"GATK_BAM", "analysis_complete":"True", "samples_id":samp})
        end_path = ["fq/reads1.fq", "fq/reads2.fq"]
        path_prefix = scratch_dir
        if source == "local":
            INPIPE = "LOCAL"
            direct = True
            ## entry[9] corresponds to fastq path
            #print(entry[0][9])
            r1 = glob.glob(entry[0][9] + "1*")[0]
            r2 = glob.glob(entry[0][9] + "2*")[0]
            #print([r1, r2])
            return([r1, r2])
        elif existing and not existing_done:
            INPIPE = "GATK_BAM_CONVERT"
        elif source == "SRA":
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
        entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]), bench = bench)
        ## Need to decide which pipeline to use - GATK best practices, or UDSeq (currently)
        ## The way we achieve this is by accessing all the runs, 1. ensuring they're all homogenous in their sequencing strategy, and 2. using that value
        runs_entries = petljakapi.select.multi_select(db, "runs", {"sample_id":entry[0][0]})
        strategies = [element[8] for element in runs_entries]
        strategies = list(set([strat for strat in strategies if strat in ["WGS", "WXS", "UDSEQ"]]))[0]
        ## TODO: If there are some WGS and some UDseq libraries for the same sample, then this won't work. Need to find a solution to loop over all ways to satisfy this endpoint and make them all for their relevant runs
        path_prefix = prod_dir
        analysis_ref=ref
        if strategies == "UDSEQ":
            ## do UDseq things
            ## Enforce hg38
            if ref != "hg38":
                raise ValueError("UDSeq requires hg38!")
            INPIPE = "UDSEQ_BAM"
            table_id = entry[0][0]
            if entry[0][5] is not None:
                end_path = [f"merge/{ref}/calls_paired/calls_snv_private.txt", f"merge/{ref}/calls_paired/calls_coverage_stats.txt"]
            else:
                end_path = [f"merge/{ref}/calls/calls_snv_private.txt", f"merge/{ref}/calls/calls_coverage_stats.txt"]
        elif strategies in ["WGS", "WXS"]:
            ## Index 0 is the id
            INPIPE = "GATK_BAM"
            table_id = entry[0][0]
            end_path = [f"merge/{ref}/merged.cram", f"merge/{ref}/merged.cram.crai", f"merge/{ref}/merged.done"]
    elif analysis_name == "SBS":
        ## Select the entry of the relevant sample table
        entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]), bench = bench)
        table_id = entry[0][0]
        cell_id = entry[0][6]
        parent_id = entry[0][5]
        path_prefix = prod_dir
        analysis_ref=ref
        if cell_id is not None and cell_id != "NULL":
            INPIPE = "MUTECT_CELLLINE"
            ## Test if this sample is used as a reference for any other samples
            daughters = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "sample_parent_id", filter_value = table_id, bench = bench)
            ## Init list of endpoints, beginning with line of origin (later will need to have another if-else block for biopsy vs cell lines to exclude this)
            end_path = [f"mutect2/{ref}/proc/line_of_origin.txt"]
            ## First, if we have a parent (ie this is a daughter line), need to add the correct endpoints for daughters
            #print(entry)
            #print(parent_id)
            if parent_id is not None and parent_id != "NULL":
                end_path.extend([f"mutect2/{ref}/proc/variants_final.done"])
            if len(daughters) > 0: ## ie. if this is a parent
                end_path.extend([f"mutect2/{ref}/parental/table_raw.txt"])
        else:
            INPIPE = "MUTECT_BIOP"
            pat_id = entry[0][8]
            germline = petljakapi.select.multi_select(db, "patients", {"id":pat_id})[0]
            germline = germline[2]
            if germline == table_id:
                if not quiet:
                    print(f"Calling germline variants is not yet supported. Skipping sample {id_dict[terminal_dep]}")
                return()
            end_path = [f"mutect2/{ref}/biop/filtered.vcf"]
            if not germline:
                end_path = [f"mutect2/{ref}/biop_noref/filtered.vcf"]
            ## If this is the germline sample, simply pass (we don't do germline calling yet)
            

    elif analysis_name == "INDEL":
        entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]), bench = bench)
        table_id = entry[0][0]
        parent_id = entry[0][5]
        analysis_ref=ref
        if parent_id is None:
            if not quiet:
                print(f"Calling indels without matched normal is not supported. Skipping sample {id_dict[terminal_dep]}")
            return()
        callers = ["mutect2", "strelka2", "varscan2"]
        callers = ["mutect2"]
        end_path = [f"{ref}/{caller}/indels.txt" for caller in callers]
        end_path.extend([f"{ref}/{caller}/indels.vcf.dbtmp" for caller in callers])
        path_prefix = prod_dir
        INPIPE = "INDEL"
    elif analysis_name == "LOAD_EXTERNAL_BAM":
        entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]), bench = bench)
        table_id = entry[0][0]
        end_path = ["db_load_runs/runs.loaded"]
        path_prefix = scratch_dir
        INPIPE = "LOAD_EXT_BAM"
    elif analysis_name == "EXTERNAL_BAM":
        entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]), bench = bench)
        table_id = entry[0][0]
        end_path = ["splitfq/split.done"]
        path_prefix = scratch_dir
        INPIPE = "EXTERNAL_BAM"
    elif analysis_name == "SV":
        entry = petljakapi.select.simple_select(db = db, table = terminal_dep, filter_column = "id", filter_value = petljakapi.translate.stringtoid(id_dict[terminal_dep]), bench = bench)
        table_id = entry[0][0]
        analysis_ref = ref
        parent_id = entry[0][5]
        if parent_id is None:
            #print(f"Calling SVs without matched normal is not supported. Skipping sample {id_dict[terminal_dep]}")
            end_path = [f"{ref}/gridss_singlesamp/somatic.filtered.gnomAD.sv.rds"]
        else:
            end_path = [f"{ref}/gridss/high_low_confidence_annotated.txt"]
        path_prefix = scratch_dir
        INPIPE = "SV"
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
        terminal_dep_string:table_id,
        'reference_genome':analysis_ref
    })
    ## Insert/get the relevant ID
    analysis_id = petljakapi.inserts.analysis_insert(analysis_entry, "analyses", db = db, dryrun = dryrun)[0][0]
    ## Construct final IDs
    out_paths = [path + petljakapi.translate.idtostring(analysis_id, "MPA") + "/" + p for p in end_path]
    return(out_paths)
