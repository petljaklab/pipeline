import os
import numpy as np
import glob
import petljakapi.cellline

MUTECT_CELLLINE_PIPELINE_VERSION = "1.0.0"
MUTECT_CELLLINE_PATH = os.path.join(basedir, "modules", "MUTECT_CELLLINE")

#configfile: os.path.join(basedir, "modules", "MUTECT_CELLLINE", "config.yaml")
#locals().update(config)

include: os.path.join(MUTECT_CELLLINE_PATH, "MUTECT2.smk")
include: os.path.join(MUTECT_CELLLINE_PATH, "PROCESS_FILTER_M2.smk")
include: os.path.join(MUTECT_CELLLINE_PATH, "MUTECT2_PARENTAL.smk")

def parent_cell(wildcards):
    sample_id = petljakapi.translate.stringtoid(wildcards.sample)
    db_line = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = sample_id, bench = config["bench"])
    #print(db_line)
    parent_id = petljakapi.translate.idtostring(db_line[0][5], "MPS")
    parent_merge = gateway("WGS_MERGE_BAM", parent_id, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)[0]
    return(parent_merge)

rule MUTECT2_CELLLINE_VCFTOTABLE:
    input:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/filtered_renamed.vcf",
    output:
        tab = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/table_raw.txt"
    resources:
        mem_mb = 3000,
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        runtime = 60 * 4,
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/table.log"
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    shell:
        """
        gatk VariantsToTable -V {input.vcf} \
            -O {output.tab} -raw --split-multi-allelic &> {log}
        """

def get_parents_for_sample(wildcards):
    ## Get parent "comprehensive" VCF paths
    ## 1. API call to get the daughter's parent ID
    ## 2. For that parent ID, get all samples from that study with the same cell line ID and parent ID
    #
    ## Get daughter sample table line
    study_id = petljakapi.translate.stringtoid(wildcards.study)
    daughter_id = petljakapi.translate.stringtoid(wildcards.sample)
    daughter_line = petljakapi.select.multi_select(db = db, table = "samples", filters = {"id":daughter_id}, headers = True, bench = config["bench"])
    ## Get the line for that daughter's annotated parent and get the cell ID from that
    parent_id = np.array(daughter_line[1])[np.where(np.array(daughter_line[0]) == "sample_parent_id")][0]
    cell_id = np.array(daughter_line[1])[np.where(np.array(daughter_line[0]) == "cell_id")][0]
    ## Now get the parent of the parent, so we can get all the possible parents
    parental_sample = petljakapi.select.multi_select(db = db, table = "samples", filters = {"id":parent_id}, headers = True, bench = config["bench"])
    parental_parent = np.array(parental_sample[1])[np.where(np.array(daughter_line[0]) == "sample_parent_id")][0]
    ## Now get all the appropriate parental lines
    all_parental_lines = petljakapi.select.multi_select(db = db, table = "samples", filters = {"study_id":study_id, "sample_parent_id":parental_parent, "cell_id":cell_id}, bench = config["bench"])
    ## Loop through to get the analysis paths and thereby the path to the parents' variant table
    parental_paths = []
    for parent in all_parental_lines:
        ## First gateway to ensure the parent is made
        gateway("MUTECT", petljakapi.translate.idtostring(parent[0], "MPS"), SCRATCH_DIR, config["PROD_DIR"], db = "petljakdb", ref = wildcards.reference)
        ## Get the analysis dir
        analysis = petljakapi.select.multi_select(db = db, table = "analyses", filters = {"samples_id":parent[0], "pipeline_name":"MUTECT_CELLLINE"}, bench = config["bench"])
        analysis = analysis[0]
        vcf_path = os.path.join(analysis[10], petljakapi.translate.idtostring(analysis[0], "MPA"), "mutect2", wildcards.reference, "parental/table_raw.txt")
        parental_paths.append(vcf_path)
    return(parental_paths)



rule MUTECT2_CELLLINE_PARENT_TABLE:
    input:
        get_parents_for_sample
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/parents.table"
    resources:
        mem_mb = 2000,
        slurm_partition = "petljaklab,cpu_dev",
        runtime = 5
    shell:
        """
        for i in {input}; do echo $i >> {output}; done;
        """

def parent_sample_id(wildcards):
    sample_id = petljakapi.translate.stringtoid(wildcards.sample)
    db_line = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = sample_id, bench = config["bench"])
    #print(db_line)
    parent_id = petljakapi.translate.idtostring(db_line[0][5], "MPS")
    return(parent_id)

def full_cell_line_cohort(wildcards):
    study_id = petljakapi.translate.stringtoid(wildcards.study)
    samples_ids = petljakapi.select.multi_select(db, "samples", {"study_id":study_id}, bench = config["bench"])
    o_list = []
    for sample in samples_ids:
        #o_list.extend(gateway("MUTECT", petljakapi.translate.idtostring(sample[0], "MPS"), scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = db))
        ## get the analysis ID for this sample
        analysis = petljakapi.select.multi_select(db, "analyses", {"samples_id":sample[0], "pipeline_name":"MUTECT_CELLLINE", "reference_genome":wildcards.reference})
        runs = petljakapi.select.multi_select(db, "runs", {"sample_id":sample[0]})
        if len(analysis) == 0 or len(runs) == 0:
            print(sample)
            continue
        analysis = analysis[0]
        p = analysis[10]
        aid = petljakapi.translate.idtostring(analysis[0], "MPA")
        ## If this sample has daughters
        daughts = petljakapi.cellline.daughter_cells(sample[0], db = db)
        ## If this sample has any parents
        parents = petljakapi.select.multi_select(db, "samples", {"id":sample[5]}, bench = config["bench"])
        if daughts or not parents:
            o_list.append(os.path.join(p, aid, "mutect2", wildcards.reference, "parental/table_raw.txt"))
        if parents:
            o_list.append(os.path.join(p, aid, "mutect2", wildcards.reference, "germ/table_raw.txt"))
            o_list.append(os.path.join(p, aid, "mutect2", wildcards.reference, "std/table_raw.txt"))
    return(o_list)
    
rule M2_SBS_TABLES:
    #input:
    #    data = full_cell_line_cohort
    output:
        daughter = SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/{reference}/daughters_table.txt",
        parental = SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/{reference}/parents_table.txt",
    threads: 1
    resources:
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        cpus = 1,
        threads = 1,
        mem_mb = 2000,
        runtime = 15,
    run:
        if not os.path.exists(f"{SCRATCH_DIR}studies/{wildcards.study}/analyses/MUTECT_CELLLINE/{wildcards.reference}"):
            os.makedirs(f"{SCRATCH_DIR}studies/{wildcards.study}/analyses/MUTECT_CELLLINE/{wildcards.reference}")
        study_id = petljakapi.translate.stringtoid(wildcards.study)
        samples_ids = petljakapi.select.multi_select(db, "samples", {"study_id":study_id})
        for sample in samples_ids:
            cell_name = petljakapi.select.multi_select(db, "cells", {"id":sample[6]})[0][1]
            #files = gateway("MUTECT", petljakapi.translate.idtostring(sample[0], "MPS"), scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = db)
            #files = [file for file in files if not file.contains("line_of_origin.txt")]
            ## Get analysis row
            analysis = petljakapi.select.multi_select(db = db, table = "analyses", filters = {"samples_id":sample[0], "pipeline_name":"MUTECT_CELLLINE", "reference_genome":wildcards.reference}, bench = config["bench"])
            if len(analysis) == 0:
                continue
            analysis = analysis[0]
            aid = petljakapi.translate.idtostring(analysis[0], "MPA")
            ap = analysis[10]
            ## Test if this sample is used as a reference for any other samples
            daughters = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "sample_parent_id", filter_value = sample[0], bench = config["bench"])
            parent_id = sample[5]
            if parent_id is not None and parent_id != "NULL": ## if this sample has a parent ie. it's a daughter
                files = [f"mutect2/{wildcards.reference}/std/table_raw.txt", f"mutect2/{wildcards.reference}/germ/table_raw.txt"]
                files = [os.path.join(ap, aid, f) for f in files]
                parent_name = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = parent_id)[0][1]
                with open(output.daughter, "a") as f:
                    f.write("%s\n" % '\t'.join([files[0], files[1], sample[1], parent_name, cell_name, sample[4], str(sample[7])]))
            if len(daughters) > 0 or parent_id == "NULL" or parent_id is None: ## ie. if this is a parent
                files = [os.path.join(ap, aid, "mutect2", wildcards.reference, "parental/table_raw.txt")]
                with open(output.parental, "a") as f:
                    f.write("%s\n" % '\t'.join([files[0], sample[1], cell_name, sample[4]]))



rule MUTECT2_CELLLINE_SBS_GROUP_FILTER:
    input:
        daughter = SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/{reference}/daughters_table.txt",
        parental = SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/{reference}/parents_table.txt",
        data = full_cell_line_cohort
    output:
        SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/{reference}/filtering_done.txt"
    threads: 4
    resources:
        mem_mb = 250000,
        slurm_partition = "fn_short",
        threads = 4,
        cpus = 4,
        runtime = 720,
    params:
        script_path = MUTECT_CELLLINE_PATH + "/scripts/filter_mutations.R"
    log:
        SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/{reference}/filtering.log"
    shell:
        """
        module load r/4.1.2;
        Rscript {params.script_path} -d {input.daughter} -p {input.parental} -t {resources.threads} 2> {log};
        echo 'Rscript {params.script_path} -d {input.daughter} -p {input.parental} -t {resources.threads}' > {output}
        """

rule PROC_FILE:
    input:
        filt = SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/{reference}/filtering_done.txt",
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/std/filtered_renamed.vcf"
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/proc/variants_final.vcf"
    threads: 1
    resources:
        mem_mb = 3000,
        slurm_partition =  config["clusters"][config["parts"]]["dev"],
        threads = 1,
        cpus = 1,
        runtime = 10,
    params:
        script_path = os.path.join(MUTECT_CELLLINE_PATH, "scripts/cellline_vcf.py"),
        proc_vcf = lambda w: SCRATCH_DIR + f"studies/{w.study}/samples/{w.sample}/analyses/MUTECT_CELLLINE/{w.analysis}/mutect2/{w.reference}/proc/variants.vcf"
    log:
         SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/proc/variants_final.log"
    shell:
        """
        cat <(python {params.script_path} -v {input.vcf} -c {input.filt}) {params.proc_vcf} > {output} &&
        singularity run -B /gpfs/ /gpfs/data/petljaklab/containers/gatk/gatk_4.4.0.0.sif gatk ValidateVariants -V {output} &> {log}
        """
        

rule CELLLINE_OF_ORIGIN:
    input:
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)[0],
        FA = lambda wildcards: FA_PATHS[wildcards.reference],
    output:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/proc/qc_genotyping.vcf",
        tbl1 =SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/proc/qc_genotyping.txt",
        tbl2 =SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/proc/line_of_origin.txt",
    singularity:
        "/gpfs/data/petljaklab/containers/bcftools/bcftools_latest.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/proc/qc_genotyping.log",
    resources:
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        runtime = 240,
        mem_mb = 3000
    shell:
        """
        bcftools mpileup \
         -R /gpfs/data/petljaklab/resources/{wildcards.reference}/pipeline_resources/somatic_celline/qc/1k_cells_genotyping/genotype_positions.txt \
         --fasta-ref {input.FA} \
         -A {input.cell_merge} 2> {log} | bcftools call -c 2>> {log} > {output.vcf};
        gatk VariantsToTable -V {output.vcf} -O {output.tbl1} 2>> {log};
        module load r/4.1.2;
        cat {output.tbl1} | Rscript scripts/qc_identity.R > {output.tbl2} 2>> {log}
        """

rule MUTECT2_CELLLINE_PAIRED_FILTER:
    input:
        germ = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/germ/table_raw.txt",
        std = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/std/table_raw.txt",
        tbl = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/parents.table"
    output:
        tab = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/variants.txt"
    resources:
        mem_mb = 100000,
        slurm_partition = "fn_short",
        runtime = 30
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/table.log"
    params:
        parent = lambda w: parent_sample_id(w)
    shell:
        "module load r/4.1.2; Rscript scripts/filter_paired.R --germline {input.germ} --input {input.std} --table {input.tbl} --parent {params.parent} --outtable {output.tab} 2> {log}"

## At the end of all this, the invocation of this module should be able to figure out if a sample is purely parental or not. If it's purely parental, then it shouldn't trigger the pipeline execution. 