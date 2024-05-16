import os
import numpy as np

MUTECT_CELLLINE_PIPELINE_VERSION = "1.0.0"
MUTECT_CELLLINE_PATH = os.path.join(basedir, "modules", "MUTECT_CELLLINE")


include: os.path.join(MUTECT_CELLLINE_PATH, "MUTECT2.smk")
include: os.path.join(MUTECT_CELLLINE_PATH, "PROCESS_FILTER_M2.smk")
include: os.path.join(MUTECT_CELLLINE_PATH, "MUTECT2_PARENTAL.smk")


def parent_cell(wildcards):
    sample_id = petljakapi.translate.stringtoid(wildcards.sample)
    db_line = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = sample_id)
    #print(db_line)
    parent_id = petljakapi.translate.idtostring(db_line[0][5], "MPS")
    parent_merge = gateway("WGS_MERGE_BAM", parent_id, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb")[0]
    return(parent_merge)

rule MUTECT2_CELLLINE_VCFTOTABLE:
    input:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/filtered.vcf",
    output:
        tab = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/{type}/table_raw.txt"
    resources:
        mem_mb = 3000,
        slurm_partition = "cpu_short",
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
    daughter_line = petljakapi.select.multi_select(db = db, table = "samples", filters = {"id":daughter_id}, headers = True)
    ## Get the line for that daughter's annotated parent and get the cell ID from that
    parent_id = np.array(daughter_line[1])[np.where(np.array(daughter_line[0]) == "sample_parent_id")][0]
    cell_id = np.array(daughter_line[1])[np.where(np.array(daughter_line[0]) == "cell_id")][0]
    ## Now get the parent of the parent, so we can get all the possible parents
    parental_sample = petljakapi.select.multi_select(db = db, table = "samples", filters = {"id":parent_id}, headers = True)
    parental_parent = np.array(parental_sample[1])[np.where(np.array(daughter_line[0]) == "sample_parent_id")][0]
    ## Now get all the appropriate parental lines
    all_parental_lines = petljakapi.select.multi_select(db = db, table = "samples", filters = {"study_id":study_id, "sample_parent_id":parental_parent, "cell_id":cell_id})
    ## Loop through to get the analysis paths and thereby the path to the parents' variant table
    parental_paths = []
    for parent in all_parental_lines:
        ## First gateway to ensure the parent is made
        gateway("MUTECT", petljakapi.translate.idtostring(parent[0], "MPS"), SCRATCH_DIR, config["PROD_DIR"], db = "petljakdb")
        ## Get the analysis dir
        analysis = petljakapi.select.multi_select(db = db, table = "analyses", filters = {"samples_id":parent[0], "pipeline_name":"MUTECT_CELLLINE"})
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
        slurm_partition = "cpu_dev",
        runtime = 5
    shell:
        """
        for i in {input}; do echo $i >> {output}; done;
        """

def parent_sample_id(wildcards):
    sample_id = petljakapi.translate.stringtoid(wildcards.sample)
    db_line = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = sample_id)
    #print(db_line)
    parent_id = petljakapi.translate.idtostring(db_line[0][5], "MPS")
    return(parent_id)

def full_cell_line_cohort(wildcards):
    study_id = petljakapi.translate.stringtoid(wildcards.study)
    samples_ids = petljakapi.select.multi_select(db, "samples", {"study_id":study_id})
    o_list = []
    for sample in samples_ids:
        o_list.extend(gateway("MUTECT", petljakapi.translate.idtostring(sample[0], "MPS"), scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = db))
    return(o_list)
    
rule M2_SBS_TABLES:
    output:
        daughter = SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/daughters_table.txt",
        parental = SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/parents_table.txt",
    threads: 1
    resources:
        slurm_partition = "cpu_dev",
        cpus = 1,
        threads = 1,
        mem_mb = 2000,
        runtime = 15,
    run:
        if not os.path.exists(f"{SCRATCH_DIR}studies/{wildcards.study}/analyses/MUTECT_CELLLINE/"):
            os.makedirs(f"{SCRATCH_DIR}studies/{wildcards.study}/analyses/MUTECT_CELLLINE/")
        study_id = petljakapi.translate.stringtoid(wildcards.study)
        samples_ids = petljakapi.select.multi_select(db, "samples", {"study_id":study_id})
        for sample in samples_ids:
            cell_name = petljakapi.select.multi_select(db, "cells", {"id":sample[6]})[0][1]
            files = gateway("MUTECT", petljakapi.translate.idtostring(sample[0], "MPS"), scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = db)
            if len(files) == 1: ## Parental
                with open(output.parental, "a") as f:
                    f.write("%s\n" % '\t'.join([files[0], sample[1], cell_name]))
            elif len(files) == 2: ## daughter
                parent_id = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = sample[0])[0][5]
                parent_name = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = parent_id)[0][1]
                with open(output.daughter, "a") as f:
                    f.write("%s\n" % '\t'.join([files[0], files[1], sample[1], parent_name, cell_name]))

rule MUTECT2_CELLLINE_SBS_GROUP_FILTER:
    input:
        daughter = SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/daughters_table.txt",
        parental = SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/parents_table.txt",
        data = full_cell_line_cohort
    output:
        SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/filtering_done.txt"
    resources:
        mem_mb = 150000,
        slurm_partition = "fn_short",
        threads = 4,
        cpus = 4,
        runtime = 480,
    params:
        script_path = MUTECT_CELLLINE_PATH + "/scripts/filter_mutations.R"
    log:
        SCRATCH_DIR + "studies/{study}/analyses/MUTECT_CELLLINE/filtering.log"
    shell:
        """
        module load r/4.1.2;
        Rscript {params.script_path} -d {input.daughter} -p {input.parental} 2> {log};
        touch {output};
        """

rule CELLLINE_OF_ORIGIN:
    input:
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb")[0],
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
        slurm_partition = "cpu_short",
        runtime = 240,
        mem_mb = 3000
    shell:
        """
        bcftools mpileup \
         -R /gpfs/data/petljaklab/resources/hg19/pipeline_resources/somatic_celline/qc/1k_cells_genotyping/genotype_positions.txt \
         --fasta-ref {input.FA} \
         -A {input.cell_merge} 2> {log} | bcftools call -c 2>> {log} > {output.vcf};
        gatk VariantsToTable -V {output.vcf} -O {output.tbl1} 2>> {log};
        module load r/4.1.2;
        cat {output.tbl1} | Rscript scripts/parse_variants.R > {output.tbl2} 2>> {log}
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