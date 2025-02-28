import logging
import re
import shutil
import os

def get_external_bam_path(wildcards):
    """
    Return the path to a bam file given the sample ID. 

    Each study ID with external_bam files must have a section here that takes in the study/sample wildcards and returns 
    an absolute path in the fileststem that points to the external bam file. 
    """
    if wildcards.study == "MPP000003":
        ## The data are all in this directory
        basepath = "/gpfs/data/petljaklab/broad_data/bams/Bams/"
        ## Get the numeric ID for the current sample by translating the string (e.g. MPS000256) to an integer ID (256)
        iden = petljakapi.translate.stringtoid(wildcards.sample)
        ## Perform a database query for this ID, return the rname as this corresponds to the file names in basepath
        res = petljakapi.select.multi_select(db = db, table = "samples", filters = {"id":iden}, bench = config["bench"])[0]
        rname = res[1]
        ## assign path to the path to the sample's bam file and return it. 
        path = os.path.join(basepath, rname + ".bam")
    ## Add sections for different study IDs

    ## Return the path to the bam file for the given sample ID
    return(path)

localrules: COUNT_RG_ASSIGN_RUNS

rule COUNT_RG_ASSIGN_RUNS:
    input:
        get_external_bam_path
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/LOAD_EXT_BAM/{analysis}/db_load_runs/runs.loaded"
    params:
        workdir = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/LOAD_EXT_BAM/{analysis}/db_load_runs"
    shell:
        """
            samtools view -H {input} > {params.workdir}/header.txt
            NRUNS=$(grep '^@RG' {params.workdir}/header.txt | wc -l)
            for ((i = 0; i < NRUNS; i++)); do
                python scripts/add_element.py -t run -r {wildcards.sample}-$i --study {wildcards.study} --sample {wildcards.sample} --db petljakdb --source external_bam
                echo "Iteration $((i + 1))"
            done
            touch {output}
        """