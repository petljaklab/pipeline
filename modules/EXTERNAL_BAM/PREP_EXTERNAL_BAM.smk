import logging
import re
import shutil
import os

def get_external_bam_path(wildcards):
    if wildcards.study == "MPP000003":
        basepath = "/gpfs/data/petljaklab/broad_data/bams/Bams/"
        iden = petljakapi.translate.stringtoid(wildcards.sample)
        res = petljakapi.select.multi_select(db = db, table = "samples", filters = {"id":iden}, bench = config["bench"])[0]
        rname = res[1]
        path = os.path.join(basepath, rname + ".bam")
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