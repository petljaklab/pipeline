def get_external_bam_path(wildcards):
    ## Run ID
    runid = wildcards.run
    runid = petljakapi.translate.stringtoid(runid)
    ## Query to get local path
    
    return(path)

rule EXTERNAL_BAM_TO_FASTQ:
    input:
        get_external_bam_path
    output:
        reads1 = PROD_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/reads1.fq",
        reads2 = PROD_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/reads2.fq",
        readsU = PROD_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/readsU.fq",
    singularity:
        "/gpfs/data/petljaklab/containers/samtools/samtools_1_18.sif"
    log:
        PROD_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/samtofastq.log",
    shell:
        "gatk SamToFastq -I {input} -F {output.reads1} -F2 {output.reads2} -FU {output.readsU} &> {log}"