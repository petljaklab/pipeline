import petljakapi.select
import petljakapi.translate


def get_prod_bam(wildcards):
    samp_id = petljakapi.translate.stringtoid(wildcards.sample)
    ret = petljakapi.select.multi_select(db = config['db'], table = "analyses", filters = {"pipeline_name":"GATK_BAM", "analysis_complete":"True", "samples_id":samp_id})
    print(ret)
    ret = ret[0]
    ret_path = glob.glob(os.path.join(ret[10], "merged.cram"))
    return(ret_path)

rule SPLIT_PROD_BAM:
    input:
        get_prod_bam
    output:
        temp(directory(SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM_CONVERT/split_bams/"))
    singularity:
        "/gpfs/data/petljaklab/containers/samtools/samtools_1.18.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM_CONVERT/split_bams/split.log"
    threads: 8
    priority: 1
    resources:
        threads = 8,
        cpus = 8,
        mem_mb = 8000,
        runtime = 240,
        slurm_partition = "petljaklab,cpu_short",
        iotasks = 5,
    shell:
        "rm -rf {output}; samtools split {input} -f {output}/{wildcards.sample}-%\!.cram -@ {resources.threads}"

rule PROD_SPLIT_BAM_TO_FQ:
    input:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM_CONVERT/split_bams/"
    output:
        reads1 = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM_CONVERT/{analysis}/fq/reads1.fq"),
        reads2 = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM_CONVERT/{analysis}/fq/reads2.fq"),
        readsS = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM_CONVERT/{analysis}/fq/readsS.fq",
        readsU = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM_CONVERT/{analysis}/fq/readsU.fq",
    singularity:
        f"/gpfs/data/petljaklab/containers/samtools/samtools_1.18.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/split.log",
    resources:
        iotasks = 2,
        mem_mb = 20000,
        runtime = 600,
        slurm_partition = "petljaklab,cpu_short"
    priority: 1
    shell:
        "samtools collate -u -O {input}/{wildcards.sample}-{wildcards.run}.cram | samtools fastq -1 {output.reads1} -2 {output.reads2} -0 {output.readsU} -s {output.readsS} - &> {log};"
