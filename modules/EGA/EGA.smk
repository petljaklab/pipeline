EGA_PIPE_VERSION = "1.0.0"
EGA_CONFIG = "/gpfs/data/petljaklab/lculibrk_prj/pipeljak/credentials_files/ega_credentials.txt"

def ega_function(wildcards):
    result = petljakapi.select.simple_select(db = db, table = "runs", filter_column = "id", filter_value = petljakapi.translate.stringtoid(wildcards.run), bench = config["bench"])
    ega_id = result[0][5]
    if not ega_id or ega_id == "NULL":
        raise ValueError(f"EGA pipeline is being used, but run {wildcards.run} has no associated EGA ID in the database")
    cram_path = SCRATCH_DIR + f"studies/{wildcards.study}/samples/{wildcards.sample}/runs/{wildcards.run}/analyses/EGA/{wildcards.analysis}/pyega/{ega_id}.cram"
    return(cram_path)


rule EGA_FETCH:
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/pyega/{ega}.cram"
    singularity:
        "/gpfs/data/petljaklab/containers/ega/ega_1.0.0.sif"
    resources:
        slurm_partition = "cpu_short",
        runtime = 720
    params:
        CREDENTIALS = EGA_CONFIG,
        out_dir = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/pyega/"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/pyega/{ega}.log"
    shell:
        """
        pyega3 -cf {params.CREDENTIALS} fetch {wildcards.ega} --output-dir {params.out_dir} &> {log}; mv {params.out_dir}/{wildcards.ega}/*.cram {output}
        if [[ $(md5sum {output} | awk '{{print $1}}') != $(cat {params.out_dir}/{wildcards.ega}/*.md5) ]]; then
            echo "md5 hash mismatch, throwing error"
            exit 69
        fi
        """
rule EGA_NAMESORT:
    input:
        ega_function
    output:
        temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/pyega/NAMESORT.bam")
    resources:
        slurm_partition = "cpu_short",
        mem_mb = 16000,
        runtime = 720,
        threads = 4,
        cpus = 4,
    threads: 4
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/pyega/NAMESORT.log"
    singularity:
        "/gpfs/data/petljaklab/containers/samtools/samtools_1.18.sif"
    params:
        ref = "/gpfs/data/petljaklab/resources/other_references/NCBI_hg19/GRCh37_alignment_nomask.fa",
        ref2 = "/gpfs/data/petljaklab/resources/other_references/hs37d5/hs37d5.fa",
        tmp = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/collate",
    shell:
        """
        set +e; \
        set +o pipefail; \
        rm -rf {params.tmp}; \
        mkdir -p {params.tmp}; \
        samtools collate -@ {threads} -u --reference {params.ref2} --output-fmt BAM {input} -o {output} -T {params.tmp}/collate &> {log}; \
        if [ ! $? == 0 ]; then \
            rm -rf {params.tmp}; \
            mkdir -p {params.tmp}; \
            samtools collate -@ {threads} -u --reference {params.ref} --output-fmt BAM {input} -o {output} -T {params.tmp}/collate &> {log}; \
            if [ ! $? == 0 ]; then \
                set -euo pipefail; \
                rm -rf {params.tmp}; \
                exit 69; \
            fi; \
        fi; \
        rm -rf {params.tmp};
        """

rule EGA_FASTQ:
    input:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/pyega/NAMESORT.bam"
    output:
        reads1 = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/fq/reads1.fq.u"),
        reads2 = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/fq/reads2.fq.u")
    singularity:
        "/gpfs/data/petljaklab/containers/samtools/samtools_1.18.sif"
    resources:
        slurm_partition = "cpu_short",
        runtime = 720,
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/fq/fastq.log"
    shell:
        """
        samtools fastq -1 {output.reads1} -2 {output.reads2} -0 /dev/null {input} 2> {log};
        """

rule EGA_COMPLETE:
    input:
        reads1 = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/fq/reads1.fq.u",
        reads2 = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/fq/reads2.fq.u"
    output:
        reads1 = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/fq/reads1.fq"),
        reads2 = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/fq/reads2.fq")
    params:
        DB = db,
    resources:
        iotasks = 2,
        slurm_partition = "cpu_short",
        runtime = 700,
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EGA/{analysis}/fq/fastq.log"
    shell:
        """
        mv {input.reads1} {output.reads1}; mv {input.reads2} {output.reads2}; \
        python scripts/mark_complete.py -i {wildcards.analysis} -d {params.DB} {output.reads1} 2>> {log}
        """


