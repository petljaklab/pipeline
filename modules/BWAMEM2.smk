BWAMEM2_PIPE_VERSION = "1.0.0"

rule ALN_BWAMEM2:
    input:
        bwa_idxbase = lambda w: ALN_REFERENCES[w.hgver],
        fq1 = SCRATCH_DIR + "data/{proj}/{samp}/{run}_1.fastq",
        fq2 = SCRATCH_DIR + "data/{proj}/{samp}/{run}_2.fastq",
    output:
        bam = temp(SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_" + ALIGNER_VERSION + "/{hgver}/{run}.bam"),
        bai = temp(SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_" + ALIGNER_VERSION + "/{hgver}/{run}.bam.bai")
    resources:
        threads = align_threads*2,
        mem_mb = 60000
    params:
        align_threads = align_threads,
        tmp = lambda wildcards: TEMP_DIR + wildcards.run

    singularity: f"/gpfs/data/petljaklab/containers/bwamem2/{ALIGNER_VERSION}/{ALIGNER}_{ALIGNER_VERSION}.sif"
    log:
        bwa = SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_" + ALIGNER_VERSION + "/{hgver}/{run}.bwa2.log",
        samtools = SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_" + ALIGNER_VERSION + "/{hgver}/{run}.samtools.log"
    benchmark:
        SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_{version}/{hgver}/{run}.resources"
    shell:
        """
        echo 'bwa2 version:' > {log.bwa}; \
        bwa-mem2 version >> {log.bwa}; \
        echo 'samtools version:' > {log.samtools}; \
        samtools --version >> {log.samtools}; \
        bwa-mem2 mem -t{params.align_threads} {input.bwa_idxbase} {input.fq1} {input.fq2}   2>> {log.bwa} \
         | samtools collate -u -@{params.align_threads} -T {params.tmp} - -O              2>> {log.samtools} \
         | samtools fixmate -u -@{params.align_threads} -m              - -               2>> {log.samtools} \
         | samtools sort    -u -@{params.align_threads} -T{params.tmp}  -                 2>> {log.samtools} \
         | samtools markdup    -@{params.align_threads} -T{params.tmp}  -   {output.bam}  2>> {log.samtools};
        samtools index -@{resources.threads} {output.bam} 2>> {log.samtools}
        """

