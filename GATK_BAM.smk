GATK_BAM_PIPELINE_VER = "1.0.0"
GATK_ALIGNER_VER = "1.0.0"
GATK_VERSION = "4.4.0.0"

rule UBAM:
    input:
        fq1 = SCRATCH_DIR + "data/{proj}/{samp}/{run}_1.fastq",
        fq2 = SCRATCH_DIR + "data/{proj}/{samp}/{run}_2.fastq",
    output:
        ubam = temp(SCRATCH_DIR + "data/{proj}/{samp}/{run}.unaln.bam")
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    log:
        SCRATCH_DIR + "data/{proj}/{samp}/{run}.unaln.bam.log"
    shell:
        """
        gatk FastqToSam \
            --FASTQ {input.fq1} \
            --FASTQ2 {input.fq2} \
            --OUTPUT {output.ubam} \
            --READ_GROUP_NAME {wildcards.run} \
            --SAMPLE_NAME {wildcards.samp} \
            --LIBRARY_NAME {wildcards.run} \
            --PLATFORM ILLUMINA > {log} 2>&1
        """

rule markadapters:
    input:
        ubam = SCRATCH_DIR + "data/{proj}/{samp}/{run}.unaln.bam"
    output:
        obam = temp(SCRATCH_DIR + "data/{proj}/{samp}/{run}.unaln.adaptmarked.bam"),
        metrics = temp(SCRATCH_DIR + "data/{proj}/{samp}/{run}.unaln.adaptmarked.metrics")
    log:
        SCRATCH_DIR + "data/{proj}/{samp}/{run}.unaln.adaptmarked.bam.log"
    params:
        tmp = lambda wildcards: TEMP_DIR + wildcards.run,
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    shell:
        """
        gatk MarkIlluminaAdapters \
            I={input.ubam} \
            O={output.obam} \
            M={output.metrics} \
            TMP_DIR={params.tmp} 2> {log}
        """

rule bwamem2_gatkmerge:
    input:
        bwa_idxbase = lambda wildcards: ALN_REFERENCES[wildcards.reference],
        fa_base = lambda wildcards: FA_PATHS[wildcards.reference],
        ubam = SCRATCH_DIR + "data/{proj}/{samp}/{run}.unaln.bam",
        bam = SCRATCH_DIR + "data/{proj}/{samp}/{run}.unaln.adaptmarked.bam"
    output:
        bam = temp(SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_{version}/{reference}/{run}.bam")
    resources:
        threads = ALIGN_THREADS,
        mem_mb = 60000
    log:
        samtofastq = SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_{version}/{reference}/{run}.samtofastq.log",
        bwa = SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_{version}/{reference}/{run}.bwa.log",
        mergebamalignment = SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_{version}/{reference}/{run}.mergebamalignment.log"
    params:
        align_threads = ALIGN_THREADS - 2,
        tmp = lambda wildcards: TEMP_DIR + wildcards.run,
    singularity: f"/gpfs/data/petljaklab/containers/gatk-alignment/{GATK_ALIGNER_VER}/gatk-alignment_{GATK_ALIGNER_VER}.sif"
    shell:
        """
        gatk SamToFastq \
            I={input.bam} \
            FASTQ=/dev/stdout \
            CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
            TMP_DIR={params.tmp} 2> {log.samtofastq} | \
        bwa-mem2 mem -M -t {params.align_threads} -p {input.bwa_idxbase} /dev/stdin 2> {log.bwa} | \
        gatk MergeBamAlignment \
            ALIGNED_BAM=/dev/stdin \
            UNMAPPED_BAM={input.ubam} \
            OUTPUT={output.bam} \
            R={input.fa_base} CREATE_INDEX=true ADD_MATE_CIGAR=true \
            CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
            INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
            TMP_DIR={params.tmp} 2> {log.mergebamalignment}
        """

rule gatk_markdups:
    input:
        bam = SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_{version}/{reference}/{run}.bam"
    output:
        bam = SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_{version}/{reference}/{run}.dupsflagged.bam",
        metrics = SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_{version}/{reference}/{run}.dupsflagged.bam.metrics"
    resources:
        threads = ALIGN_THREADS,
        mem_mb = 20000
    log:
        SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_{version}/{reference}/{run}.dupsflagged.bam.log"
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    shell:
        """
        gatk MarkDuplicatesSpark \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --conf 'spark.executor.cores={resources.threads}' 2> {log}
        """
rule gatkmerge:
    input:
        lambda wildcards: expand(SCRATCH_DIR + "data/{{proj}}/{{samp}}/bwa-mem2_{{version}}/{{reference}}/{runid}.dupsflagged.bam", runid = accessions[wildcards.proj][wildcards.samp])
    output:
        merge = SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_{version}/merge/{samp}.merged.bam"
    log:
        SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_{version}/{samp}.merged.bam.log"
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{MUTECT_VERSION}.sif"
    resources:
        threads = 2
    params:
        inputlist = lambda wildcards: expand("-I {SCRATCH_DIR}data/{proj}/{samp}/bwa-mem2_{version}/{reference}/{runid}.dupsflagged.bam ", SCRATCH_DIR = SCRATCH_DIR, proj = wildcards.proj, samp = wildcards.samp, version = wildcards.version, reference = wildcards.reference, runid = accessions[wildcards.proj][wildcards.samp])
    shell:
        """
        gatk MergeSamFiles \
            {params.inputlist} -O {output.merge} \
            --USE_THREADING true \
             2> {log}
        """

rule COPY_MERGE:
    input:
        SCRATCH_DIR + "data/{proj}/{samp}/bwa-mem2_{version}/merge/{samp}.merged.bam"
    output:
        "data/{proj}/{samp}/bwa-mem2_{version}/merge/{samp}.merged.bam"