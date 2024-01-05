GATK_BAM_PIPELINE_VER = "1.0.0"

GATK_ALIGNER_VER = "1.0.0"
GATK_VERSION = "4.4.0.0"
SAMTOOLS_VERSION = "1.18"

def aggregate_runs(wildcards) -> list:
    """
    Takes wildcards as input and generates a list of expected filepaths for the rule GATK_MERGE

    Parameters:
    wildcards (snakemake.wildcards): Wildcards object from snakemake

    Returns:
    list: List of aligned.dupsflagged.bam files associated with each run of the sample in wildcards.
    """
    ## First get the numeric ID for sample
    sample_id = petljakapi.translate.stringtoid(wildcards.sample)
    ## API call for all the runs with that sample id
    result = petljakapi.select.simple_select(db = db, table = "runs", filter_column = "sample_id", filter_value = sample_id)
    ## If the sample isn't in the db, throw an error
    if not result:
        raise ValueError(f"No runs found for the requested sample {wildcards.sample}")
    ## Create list of MPR IDs
    run_ids = [petljakapi.translate.idtostring(line[0], "MPR") for line in result]
    ## Create list of paths that will be the input for this rule
    runpaths = expand(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam", run = run_ids, study = wildcards.study, sample = wildcards.sample, analysis = wildcards.analysis, reference = wildcards.reference)
    return(runpaths)

rule UBAM:
    ## Generates unaligned bams for GATK pipeline
    input:
        lambda wildcards: gateway("FASTQ", wildcards.run, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb_devel")
    output:
        ubam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.bam")
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.bam.log"
    shell:
        """
        gatk FastqToSam \
            --FASTQ {input[0]} \
            --FASTQ2 {input[1]} \
            --OUTPUT {output.ubam} \
            --READ_GROUP_NAME {wildcards.run} \
            --SAMPLE_NAME {wildcards.sample} \
            --LIBRARY_NAME {wildcards.run} \
            --PLATFORM ILLUMINA > {log} 2>&1
        """

rule MARK_ADAPTERS:
    ## Marks reads with adapter sequences
    input:
        ubam = rules.UBAM.output
    output:
        obam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.adaptmarked.bam"),
        metrics = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.adaptmarked.metrics"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.adaptmarked.bam.log"
    params:
        tmp = lambda wildcards: TEMP_DIR + wildcards.run,
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    shell:
        """
        gatk MarkIlluminaAdapters \
            -I {input.ubam} \
            -O {output.obam} \
            -M {output.metrics} \
            -TMP_DIR {params.tmp} 2> {log}
        """

rule BWAMEM2_GATKMERGE:
    ## Does a few things. First converts the adapter-marked bam to fastq for alignment, then aligns them to the reference, then performs MergeBamAlignment to merge
    ## aligned/unaligned reads (I think) from the unaligned bam and the newly aligned sam and output a bam
    input:
        bwa_idxbase = lambda wildcards: ALN_REFERENCES[wildcards.reference],
        fa_base = lambda wildcards: FA_PATHS[wildcards.reference],
        ubam = rules.UBAM.output,
        bam = rules.MARK_ADAPTERS.output.obam,
    output:
        bam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.bam")
    resources:
        threads = ALIGN_THREADS,
        mem_mb = 60000
    log:
        samtofastq = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.samtofastq.log",
        bwa = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.bwa.log",
        mergebamalignment = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.mergebamalignment.log"
    params:
        align_threads = ALIGN_THREADS - 2,
        tmp = lambda wildcards: TEMP_DIR + wildcards.run,
    singularity: f"/gpfs/data/petljaklab/containers/gatk-alignment/{GATK_ALIGNER_VER}/gatk-alignment_{GATK_ALIGNER_VER}.sif"
    shell:
        """
        gatk SamToFastq \
            -I {input.bam} \
            -FASTQ /dev/stdout \
            -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true \
            -TMP_DIR {params.tmp} 2> {log.samtofastq} | \
        bwa-mem2 mem -M -t {params.align_threads} -p {input.bwa_idxbase} /dev/stdin 2> {log.bwa} | \
        gatk MergeBamAlignment \
            -ALIGNED_BAM /dev/stdin \
            -UNMAPPED_BAM {input.ubam} \
            -OUTPUT {output.bam} \
            -R {input.fa_base} -CREATE_INDEX true -ADD_MATE_CIGAR true \
            -CLIP_ADAPTERS false -CLIP_OVERLAPPING_READS true \
            -INCLUDE_SECONDARY_ALIGNMENTS true -MAX_INSERTIONS_OR_DELETIONS -1 \
            -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS \
            -TMP_DIR {params.tmp} 2> {log.mergebamalignment}
        """

rule GATK_MARKDUPS:
    ## Mark duplicate reads/alignments
    input:
        bam = rules.BWAMEM2_GATKMERGE.output,
    output:
        bam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam"),
        metrics = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam.metrics"
    resources:
        threads = ALIGN_THREADS,
        mem_mb = 20000
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam.log"
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    shell:
        """
        gatk MarkDuplicatesSpark \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --conf 'spark.executor.cores={resources.threads}' 2> {log}
        """


rule GATK_MERGE:
    ## Merge multi-run samples into single bams
    ## Currently merges ALL bam runs from a sample into a single bam. If we want subsets, will need to add that functionality later
    input:
        aggregate_runs
    output:
        merge = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam")
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam.log"
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{MUTECT_VERSION}.sif"
    resources:
        threads = 2
    params:
        inputlist = lambda wildcards, input: "-I ".join([input]) if isinstance(input, str) else "-I ".join(input)
    shell:
        """
        gatk MergeSamFiles \
            {params.inputlist} -O {output.merge} \
            --USE_THREADING true \
             2> {log}
        """

rule INDEX_GATKMERGE:
    ## Index the merge
    input:
        rules.GATK_MERGE.output
    output:
        temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam.bai")
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam.log"
    singularity: f"/gpfs/data/petljaklab/containers/samtools/samtools_{SAMTOOLS_VERSION}.sif"
    shell:
        "samtools index {input} 2> {output}"

rule copy_files:
    ## Copy to production path
    input:
        bam = rules.GATK_MERGE.output,
        bai = rules.INDEX_GATKMERGE.output,
    output:
        bam = PROD_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam",
        bai = PROD_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam.bai"
    shell:
        "rsync {input.bam} {output.bam}; rsync {input.bai} {output.bai} "