import os

GATK_BAM_PIPELINE_VER = "1.0.0"

GATK_ALIGNER_VER = "1.0.1"
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
        lambda wildcards: gateway("FASTQ", wildcards.run, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = config["db"])
    output:
        ubam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.bam")
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.bam.log"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.benchmark"
    resources:
        runtime = 120,
        iotasks = 2,
        slurm_partition = "cpu_dev,cpu_short,fn_short"
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
        tmp = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/markadapts/",
    benchmark: 
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.adaptmarked.benchmark"
    resources:
        runtime = 120,
        iotasks = 2,
        slurm_partition = "cpu_dev,cpu_short,fn_short"
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    shell:
        """
        mkdir -p {params.tmp}
        gatk MarkIlluminaAdapters \
            -I {input.ubam} \
            -O {output.obam} \
            -M {output.metrics} \
            -TMP_DIR {params.tmp} 2> {log}
        rm -rf {params.tmp}
        """

rule MARKED_SAM_TO_FASTQ:
    input:
        bam = rules.MARK_ADAPTERS.output.obam,
    output:
        ifastq = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/interleaved.fastq")
    log:
        samtofastq = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.samtofastq.log",
    params:
        align_threads = ALIGN_THREADS - 2,
        tmp = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/marktofq/",
    singularity: f"/gpfs/data/petljaklab/containers/gatk-alignment/gatk-alignment_{GATK_ALIGNER_VER}.sif"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.samtofastq.benchmark"
    resources:
        runtime = 480,
        iotasks = 2,
        slurm_partition = "cpu_dev,cpu_short,fn_short"
    shell:
        """
        mkdir -p {params.tmp}
        gatk SamToFastq \
            -I {input.bam} \
            -FASTQ {output.ifastq} \
            -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true \
            -TMP_DIR {params.tmp} 2> {log.samtofastq}
        rm -rf {params.tmp}
        """

rule BWA:
    input:
        ifastq = rules.MARKED_SAM_TO_FASTQ.output.ifastq,
        bwa_idxbase = lambda wildcards: ALN_REFERENCES[wildcards.reference],
    output:
        rawbam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.raw.sam")
    threads: ALIGN_THREADS
    resources:
        threads = ALIGN_THREADS,
        cpus = ALIGN_THREADS,
        mem_mb = 30000,
        iotasks = 4,
        slurm_partition = "cpu_medium,fn_medium"
    log:
        bwa = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.bwa.log",
    singularity: f"/gpfs/data/petljaklab/containers/gatk-alignment/gatk-alignment_{GATK_ALIGNER_VER}.sif"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.bwa.benchmark"
    params:
        tmpfastqpath = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/fastq/",
        tmpfastq = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/fastq/input.fastq"
    shell:
        """
            mkdir -p {params.tmpfastqpath}
            cp {input.ifastq} {params.tmpfastq};
            bwa mem -M -t {threads} -p {input.bwa_idxbase} {params.tmpfastq} > {output} 2> {log.bwa};
            rm {params.tmpfastq}
        """

rule MERGEBAMALIGNMENT:
    input:
        bwa = rules.BWA.output.rawbam,
        fa_base = lambda wildcards: FA_PATHS[wildcards.reference],
        ubam = rules.UBAM.output,
    output:
        bam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.bam")
    log:
        mergebamalignment = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.mergebamalignment.log"
    params:
        tmp = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/mergebamalignment/",
    resources:
        threads = 1,
        cpus = 1,
        mem_mb = 10000,
        iotasks = 4,
        slurm_partition = "cpu_medium,fn_medium"
    singularity: f"/gpfs/data/petljaklab/containers/gatk-alignment/gatk-alignment_{GATK_ALIGNER_VER}.sif"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.mergebamalignment.benchmark"
    shell:
        """
        mkdir -p {params.tmp}
        gatk MergeBamAlignment \
            -ALIGNED_BAM {input.bwa} \
            -UNMAPPED_BAM {input.ubam} \
            -OUTPUT {output.bam} \
            -R {input.fa_base} -CREATE_INDEX true -ADD_MATE_CIGAR true \
            -CLIP_ADAPTERS false -CLIP_OVERLAPPING_READS true \
            -INCLUDE_SECONDARY_ALIGNMENTS true -MAX_INSERTIONS_OR_DELETIONS -1 \
            -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS \
            -TMP_DIR {params.tmp} 2> {log.mergebamalignment}
        rm -rf {params.tmp}
        """


rule BWAMEM2_GATKMERGE:
    ## Does a few things. First converts the adapter-marked bam to fastq for alignment, then aligns them to the reference, then performs MergeBamAlignment to merge
    ## aligned/unaligned reads (I think) from the unaligned bam and the newly aligned sam and output a bam
    ## Deprecated because it doesn't play nice with memory
    input:
        bwa_idxbase = lambda wildcards: ALN_REFERENCES[wildcards.reference],
        fa_base = lambda wildcards: FA_PATHS[wildcards.reference],
        ubam = rules.UBAM.output,
        bam = rules.MARK_ADAPTERS.output.obam,
    output:
        bam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.bama")
    threads: ALIGN_THREADS
    resources:
        threads = ALIGN_THREADS,
        cpus = ALIGN_THREADS,
        mem_mb = 20000,
        iotasks = 4,
        slurm_partition = "cpu_medium,fn_medium"
    log:
        samtofastq = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.samtofastq.log",
        bwa = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.bwa.log",
        mergebamalignment = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.mergebamalignment.log"
    params:
        align_threads = ALIGN_THREADS - 2,
        tmp = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis,
    singularity: f"/gpfs/data/petljaklab/containers/gatk-alignment/gatk-alignment_{GATK_ALIGNER_VER}.sif"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.benchmark"
    shell:
        """
        gatk SamToFastq \
            -I {input.bam} \
            -FASTQ /dev/stdout \
            -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true \
            -TMP_DIR {params.tmp} 2> {log.samtofastq} | \
        bwa mem -M -t {params.align_threads} -p {input.bwa_idxbase} /dev/stdin 2> {log.bwa} | \
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
        bam = rules.MERGEBAMALIGNMENT.output,
    output:
        bam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam"),
#        metrics = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam.metrics"
    threads: 1
    resources:
        threads = 1,
        cpus = ALIGN_THREADS,
        mem_mb = 20000,
        iotasks = 4,
        slurm_partition = "cpu_short,fn_short",
        runtime = 720
    params:
        TEMP_DIR = TEMP_DIR,
        TEMP_BAM_PATH = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis +  "/markdups/",
        TEMP_BAM_FILE = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/markdups/input.bam"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam.log"
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.benchmark"
    shell:
        """
        rm -rf {output.bam}.parts/;
        mkdir -p {params.TEMP_BAM_PATH};
        cp {input.bam} {params.TEMP_BAM_FILE};
        gatk MarkDuplicatesSpark \
            -I {params.TEMP_BAM_FILE} \
            -O {output.bam} \
            --conf 'spark.executor.cores={resources.threads}' \
            --conf 'spark.local.dir={params.TEMP_DIR}' &> {log};
        rm {params.TEMP_BAM_FILE}
        """

rule GATK_MERGE:
    ## Merge multi-run samples into single bams
    ## Currently merges ALL bam runs from a sample into a single bam. If we want subsets, will need to add that functionality later
    input:
        aggregate_runs
    output:
        merge = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam.log"
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{MUTECT_VERSION}.sif"
    threads: 2
    resources:
        threads = 2,
        cpus = 2,
        runtime = 480,
        iotasks = 4,
        slurm_partition = "cpu_short,fn_short,cpu_dev"
    params:
        inputlist = lambda wildcards, input: f"-I {input}" if isinstance(input, str) else "-I " + " -I ".join(input)
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.benchmark"
    shell:
        """
        gatk MergeSamFiles \
            {params.inputlist} -O {output.merge} \
            --USE_THREADING true \
             &>> {log}
        """

rule INDEX_GATKMERGE:
    ## Index the merge
    input:
        rules.GATK_MERGE.output
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam.bai"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam.log"
    singularity: f"/gpfs/data/petljaklab/containers/samtools/samtools_{SAMTOOLS_VERSION}.sif"
    resources:
        runtime = 240,
        slurm_partition = "cpu_short"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bai.benchmark"
    shell:
        "samtools index {input} &>> {log}"

rule GATK_BAM_DONE:
    input:
        rules.GATK_MERGE.output,
        rules.INDEX_GATKMERGE.output,
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.done"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam.log"
    params:
        db = config["db"]
    resources:
        runtime = 10,
        slurm_partition = "cpu_short",
    shell:
        "python scripts/mark_complete.py --id {wildcards.analysis} --db {params.db} {input} >> {log}; touch {output}"

rule all_at_once:
    input:
        lambda wildcards: gateway("FASTQ", wildcards.run, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = config["db"]),
        bwa_idxbase = lambda wildcards: ALN_REFERENCES["hg19"],
        fa_base = lambda wildcards: FA_PATHS["hg19"],
    output:
        ubam = temp(TEMP_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.bam"),
        obam = temp(TEMP_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.adaptmarked.bam"),
        adaptmetrics = temp(TEMP_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.adaptmarked.metrics"),
        abam = temp(TEMP_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/hg19/aligned.bam"),
        dupsflagged = temp(TEMP_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/hg19/aligned.dupsflagged.bam"),
        dupmetrics = temp(TEMP_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/hg19/aligned.dupsflagged.bam.metrics"),
        finalbam = "/gpfs/data/petljaklab/lculibrk_prj/misc/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/hg19/aligned.dupsflagged.bam",
    singularity: f"/gpfs/data/petljaklab/containers/gatk-alignment/gatk-alignment_{GATK_ALIGNER_VER}.sif"
    benchmark: "/gpfs/data/petljaklab/lculibrk_prj/misc/studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/hg19/benchmark.file"
    resources:
        threads = ALIGN_THREADS,
        cpus = ALIGN_THREADS,
        mem_mb = 5000,
        runtime = 10,
        slurm_partition = "cpu_dev"
    params:
        tmpdir = TEMP_DIR,
        tmp = lambda wildcards: os.path.join(TEMP_DIR, wildcards.analysis),
        align_threads = ALIGN_THREADS - 2,
    log:
        ubam = TEMP_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.bam.log",
        markadapts = TEMP_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/unaln.adaptmarked.bam.log",
        samtofastq = TEMP_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/hg19/aligned.samtofastq.log",
        bwa = TEMP_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/hg19/aligned.bwa.log",
        mergebamalignment = TEMP_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/hg19/aligned.mergebamalignment.log",
        markdups = TEMP_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/hg19/aligned.dupsflagged.log",
    shell:
        """
        mkdir -p {params.tmp}
        gatk FastqToSam \
            --FASTQ {input[0]} \
            --FASTQ2 {input[1]} \
            --OUTPUT {output.ubam} \
            --READ_GROUP_NAME {wildcards.run} \
            --SAMPLE_NAME {wildcards.sample} \
            --LIBRARY_NAME {wildcards.run} \
            --PLATFORM ILLUMINA &> {log.ubam};
        gatk MarkIlluminaAdapters \
            -I {output.ubam} \
            -O {output.obam} \
            -M {output.adaptmetrics} \
            -TMP_DIR {params.tmp} &> {log.markadapts}
        gatk SamToFastq \
            -I {output.ubam} \
            -FASTQ /dev/stdout \
            -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true \
            -TMP_DIR {params.tmp} 2> {log.samtofastq} | \
        bwa mem -M -t {params.align_threads} -p {input.bwa_idxbase} /dev/stdin 2> {log.bwa} | \
        gatk MergeBamAlignment \
            -ALIGNED_BAM /dev/stdin \
            -UNMAPPED_BAM {output.ubam} \
            -OUTPUT {output.abam} \
            -R {input.fa_base} -CREATE_INDEX true -ADD_MATE_CIGAR true \
            -CLIP_ADAPTERS false -CLIP_OVERLAPPING_READS true \
            -INCLUDE_SECONDARY_ALIGNMENTS true -MAX_INSERTIONS_OR_DELETIONS -1 \
            -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS \
            -TMP_DIR {params.tmp} 2> {log.mergebamalignment};
        gatk MarkDuplicatesSpark \
            -I {output.abam} \
            -O {output.dupsflagged} \
            -M {output.dupmetrics} \
            --conf 'spark.executor.cores={resources.threads}' \
            --conf 'spark.local.dir={params.tmpdir}' &> {log.markdups};
        cp {output.dupsflagged} {output.finalbam};
        cp {params.tmpdir}studies/*/samples/*/runs/*/analyses/GATK_BAM/*/bam/*.log /gpfs/data/petljaklab/lculibrk_prj/misc/studies/{wildcards.study}/samples/{wildcards.sample}/runs/{wildcards.run}/analyses/GATK_BAM/{wildcards.analysis}/;
        cp {params.tmpdir}studies/*/samples/*/runs/*/analyses/GATK_BAM/*/bam/hg19/*.log /gpfs/data/petljaklab/lculibrk_prj/misc/studies/{wildcards.study}/samples/{wildcards.sample}/runs/{wildcards.run}/analyses/GATK_BAM/{wildcards.analysis}/;
        """
#rule copy_files:
#    ## Copy to production path
#    input:
#        bam = rules.GATK_MERGE.output,
#        bai = rules.INDEX_GATKMERGE.output,
#    output:
#        bam = PROD_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam",
#        bai = PROD_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.bam.bai"
#    shell:
#        "rsync {input.bam} {output.bam}; rsync {input.bai} {output.bai}"