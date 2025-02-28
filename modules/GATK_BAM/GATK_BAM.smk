import os
import numpy as np

GATK_BAM_PIPELINE_VER = "1.0.0"

GATK_ALIGNER_VER = "1.0.1"
GATK_VERSION = "4.4.0.0"
SAMTOOLS_VERSION = "1.18"

ruleorder: UBAM > EXTERNAL_BAM_TO_FASTQ > SPLIT_BAM

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
    result = petljakapi.select.simple_select(db = db, table = "runs", filter_column = "sample_id", filter_value = sample_id, bench = config["bench"])
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
        runtime = 480,
        iotasks = 2,
        mem_mb = 4000,
        jv_mem = 3750,
        slurm_partition = config["clusters"][config["parts"]]["short"]
    params:
        TEMP_DIR = TEMP_DIR,
        TEMP_BAM_PATH = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis +  "/markdups/",
    priority: 9999
    shell:
        """
        gatk --java-options '-Xmx{resources.jv_mem}M -Dsamjdk.compression_level=6' \
            FastqToSam \
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
        runtime = 1200,
        mem_mb = 2000,
        jv_mem = 1900,
        iotasks = 2,
        slurm_partition = config["clusters"][config["parts"]]["short"]
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    priority: 1000
    shell:
        """
        mkdir -p {params.tmp}
        gatk --java-options '-Xmx{resources.jv_mem}M -Dsamjdk.compression_level=6' \
            MarkIlluminaAdapters \
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
        runtime = 1200,
        mem_mb = 4000,
        jv_mem = 3750,
        iotasks = 2,
        slurm_partition = config["clusters"][config["parts"]]["short"]
    priority: 1001
    shell:
        """
        mkdir -p {params.tmp}
        gatk --java-options '-Xmx{resources.jv_mem}M' \
            SamToFastq \
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
        slurm_partition = config["clusters"][config["parts"]]["med"],
        tmpdisk = lambda wc, input: int(np.round(input.size_mb))
    log:
        bwa = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.bwa.log",
    singularity: f"/gpfs/data/petljaklab/containers/gatk-alignment/gatk-alignment_{GATK_ALIGNER_VER}.sif"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.bwa.benchmark"
    params:
        tmpfastqpath = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/fastq/",
        tmpfastq = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/fastq/input.fastq"
    priority: 1002
    shell:
        """
            mkdir -p {params.tmpfastqpath};
            trap "rm -rf {params.tmpfastqpath}" EXIT;
            cp {input.ifastq} {params.tmpfastq};
            bwa mem -M -t {threads} -p {input.bwa_idxbase} {params.tmpfastq} > {output} 2> {log.bwa};
            rm -rf {params.tmpfastqpath}
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
        runtime = 1200,
        threads = 1,
        cpus = 1,
        mem_mb = lambda wildcards, attempt: 5000 * (1 + ((attempt-1)/2)),
        jv_mem = lambda wildcards, attempt: 4500 * (1 + ((attempt-1)/2)),
        iotasks = 4,
        slurm_partition = config["clusters"][config["parts"]]["med"]
    singularity: f"/gpfs/data/petljaklab/containers/gatk-alignment/gatk-alignment_{GATK_ALIGNER_VER}.sif"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.mergebamalignment.benchmark"
    priority: 9999
    shell:
        """
        mkdir -p {params.tmp}
        gatk --java-options '-Xmx{resources.jv_mem}M -Dsamjdk.compression_level=6' \
            MergeBamAlignment \
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

rule GATK_MARKDUPS:
    ## Mark duplicate reads/alignments
    input:
        rules.MERGEBAMALIGNMENT.output,
    output:
        bam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam"),
#        metrics = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam.metrics"
    threads: 1    
    params:
        TEMP_DIR = TEMP_DIR,
        TEMP_BAM_PATH = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis +  "/markdups/",
        TEMP_BAM_FILE = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/markdups/input.bam"
    resources:
        threads = 8,
        cpus = ALIGN_THREADS,
        mem_mb = 35000,
        iotasks = 4,
        slurm_partition = config["clusters"][config["parts"]]["med"],
        runtime = 60*24*4,
        tmpdisk = lambda wc, input: int(np.round(input.size_mb))
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam.log"
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/GATK_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.benchmark"
    priority: 1003
    shell:
        """
        set +eo pipefail;
        trap "rm -rf {output.bam}.parts/ {params.TEMP_BAM_PATH}" EXIT;
        mkdir -p {params.TEMP_BAM_PATH};
        cp {input} {params.TEMP_BAM_FILE};
        gatk --java-options '-Dsamjdk.compression_level=6' \
            MarkDuplicatesSpark \
            -I {params.TEMP_BAM_FILE} \
            -O {output.bam} \
            --conf 'spark.executor.cores={resources.threads}' \
            --conf 'spark.local.dir={params.TEMP_DIR}' &> {log};
        rm -rf {params.TEMP_BAM_PATH}
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
    threads: 2
    resources:
        threads = 2,
        cpus = 2,
        runtime = 960,
        iotasks = 4,
        slurm_partition = config["clusters"][config["parts"]]["short"]
    params:
        inputlist = lambda wildcards, input: f"-I {input}" if isinstance(input, str) else "-I " + " -I ".join(input)
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.benchmark"
    priority: 1004
    shell:
        """
        gatk --java-options '-Dsamjdk.compression_level=6' \
            MergeSamFiles \
            {params.inputlist} -O {output.merge} \
            --USE_THREADING true \
             &>> {log}
        """

rule BAM2CRAM:
    input:
        merge = rules.GATK_MERGE.output,
        bwa_idxbase = lambda wildcards: ALN_REFERENCES[wildcards.reference],
    output:
        cram = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.cram",
        crai = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.cram.crai",
        ref_md5 = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.ref.md5",
        readme = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/CRAM_README.txt",
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.cram.log"
    singularity: f"/gpfs/data/petljaklab/containers/samtools/samtools_{SAMTOOLS_VERSION}.sif"
    threads: 1
    resources:
        runtime = 600,
        slurm_partition = config["clusters"][config["parts"]]["short"],
        mem_mb = 4000,
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.cram.bench"
    priority: 1005
    shell:
        """
        md5sum {input.bwa_idxbase} > {output.ref_md5};
        echo    'This CRAM was generated using the fasta file located at {input.bwa_idxbase}. The md5 hash of the reference is at {output.ref_md5}. \
                The matching reference fasta at the specified directory is REQUIRED for proper decompression of the file' > {output.readme};
        samtools view -@ {threads} -C -T {input.bwa_idxbase} {input.merge} > {output.cram} 2> {log}
        samtools index {output.cram} &>> {log}
        chmod 750 {output};
        """

rule GATK_BAM_DONE:
    input:
        rules.BAM2CRAM.output,
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.done"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/GATK_BAM/{analysis}/merge/{reference}/merged.cram.log"
    params:
        db = config["db"]
    resources:
        runtime = 10,
        slurm_partition = config["clusters"][config["parts"]]["short"],
    priority: 1007
    shell:
        "python scripts/mark_complete.py --id {wildcards.analysis} --db {params.db} {input} >> {log}; touch {output}"

