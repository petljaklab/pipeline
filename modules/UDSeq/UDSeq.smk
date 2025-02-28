UDSEQ_BASE_PATH = os.path.join(basedir, "modules", "UDSeq")

def aggregate_udseq_runs(wildcards) -> list:
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
    result = petljakapi.select.multi_select(db, "runs", {"sample_id":sample_id, "sequencing_strategy":"UDSEQ"})
    ## If the sample isn't in the db, throw an error
    if not result:
        raise ValueError(f"No runs found for the requested sample {wildcards.sample}")
    ## Create list of MPR IDs
    run_ids = [petljakapi.translate.idtostring(line[0], "MPR") for line in result]
    ## Create list of paths that will be the input for this rule
    runpaths = expand(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam", run = run_ids, study = wildcards.study, sample = wildcards.sample, analysis = wildcards.analysis, reference = wildcards.reference)
    return(runpaths)

rule UDSEQ_TRIM:
    ## Trims fastq files
    input:
        lambda wildcards: gateway("FASTQ", wildcards.run, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = config["db"], ref = "hg38")
    output:
        trim1 = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/{run}_trimmed_1.fastq"),
        trim2 = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/{run}_trimmed_2.fastq")
    conda: f"{UDSEQ_BASE_PATH}/environment.yaml"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/unaln.bam.log"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/unaln.benchmark"
    resources:
        runtime = 240,
        iotasks = 2,
        mem_mb = 4000,
        jv_mem = 3750,
        slurm_partition = config["clusters"][config["parts"]]["short"],
    params:
        out_path = lambda wildcards: SCRATCH_DIR + f"studies/{wildcards.study}/samples/{wildcards.sample}/runs/{wildcards.run}/analyses/UDSEQ_BAM/{wildcards.analysis}/fastq/{wildcards.run}_trimmed"
    priority: 9999
    shell:
        """
        DupCaller.py trim -i {input[0]} -i2 {input[1]} -p "NNNXXXX" -o {params.out_path}
        """

rule UDSEQ_FASTP:
    input:
        r1 = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/{run}_trimmed_1.fastq",
        r2 = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/{run}_trimmed_2.fastq"
    output:
        trim1 = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/{run}_fastp_1.fastq"),
        trim2 = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/{run}_fastp_2.fastq"),
        json = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/fastp.json",
        html = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/fastp.html",
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/{run}_fastp.log"
    singularity:
        "docker://biocontainers/fastp:v0.20.1_cv1"
    threads: 4
    resources:
        runtime = 240,
        mem_mb = 20000,
        slurm_partition = config["clusters"][config["parts"]]["short"],
        threads = 4,
    priority: 9999
    shell:
        """
        fastp -w {threads} -i {input.r1} -I {input.r2} -o {output.trim1} -O {output.trim2} -j {output.json} -h {output.html} -g 2> {log}
        """

rule UDSEQ_UBAM:
    ## Generates unaligned bam
    input:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/{run}_fastp_1.fastq",
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/fastq/{run}_fastp_2.fastq"
    output:
        ubam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/unaln.bam")
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/unaln.bam.log"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/unaln.benchmark"
    resources:
        runtime = 240,
        iotasks = 2,
        mem_mb = 4000,
        jv_mem = 3750,
        slurm_partition = config["clusters"][config["parts"]]["short"],
    params:
        TEMP_DIR = TEMP_DIR,
        TEMP_BAM_PATH = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis +  "/markdups/",
    priority: 9999
    shell:
        """
        gatk --java-options '-Xmx{resources.jv_mem}M' \
            FastqToSam \
            --FASTQ {input[0]} \
            --FASTQ2 {input[1]} \
            --OUTPUT {output.ubam} \
            --READ_GROUP_NAME {wildcards.run} \
            --SAMPLE_NAME {wildcards.sample} \
            --LIBRARY_NAME {wildcards.run} \
            --PLATFORM ILLUMINA > {log} 2>&1
        """

rule UDSEQ_MARK_ADAPTERS:
    ## Marks reads with adapter sequences
    input:
        ubam = rules.UDSEQ_UBAM.output
    output:
        obam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/unaln.adaptmarked.bam"),
        metrics = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/unaln.adaptmarked.metrics"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/unaln.adaptmarked.bam.log"
    params:
        tmp = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/markadapts/",
    benchmark: 
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/unaln.adaptmarked.benchmark"
    resources:
        runtime = 240,
        mem_mb = 2000,
        jv_mem = 1900,
        iotasks = 2,
        slurm_partition = config["clusters"][config["parts"]]["short"],
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    priority: 1000
    shell:
        """
        mkdir -p {params.tmp}
        gatk --java-options '-Xmx{resources.jv_mem}M' \
            MarkIlluminaAdapters \
            -I {input.ubam} \
            -O {output.obam} \
            -M {output.metrics} \
            -TMP_DIR {params.tmp} 2> {log}
        rm -rf {params.tmp}
        """

rule UDSEQ_MARKED_SAM_TO_FASTQ:
    input:
        bam = rules.UDSEQ_MARK_ADAPTERS.output.obam,
    output:
        ifastq = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/interleaved.fastq")
    log:
        samtofastq = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.samtofastq.log",
    params:
        align_threads = ALIGN_THREADS - 2,
        tmp = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/marktofq/",
    singularity: f"/gpfs/data/petljaklab/containers/gatk-alignment/gatk-alignment_{GATK_ALIGNER_VER}.sif"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.samtofastq.benchmark"
    resources:
        runtime = 480,
        mem_mb = 4000,
        jv_mem = 3750,
        iotasks = 2,
        slurm_partition = config["clusters"][config["parts"]]["short"],
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


rule UDSEQ_BWA:
    input:
        ifastq = rules.UDSEQ_MARKED_SAM_TO_FASTQ.output.ifastq,
        bwa_idxbase = lambda wildcards: ALN_REFERENCES[wildcards.reference],
    output:
        rawbam = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.raw.sam")
    threads: ALIGN_THREADS
    resources:
        threads = ALIGN_THREADS,
        cpus = ALIGN_THREADS,
        mem_mb = 30000,
        iotasks = 4,
        slurm_partition = config["clusters"][config["parts"]]["med"],
        tmpdisk = lambda wc, input: int(np.round(input.size_mb))
    log:
        bwa = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.bwa.log",
    singularity: f"/gpfs/data/petljaklab/containers/gatk-alignment/gatk-alignment_{GATK_ALIGNER_VER}.sif"
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.bwa.benchmark"
    params:
        tmpfastqpath = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/fastq/",
        tmpfastq = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/fastq/input.fastq",
    priority: 1002
    shell:
        """
            set +eo pipefail;
            mkdir -p {params.tmpfastqpath};
            cp {input.ifastq} {params.tmpfastq};
            bwa mem -C -M -p -t {threads} -R "@RG\\tID:{wildcards.run}\\tSM:{wildcards.sample}\\tPL:ILLUMINA" {input.bwa_idxbase} {params.tmpfastq} > {output} 2> {log.bwa};
            if [ $? -ne 0 ]; then
                rm {params.tmpfastq};
                set -eo pipefail;
                exit 69;
            fi;
            rm {params.tmpfastq}
        """

rule UDSEQ_SORT:
    input:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.raw.sam"
    output:
        temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.sorted.bam")
    threads: ALIGN_THREADS
    resources:
        threads = ALIGN_THREADS,
        cpus = ALIGN_THREADS,
        mem_mb = 30000,
        iotasks = 4,
        slurm_partition = config["clusters"][config["parts"]]["med"],
        tmpdisk = lambda wc, input: int(np.round(input.size_mb))
    singularity: f"/gpfs/data/petljaklab/containers/samtools/samtools_{SAMTOOLS_VERSION}.sif"
    log:
        temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.sorted.log")
    params:
        tmp_path = lambda wildcards: TEMP_DIR + wildcards.run + "/" + wildcards.analysis + "/sortbam/sorting",
    shell:
        "mkdir -p {params.tmp_path} && samtools sort -@ {threads} -T {params.tmp_path} {input} > {output} 2> {log} || rm -rf {params.tmp_path}; rm -rf {params.tmp_path}"

rule UDSEQ_MARKDUPS:
    input:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.sorted.bam"
    output:
        bam = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.bam",
        metrics = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.metrics"
    threads: 1
    resources:
        threads = 1,
        cpus = 1,
        mem_mb = 10000,
        runtime = 60*24*3,
        jv_mem = 8000,
        iotasks = 4,
        slurm_partition = config["clusters"][config["parts"]]["med"],
        tmpdisk = lambda wc, input: int(np.round(input.size_mb))
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_4.4.0.0.sif"
    log: 
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/UDSEQ_BAM/{analysis}/bam/{reference}/aligned.dupsflagged.log",
    shell:
        """
        gatk --java-options '-Xmx{resources.jv_mem}M' \
            MarkDuplicates \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics} \
            --ASSUME_SORT_ORDER coordinate \
            --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$" \
            --DUPLEX_UMI \
            --TAGGING_POLICY OpticalOnly 2> {log}
        """

rule UDSEQ_MERGE:
    ## Merge multi-run samples into single bams
    ## Currently merges ALL bam runs from a sample into a single bam. If we want subsets, will need to add that functionality later
    input:
        aggregate_udseq_runs
    output:
        merge = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/merged.bam")
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/merged.bam.log"
    singularity: f"/gpfs/data/petljaklab/containers/gatk/gatk_{MUTECT_VERSION}.sif"
    threads: 2
    resources:
        threads = 2,
        cpus = 2,
        runtime = 480,
        iotasks = 4,
        slurm_partition = config["clusters"][config["parts"]]["short"],
    params:
        inputlist = lambda wildcards, input: f"-I {input}" if isinstance(input, str) else "-I " + " -I ".join(input)
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/merged.benchmark"
    priority: 1004
    shell:
        """
        gatk MergeSamFiles \
            {params.inputlist} -O {output.merge} \
            --USE_THREADING true \
             &>> {log}
        """

rule UDSEQ_BAM2CRAM:
    input:
        merge = rules.UDSEQ_MERGE.output,
        bwa_idxbase = lambda wildcards: ALN_REFERENCES[wildcards.reference],
    output:
        cram = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/merged.cram",
        crai = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/merged.cram.crai",
        ref_md5 = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/merged.ref.md5",
        readme = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/CRAM_README.txt",
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/merged.cram.log"
    singularity: f"/gpfs/data/petljaklab/containers/samtools/samtools_{SAMTOOLS_VERSION}.sif"
    threads: 1
    resources:
        runtime = 480,
        slurm_partition = config["clusters"][config["parts"]]["short"],
        mem_mb = 4000,
    benchmark: SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/merged.cram.bench"
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

rule UDSEQ_CALL:
    input:
        cram = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/merged.cram",
        BWA_IDXBASE = "/gpfs/data/petljaklab/resources/{reference}/pipeline_resources/genome/udseq/{reference}.fa",
    output:
        variants = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/calls/calls_snv.vcf",
        tmp_cov_files = temp(directory("tmp/" + SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/calls/")),
    conda: f"{UDSEQ_BASE_PATH}/environment.yaml"
    params:
        out_prefix = lambda wildcards: SCRATCH_DIR + f"studies/{wildcards.study}/samples/{wildcards.sample}/analyses/UDSEQ_BAM/{wildcards.analysis}/merge/{wildcards.reference}/calls/calls",
        gnomad = lambda wildcards: f"/gpfs/data/petljaklab/resources/{wildcards.reference}/pipeline_resources/udseq/af-only-gnomad.4.0.0.hg38.vcf.gz",
        noise = lambda wildcards:  f"/gpfs/data/petljaklab/resources/{wildcards.reference}/pipeline_resources/udseq/NOISE.hg38.sorted.withchr.bed.gz",
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/calls/calls.log"
    threads: 12
    resources:
        runtime = 600,
        slurm_partition = config["clusters"][config["parts"]]["short"],
        mem_mb = 12*4000,
    shell:
        """
        DupCaller.py call -b {input.cram} -f {input.BWA_IDXBASE} -o {params.out_prefix} -p {threads} -g {params.gnomad} -m {params.noise} -maf 1 2> {log}
        """


rule UDSEQ_PAIRED_CALL:
    input:
        cram = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/merged.cram",
        normal = lambda wildcards: parent_cell(wildcards),
        BWA_IDXBASE = "/gpfs/data/petljaklab/resources/{reference}/pipeline_resources/genome/udseq/{reference}.fa",
    output:
        variants = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/calls_paired/calls_snv.vcf",
        tmp_cov_files = temp(directory("tmp/" + SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/calls_paired/")),
    conda: f"{UDSEQ_BASE_PATH}/environment.yaml"
    params:
        out_prefix = lambda wildcards: SCRATCH_DIR + f"studies/{wildcards.study}/samples/{wildcards.sample}/analyses/UDSEQ_BAM/{wildcards.analysis}/merge/{wildcards.reference}/calls_paired/calls",
        gnomad = lambda wildcards: f"/gpfs/data/petljaklab/resources/{wildcards.reference}/pipeline_resources/udseq/af-only-gnomad.4.0.0.hg38.vcf.gz",
        noise = lambda wildcards:  f"/gpfs/data/petljaklab/resources/{wildcards.reference}/pipeline_resources/udseq/NOISE.hg38.sorted.withchr.bed.gz",
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/calls_paired/calls.log"
    threads: 12
    resources:
        runtime = 600,
        slurm_partition = config["clusters"][config["parts"]]["short"],
        mem_mb = 12*4000,
    shell:
        """
        DupCaller.py call -b {input.cram} -n {input.normal} -f {input.BWA_IDXBASE} -o {params.out_prefix} -p {threads} -g {params.gnomad} -m {params.noise} -maf 1 2> {log}
        """

rule UDSEQ_VCFTOTABLE:
    input:
        variants = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/{calls}/calls_snv.vcf"
    output:
        variants = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/{calls}/calls_snv.tbl"
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{MUTECT_VERSION}.sif"
    threads: 1
    resources:
        runtime = 10,
        slurm_partition = config["clusters"][config["parts"]]["dev"],
    shell:
        "gatk VariantsToTable -V {input} -O {output} -raw"


## Function that gets all the samples for the study/cell line
def udseq_all_samples(wildcards, ext):
    study_id = petljakapi.translate.stringtoid(wildcards.study)
    sample_id = petljakapi.translate.stringtoid(wildcards.sample)
    entries = petljakapi.select.multi_select(db, "samples", {"study_id":study_id})
    samp = [entry for entry in entries if entry[0] == sample_id][0]
    relevant_samples = [entry for entry in entries if entry[6] == samp[6]]
    rel_samp_ids = [petljakapi.translate.idtostring(entry[0], "MPS") for entry in relevant_samples]
    analyses = [petljakapi.select.multi_select(db, "analyses", {"samples_id":entry[0], "pipeline_name":"UDSEQ_BAM"}) for entry in relevant_samples]
    analyses = [petljakapi.translate.idtostring(a[0][0], "MPA") for a in analyses if a]
    analysis_paths = [f"{SCRATCH_DIR}studies/{wildcards.study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{wildcards.reference}/{wildcards.calls}/{ext}" for analysis, sample in zip(analyses, rel_samp_ids)]
    return(analysis_paths)


rule UDSEQ_PREP_GENOTYPING:
    input:
        lambda w: udseq_all_samples(w, "calls_snv.tbl")
    output:
        variants = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/{calls}/calls_agg.txt"
    threads: 1
    resources:
        runtime = 10,
        slurm_partition = config["clusters"][config["parts"]]["dev"],
    shell:
        "cat {input} | cut -f1,2 | grep -v CHROM | sort -k1,1 -k2,2n | uniq > {output}"

rule UDSEQ_GENOTYPE:
    input:
        cram = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/merged.cram",
        positions = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/{calls}/calls_agg.txt",
        bwa_idxbase = lambda wildcards: ALN_REFERENCES[wildcards.reference],
    output:
        variants = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/{calls}/variant_genotyping.txt"
    singularity: "/gpfs/data/petljaklab/containers/bcftools/bcftools_latest.sif"
    threads: 1
    resources:
        runtime = 60,
        slurm_partition = config["clusters"][config["parts"]]["dev"],
    shell:
        """
        bcftools mpileup -R {input.positions} \
            --fasta-ref {input.bwa_idxbase} \
            {input.cram} | \
            bcftools call -c | \
            grep -v '#' | \
            awk '$5 != "."{{print $0}}' > {output}
        """

rule UDSEQ_AGGREGATE:
    input:
        lambda w: udseq_all_samples(w, "variant_genotyping.txt")
    output:
        variants = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/{calls}/genotyped_variant_positions.txt"
    threads: 1
    resources:
        runtime = 10,
        slurm_partition = config["clusters"][config["parts"]]["dev"],
    shell:
        "cat {input} | cut -f1,2 | sort -k1,1 -k2,2n > {output}"

rule UDSEQ_FILTER_SHAREDNESS:
    input:
        calls = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/{calls}/calls_snv.tbl",
        genos = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/{calls}/genotyped_variant_positions.txt"
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/{calls}/calls_snv_private.txt"
    threads: 1
    resources:
        runtime = 10,
        slurm_partition = config["clusters"][config["parts"]]["dev"],
    shell:
        """
        module load r/4.1.2; Rscript modules/UDSeq/scripts/filter_private.R -v {input.calls} -g {input.genos} -o {output}
        """

rule UDSEQ_COV_FILE:
    input:
        "tmp/" + SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/{calls}/"
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/UDSEQ_BAM/{analysis}/merge/{reference}/{calls}/calls_coverage_stats.txt"
    threads: 1
    resources:
        runtime = 30,
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        mem_mb = 4000,
    shell:
        "python modules/UDSeq/scripts/compute_mol_cov.py -i {input} -o {output}"