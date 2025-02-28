rule PILEUP_SUMMARIES:
    input:
        BAM = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)[0],
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
    output:
        SUMMARY = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/pileup.txt",
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/pileup.log",
    benchmark:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/pileup.benchmark",
    resources:
        mem_mb = 3000,
        jv_mem = 2900,
        slurm_partition = config["clusters"][config["parts"]]["short"],
        runtime = 60 * 4,
    shell:
        """
        gatk --java-options '-Xmx{resources.jv_mem}M' \
            GetPileupSummaries \
            -I {input.BAM} \
            -R {input.fa} \
            -V /gpfs/data/petljaklab/resources/{wildcards.reference}/pipeline_resources/somatic_celline/reference_vcf/small_exac_common_3.vcf \
            -L /gpfs/data/petljaklab/resources/{wildcards.reference}/pipeline_resources/somatic_celline/reference_vcf/small_exac_common_3.vcf \
            -O {output.SUMMARY} &> {log}
        """

rule COMBINE_MUTECT2_VCFS:
    input:
        FILES = lambda wildcards: expand(SCRATCH_DIR + "studies/{{study}}/samples/{{sample}}/analyses/{{pipe}}/{{analysis}}/mutect2/{{reference}}/{{type}}/raw_{chrom}.vcf", chrom = chromosomes[wildcards.reference])
    output:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/{type}/merged.vcf"
    params:
        inputlist = lambda wildcards, input: " -I ".join([input]) if isinstance(input, str) else " -I ".join(input.FILES)
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/{type}/merge_vcf.log"
    resources:
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        jv_mem = 2900,
        runtime = 240,
        mem_mb = 3000
    shell:
        "gatk --java-options '-Xmx{resources.jv_mem}M' MergeVcfs -I {params.inputlist} -O {output.vcf} &> {log}"

rule MERGE_MUTECT2_STATS:
    input:
        FILES = lambda wildcards: expand(SCRATCH_DIR + "studies/{{study}}/samples/{{sample}}/analyses/{{pipe}}/{{analysis}}/mutect2/{{reference}}/{{type}}/raw_{chrom}.vcf.stats", chrom = chromosomes[wildcards.reference])
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/{type}/merged.vcf.stats",
    params:
        inputlist = lambda wildcards, input: " --stats ".join([input]) if isinstance(input, str) else " --stats ".join(input.FILES)
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/{type}/stats.log",
    resources:
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        runtime = 240,
        jv_mem = 2900,
        mem_mb = 3000
    shell:
        """
        gatk --java-options '-Xmx{resources.jv_mem}M' \
            MergeMutectStats \
            -stats {params.inputlist} -O {output}
        """


rule MUTECT2_READ_ORIENTATION:
    input:
        files = lambda wildcards: expand(SCRATCH_DIR + "studies/{{study}}/samples/{{sample}}/analyses/{{pipe}}/{{analysis}}/mutect2/{{reference}}/{{type}}/f1r2_{chrom}.tar.gz", chrom = chromosomes[wildcards.reference])
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/{type}/read_orientation_model.tar.gz",
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/{type}/read_orientation_model.log",
    resources:
        mem_mb = lambda wildcards, attempt: 5000 * attempt,
        jv_mem = lambda wildcards, attempt: 4500 * attempt,
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        runtime = 30,
    params:
        inputlist = lambda wildcards, input: " -I ".join([input]) if isinstance(input, str) else " -I ".join(input.files)
    shell:
        "gatk --java-options '-Xmx{resources.jv_mem}M' LearnReadOrientationModel -I {params.inputlist} -O {output} --verbosity DEBUG &> {log}"


rule MUTECT2_CALCULATE_CONTAMINATION:
    input:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/pileup.txt"
    output:
        contam = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/contamination.table",
        segments = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/segments.table",
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/contamination.log"
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    resources:
        mem_mb = 3000,
        jv_mem = 2900,
        slurm_partition = config["clusters"][config["parts"]]["short"],
        runtime = 60 * 4,
    shell:
        """
        gatk --java-options '-Xmx{resources.jv_mem}M' \
            CalculateContamination -I {input} \
            -O {output.contam} --tumor-segmentation {output.segments} &> {log}
        """    


rule MUTECT2_FILTERMUTECTCALLS:
    input:
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
        raw = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/{type}/merged.vcf",
        contam = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/contamination.table",
        segments = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/segments.table",
        mergedstats = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/{type}/merged.vcf.stats",
        orientation = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/{type}/read_orientation_model.tar.gz",
    output:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/{type}/filtered.vcf",
        stats = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/{type}/filter.stats",
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/{pipe}/{analysis}/mutect2/{reference}/{type}/filtered.log"
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * (1 + ((attempt-1)/2)),
        jv_mem = lambda wildcards, attempt: 3500 * (1 + ((attempt-1)/2)),
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        runtime = 240,
    shell:
        """
        gatk --java-options '-Xmx{resources.jv_mem}M' \
            FilterMutectCalls -V {input.raw} \
            -R {input.fa} \
            -O {output.vcf} \
            --contamination-table {input.contam} \
            --tumor-segmentation {input.segments} \
            --ob-priors {input.orientation} \
            --filtering-stats {output.stats} &> {log}
        """