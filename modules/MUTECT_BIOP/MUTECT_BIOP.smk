import petljakapi.select
import petljakapi.translate

def matched_norm(wildcards):
    ## Sample ID
    samid = petljakapi.translate.stringtoid(wildcards.sample)
    patient_id = petljakapi.select.multi_select(db, "samples", {"id":samid})[0][8]
    germ_sample = petljakapi.select.multi_select(db, "patients", {"id":patient_id})[0][2]
    return(petljakapi.translate.idtostring(germ_sample, "MPS"))


rule MUTECT2_BIOP:
    input:
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
        biop_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)[0],
        germ_merge = lambda wildcards: gateway("WGS_MERGE_BAM", matched_norm(wildcards), scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)[0]
    output:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop/raw_{chrom}.vcf",
        stats= SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop/raw_{chrom}.vcf.stats",
        tgz = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop/f1r2_{chrom}.tar.gz"
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    resources:
        threads = 1,
        mem_mb = lambda wildcards, attempt: 4500 * (1 + ((attempt-1)/2)),
        jv_mem = lambda wildcards, attempt: 4000 * (1 + ((attempt-1)/2)),
        iotasks = 2,
        runtime = 24*60,
        slurm_partition = "petljaklab,cpu_medium",
        att = lambda wildcards, attempt: attempt
    benchmark:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop/variants_{chrom}.resources",
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop/variants_{chrom}.log"
    shell:
        """
            NORMAL_NAME=$(gatk GetSampleName -I {input.germ_merge} -R {input.fa} -O /dev/stdout 2> {log}); \
            gatk --java-options '-Xmx{resources.jv_mem}M' \
                Mutect2 -R {input.fa} \
                -L {wildcards.chrom} \
                -I {input.biop_merge} \
                -I {input.germ_merge} \
                --normal $NORMAL_NAME \
                --panel-of-normals /gpfs/data/petljaklab/resources/hg19/pipeline_resources/somatic_celline/reference_vcf/pon.vcf \
                --germline-resource /gpfs/data/petljaklab/resources/hg19/pipeline_resources/somatic_celline/reference_vcf/gnomad.vcf \
                --native-pair-hmm-threads 1 \
                --f1r2-tar-gz {output.tgz} \
                --output {output.vcf} &>> {log}.{resources.att}
        """

rule MUTECT2_BIOP_NONORM:
    input:
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
        biop_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)[0],
    output:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop_noref/raw_{chrom}.vcf",
        stats= SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop_noref/raw_{chrom}.vcf.stats",
        tgz = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop_noref/f1r2_{chrom}.tar.gz"
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    resources:
        threads = 1,
        mem_mb = lambda wildcards, attempt: 4500 * (1 + ((attempt-1)/2)),
        jv_mem = lambda wildcards, attempt: 4000 * (1 + ((attempt-1)/2)),
        iotasks = 2,
        runtime = 24*60,
        slurm_partition = "petljaklab,cpu_medium",
        att = lambda wildcards, attempt: attempt
    benchmark:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop_noref/variants_{chrom}.resources",
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop_noref/variants_{chrom}.log"
    shell:
        """
            gatk --java-options '-Xmx{resources.jv_mem}M' \
                Mutect2 -R {input.fa} \
                -L {wildcards.chrom} \
                -I {input.biop_merge} \
                --panel-of-normals /gpfs/data/petljaklab/resources/hg19/pipeline_resources/somatic_celline/reference_vcf/pon.vcf \
                --germline-resource /gpfs/data/petljaklab/resources/hg19/pipeline_resources/somatic_celline/reference_vcf/gnomad.vcf \
                --native-pair-hmm-threads 1 \
                --f1r2-tar-gz {output.tgz} \
                --output {output.vcf} &>> {log}.{resources.att}
        """

rule MUTECT2_BIOP_DONE:
    input:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop/filtered.vcf"
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop/filtered.done"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_BIOP/{analysis}/mutect2/{reference}/biop/filtered.done.log"
    params:
        db = config["db"]
    resources:
        runtime = 10,
        slurm_partition = config["clusters"][config["parts"]]["short"],
    priority: 1007
    shell:
        "python scripts/mark_complete.py --id {wildcards.analysis} --db {params.db} {input} >> {log}; touch {output}"