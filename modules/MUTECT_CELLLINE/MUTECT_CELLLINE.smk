def parent_cell(wildcards):
    sample_id = petljakapi.translate.stringtoid(wildcards.sample)
    db_line = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = sample_id)
    parent_id = petljakapi.translate.idtostring(db_line[0][5], "MPS")
    parent_merge = gateway("WGS_MERGE_BAM", parent_id, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb_devel")[0]
    return(parent_merge)


rule MUTECT2:
    input:
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb_devel")[0],
        parent_merge = lambda wildcards: parent_cell(wildcards)
    output:
        PROD_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/variants.vcf",
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    resources:
        threads = 4,
        mem_mb = 50000
    benchmark:
        PROD_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/variants.resources",
    log:
        "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/variants.log"
    shell:
        """
            gatk Mutect2 -R {input.fa} \
            -I {input.cell_merge} \
            -I {input.parent_merge} \
            --normal {wildcards.sample} \
            --panel-of-normals /gpfs/data/petljaklab/resources/hg19/pipeline_resources/somatic_celline/reference_vcf/pon.vcf \
            -O {output} &> {log}
        """