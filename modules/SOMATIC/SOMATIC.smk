rule GET_SBS:
    input:
        lambda wildcards: gateway(analysis_name="SBS", given_id=wildcards.sample, scratch_dir=SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)
    output:
        vcf = PROD_DIR + "studies/{study}/samples/{sample}/analyses/SOMATIC/{analysis}/mutations/{reference}/small_mutations/sbs.vcf",
    threads: 1
    resources:
        mem_mb = 500,
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        runtime = 240,
    shell:
        "ln -s --relative {input[0]} {output}"

rule GET_ID:
    input:
        lambda wildcards: gateway(analysis_name="INDEL", given_id=wildcards.sample, scratch_dir=SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)
    output:
        vcf = PROD_DIR + "studies/{study}/samples/{sample}/analyses/SOMATIC/{analysis}/mutations/{reference}/small_mutations/indel.vcf",
    threads: 1
    resources:
        mem_mb = 500,
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        runtime = 240,
    shell:
        "ln -s --relative {input[0]} {output}"