SRA_PIPE_VERSION = "1.0.0"

rule PREFETCH:
    output:
        sra = temp(SCRATCH_DIR + "data/sra_" + SRA_PIPE_VERSION + "/{proj}/{samp}/{run}/{run}.sra"),
    retries: 3
    params:
        SCRATCH_DIR = SCRATCH_DIR,
        SRA_PIPE_VERSION = SRA_PIPE_VERSION
    log:
        SCRATCH_DIR + "data/sra_" + SRA_PIPE_VERSION + "/{proj}/{samp}/{run}/{run}_prefetch.log"
    shell:
        "prefetch {wildcards.run} -O {params.SCRATCH_DIR}data/sra_{params.SRA_PIPE_VERSION}/{wildcards.proj}/{wildcards.samp}/ -f ALL -r no &> {log};"

rule DUMP:
    input:
        sra = SCRATCH_DIR + "data/sra_" + SRA_PIPE_VERSION + "/{proj}/{samp}/{run}/{run}.sra",
    output:
        temp(SCRATCH_DIR + "data/sra_" + SRA_PIPE_VERSION + "/{proj}/{samp}/{run}_1.fastq"),
        temp(SCRATCH_DIR + "data/sra_" + SRA_PIPE_VERSION + "/{proj}/{samp}/{run}_2.fastq")
    params:
        SCRATCH_DIR = SCRATCH_DIR,
        SRA_PIPE_VERSION = SRA_PIPE_VERSION
    log:
        SCRATCH_DIR + "data/sra_" + SRA_PIPE_VERSION + "/{proj}/{samp}/{run}_dump.log"
    shell:
        "fasterq-dump {input.sra} --split-files -t temp/{wildcards.run} -O {params.SCRATCH_DIR}data/sra_{params.SRA_PIPE_VERSION}/{wildcards.proj}/{wildcards.samp}/ &> {log}; rm -rf ./{wildcards.run}/"

rule SRA_DB:
    input:
        SCRATCH_DIR + "data/sra_" + SRA_PIPE_VERSION + "/{proj}/{samp}/{run}_1.fastq",
        SCRATCH_DIR + "data/sra_" + SRA_PIPE_VERSION + "/{proj}/{samp}/{run}_2.fastq",
    output:
        SCRATCH_DIR + "data/sra_" + SRA_PIPE_VERSION + "/{proj}/{samp}/{run}.stamp.SRAPIPE