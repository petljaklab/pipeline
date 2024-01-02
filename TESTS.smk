## Fully synthetic testing data based on the TP53 locus of the human genome. 
## 
## We generate 40x coverage of the TP53 locus with some SNPs to test alignment/processing pipelines

TESTS_PIPE_VERSION = "1.0.0"

## Cut out the tp53 locus & index
## Recall TP53 locus is at position chr17:7,571,739-7,590,808 on hg19
rule CUTFA:
    input:
        ALN_REFERENCES["hg19"]
    output:
        temp(SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fa/cut.fa")
    singularity:
        "/gpfs/data/petljaklab/containers/bedtools/bedtools_2.31.0.sif"
    shell:
        "bedtools getfasta -fi {input} -bed <(echo '17\t7571739\t7590808') > {output}"

rule INDFA:
    input:
        SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fa/cut.fa"
    output:
        temp(SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fa/cut.fa.fai")
    singularity:
        "/gpfs/data/petljaklab/containers/samtools/samtools_1.18.sif"
    shell:
        "samtools faidx {input}"
## Simulate 40x coverage of this gene
rule GEN_READS:
    input:
        fa = SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fa/cut.fa",
        fai = SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fa/cut.fa.fai",
    output:
        reads1 = SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fq/readstmp1.fq",
        reads2 = SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fq/readstmp2.fq",
    singularity:
        "/gpfs/data/petljaklab/containers/art/art_latest.sif"
    params:
        SCRATCH_DIR = SCRATCH_DIR,
        outbase = SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fq/readstmp"
    log:
        SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fq/art.log"
    shell:
        """
        art_illumina -ss HS25 -sam -i {input.fa} -p -l 150 -f 40 -m 200 -s 10 -o {params.outbase} 2> {log}; \
        """

rule MARK_COMPLETE:
    input:
        reads1 = SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fq/readstmp1.fq",
        reads2 = SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fq/readstmp2.fq",
    output:
        reads1 = SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fq/reads1.fq",
        reads2 = SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fq/reads2.fq",
    log:
        SCRATCH_DIR + "studies/{study_id}/samples/{sample_id}/runs/{run_id}/analyses/TESTS/{aid}/fq/api_completion.txt",
    shell:
        "python scripts/mark_complete.py -i {wildcards.aid} -d petljakdb_devel {output} > {log} 2>&1; mv {input.reads1} {output.reads1}; mv {input.reads2} {output.reads2}"