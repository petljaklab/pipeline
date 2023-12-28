## Fully synthetic testing data based on the TP53 locus of the human genome. 
## 
## We generate 40x coverage of the TP53 locus with some SNPs to test alignment/processing pipelines
## Recall TP53 locus is at position chr17:7,571,739-7,590,808 on hg19
TESTS_PIPE_VERSION = "1.0.0"

## Cut out the tp53 locus & index
rule CUTFA:
    input:
        ALN_REFERENCES["hg19"]
    output:
        temp(SCRATCH_DIR + "pipetest/fa/cut.fa")
    singularity:
        "/gpfs/data/petljaklab/containers/bedtools/bedtools_2.31.0.sif"
    shell:
        "bedtools getfasta -fi {input} -bed <(echo '17\t7571739\t7590808') > {output}"

rule INDFA:
    input:
        SCRATCH_DIR + "pipetest/fa/cut.fa"
    output:
        temp(SCRATCH_DIR + "pipetest/fa/cut.fa.fai")
    singularity:
        "/gpfs/data/petljaklab/containers/samtools/samtools_1.18.sif"
    shell:
        "samtools faidx {input}"
## Simulate 40x coverage of this gene
rule gen_reads:
    input:
        fa = SCRATCH_DIR + "pipetest/fa/cut.fa",
        fai = SCRATCH_DIR + "pipetest/fa/cut.fa.fai",
    output:
        reads1 = SCRATCH_DIR + "pipetest/fq/reads1.fq",
        reads2 = SCRATCH_DIR + "pipetest/fq/reads2.fq",
    singularity:
        "/gpfs/data/petljaklab/containers/art/art_latest.sif"
    params:
        SCRATCH_DIR = SCRATCH_DIR
    log:
        SCRATCH_DIR + "pipetest/fq/art.log"
    shell:
        "art_illumina -ss HS25 -sam -i {input.fa} -p -l 150 -f 40 -m 200 -s 10 -o {params.SCRATCH_DIR}pipetest/fq/reads 2> {log}"

rule TEST_DB_ENTRY:
    input:
        reads1 = SCRATCH_DIR + "pipetest/fq/reads1.fq",
        reads2 = SCRATCH_DIR + "pipetest/fq/reads2.fq",
    output:
        SCRATCH_DIR + "pipetest/fq/reads.stamp.TESTPIPE"
    log:
        SCRATCH_DIR + "pipetest/fq/db_log.txt"
    params:
        script_dir = os.path.join(workflow.basedir, "modules/TESTS/scripts/db.py")
    shell:
        "python {params.script_dir} -s TESTCASES -m TP53TST1 -r TP53TST1_1 -d petljakdb_devel -p {input.reads1} -o {output} > {log} 2>&1"