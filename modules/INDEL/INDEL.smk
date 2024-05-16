import re

def mutect_std_calls(wildcards):
    ## Run the gateway function to ensure the mutect run has a path
    final_mutect_path = gateway("MUTECT", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb")
    ## Get the analysis ID
    mutect_path = os.path.join(re.sub(r"(.*/MPA[0-9]*/mutect2/[^/]*/).*", "\\1", final_mutect_path[0]), "std/filtered.vcf")
    if isinstance(mutect_path, list):
        raise ValueError("mutect_std_calls tried to return a list! Something is wrong and you should debug this.")
    return(mutect_path)

rule copy_m2_calls:
## We already have code to run mutect2 for SBS, so just reuse those results
    input:
        mutect_std_calls
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/mutect2/indels.vcf"
    threads: 1
    resources:
        slurm_partition = "cpu_dev",
        runtime = 15,
        cpus = 1,
        threads = 1,
        mem_mb = 2000,
    log: 
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/mutect2/indels.log"
    shell:
        """
        ## First copy over the header
        grep '#' {input} > {output} 2> {log}
        ## Now use awk to get indels
        grep -v '#' {input} | awk -v OFS='\\t' -v FS='\\t' 'length($4) != 1 || length($5) != 1{{print $0}}' >> {output} 2>> {log}
        echo 'Mutect2 calls copied from {input}'
        """

rule varscan2:
## TODO: change the parent/cell stuff to something that works for either cells or patients
    input:
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb")[0],
        parent_merge = lambda wildcards: parent_cell(wildcards),
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
    output:
        tmpfile = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/varscan2/output.indel.vcf"),
        realout = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/varscan2/indels.vcf"
    threads: 1
    resources:
        threads = 1,
        cpus = 1,
        runtime = 60*24*5,
        slurm_partition = "cpu_medium",
        mem_mb = 2000
    singularity: "/gpfs/data/petljaklab/containers/varscan2/varscan2_2.4.6.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/varscan2/indels.log"
    params:
        outpath = lambda wildcards: SCRATCH_DIR + "studies/{wildcards.study}/samples/{wildcards.sample}/analyses/INDEL/{wildcards.analysis}/{wildcards.reference}/varscan2/output",
    shell:
        """
        java -jar /varscan2/VarScan.v2.4.6.jar somatic <(samtools mpileup -q 10 -f {input.fa} {input.parent_merge}) \
                                                       <(samtools mpileup -q 10 -f {input.fa} {input.cell_merge}) \
                                                       {params.outpath} \
                                                       --output-vcf \
                                                       2> {log}
        cp {output.tmpfile} {output.realout} 2>> {log}
        """

rule makeStrelkaWorkflow:
    input:
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb")[0],
        parent_merge = lambda wildcards: parent_cell(wildcards),
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/strelka2/runWorkflow.py"
    singularity: "/gpfs/data/petljaklab/containers/strelka2/strelka2_2.9.10.sif"
    resources:
        runtime = 10,
        slurm_partition = "cpu_dev",
        mem_mb = 8000
    params:
        odir = lambda wildcards: SCRATCH_DIR + "studies/{wildcards.study}/samples/{wildcards.sample}/analyses/INDEL/{wildcards.analysis}/{wildcards.reference}/strelka2/"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/strelka2/log.txt"
    shell:
        """
        rm -rf {params.odir}
        configureStrelkaSomaticWorkflow.py \
            --normalBam={input.parent_merge} \
            --tumorBam={input.cell_merge} \
            --referenceFasta={input.fa} \
            --runDir={params.odir} 2> {log}
        """

rule runStrelka2:
    input:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/strelka2/runWorkflow.py"
    output:
        tmpout = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/strelka2/results/variants/somatic.indels.vcf.gz"),
        realout = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/strelka2/indels.vcf",
    threads: 12
    resources:
        threads = 12,
        cpus = 12,
        runtime = 60*6,
        slurm_partition = "cpu_medium",
        mem_mb = 10000
    singularity: 
        "/gpfs/data/petljaklab/containers/strelka2/strelka2_2.9.10.sif"
    benchmark: 
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/strelka2/benchmark.txt"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/strelka2/log.txt"
    shell:
        """
        ./{input} -j {threads} -m local 2>> {log};
        gunzip -c {output.tmpout} > {output.realout} 2>> {log}
        """