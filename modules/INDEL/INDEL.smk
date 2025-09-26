import re

INDEL_PATH = os.path.join(basedir, "modules", "INDEL")

def mutect_std_calls(wildcards):
    ## Run the gateway function to ensure the mutect run has a path
    final_mutect_path = gateway("SBS", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)
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
        slurm_partition = "petljaklab,cpu_dev",
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
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)[0],
        parent_merge = lambda wildcards: parent_cell(wildcards),
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
    output:
        tmpfile = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/varscan2/output.indel.vcf"),
        #realout = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/varscan2/indels.vcf"
    threads: 1
    resources:
        threads = 1,
        cpus = 1,
        runtime = 60*24*5,
        slurm_partition = "petljaklab,cpu_medium",
        mem_mb = 2000
    singularity: "/gpfs/data/petljaklab/containers/varscan2/varscan2_2.4.6.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/varscan2/indels.log"
    params:
        outpath = lambda wildcards: SCRATCH_DIR + f"studies/{wildcards.study}/samples/{wildcards.sample}/analyses/INDEL/{wildcards.analysis}/{wildcards.reference}/varscan2/output",
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
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)[0],
        parent_merge = lambda wildcards: parent_cell(wildcards),
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/strelka2/runWorkflow.py"
    singularity: "/gpfs/data/petljaklab/containers/strelka2/strelka2_2.9.10.sif"
    resources:
        runtime = 10,
        slurm_partition = "petljaklab,cpu_dev",
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    params:
        odir = lambda wildcards: SCRATCH_DIR + f"studies/{wildcards.study}/samples/{wildcards.sample}/analyses/INDEL/{wildcards.analysis}/{wildcards.reference}/strelka2/"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/strelka_log.txt"
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
        runtime = 60*12,
        slurm_partition = "petljaklab,cpu_short",
        mem_mb =  lambda wildcards, attempt: 20000 * attempt
    singularity: 
        "/gpfs/data/petljaklab/containers/strelka2/strelka2_2.9.10.sif"
    benchmark: 
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/strelka2/benchmark.txt"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/strelka2/log.txt"
    params:
        workspacedir = lambda wildcards: SCRATCH_DIR + f"studies/{wildcards.study}/samples/{wildcards.sample}/analyses/INDEL/{wildcards.analysis}/{wildcards.reference}/strelka2/workspace/"
    shell:
        """
        rm -rf {params.workspacedir};
        {input} -j {threads} -m local 2>> {log};
        gunzip -c {output.tmpout} > {output.realout} 2>> {log}
        """

rule varscan2_split:
## TODO: change the parent/cell stuff to something that works for either cells or patients
    input:
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)[0],
        parent_merge = lambda wildcards: parent_cell(wildcards),
        fa = lambda wildcards: FA_PATHS[wildcards.reference],
    output:
        tmpfile = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/varscan2/{chrom}/output.indel.vcf",
    threads: 1
    resources:
        threads = 1,
        cpus = 1,
        runtime = 60*12,
        slurm_partition = "petljaklab,cpu_short",
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    singularity: "/gpfs/data/petljaklab/containers/varscan2/varscan2_2.4.6.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/varscan2/{chrom}/indels.log"
    params:
        outpath = lambda wildcards: SCRATCH_DIR + f"studies/{wildcards.study}/samples/{wildcards.sample}/analyses/INDEL/{wildcards.analysis}/{wildcards.reference}/varscan2/{wildcards.chrom}/output",
    shell:
        """
        java -jar /varscan2/VarScan.v2.4.6.jar somatic <(samtools mpileup -r {wildcards.chrom} -q 10 -f {input.fa} {input.parent_merge}) \
                                                       <(samtools mpileup -r {wildcards.chrom} -q 10 -f {input.fa} {input.cell_merge}) \
                                                       {params.outpath} \
                                                       --output-vcf \
                                                       2> {log}
        """

rule merge_varscan_runs:
    input:
        lambda wildcards: expand(SCRATCH_DIR + "studies/{{study}}/samples/{{sample}}/analyses/INDEL/{{analysis}}/{{reference}}/varscan2/{chrom}/output.indel.vcf", chrom = chromosomes[wildcards.reference])
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/varscan2/indels.vcf"
    threads: 1
    resources:
        threads = 1,
        cpus = 1,
        runtime = 60,
        slurm_partition = "petljaklab,cpu_dev",
        mem_mb = lambda wildcards, attempt: 4000 * attempt,
    log: 
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/varscan2/indels.log"
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    params:
        dic = lambda wildcards: re.sub("\\.fa$", ".dict", FA_PATHS[wildcards.reference]),
        inputlist = lambda wildcards, input: " -I ".join([input]) if isinstance(input, str) else " -I ".join(input)
    shell:
        "gatk MergeVcfs -D {params.dic} -I {params.inputlist} -O {output} &> {log}"

rule ID_VCF_TO_TABLE:
    input:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/{tool}/indels.vcf",
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/{tool}/indels.txt",
    threads: 1
    resources:
        threads = 1,
        cpus = 1,
        runtime = 5,
        slurm_partition = "petljaklab,cpu_dev",
        mem_mb = lambda wildcards, attempt: 4000 * attempt,
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    log: 
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/{tool}/indels.log"
    shell:
        "gatk VariantsToTable -V {input} -O {output} -raw &>> {log}"



rule ID_PARENT_DAUGHTER_TABLES:
    output:
        daughter = SCRATCH_DIR + "studies/{study}/analyses/INDEL/{reference}/daughters_table.txt",
        #parental = SCRATCH_DIR + "studies/{study}/analyses/INDEL/{reference}/parents_table.txt",
    threads: 1
    resources:
        slurm_partition = "petljaklab,cpu_dev",
        cpus = 1,
        threads = 1,
        mem_mb = 2000,
        runtime = 15,
    run:
        if not os.path.exists(f"{SCRATCH_DIR}studies/{wildcards.study}/analyses/INDEL/{wildcards.reference}"):
            os.makedirs(f"{SCRATCH_DIR}studies/{wildcards.study}/analyses/INDEL/{wildcards.reference}")
        study_id = petljakapi.translate.stringtoid(wildcards.study)
        samples_ids = petljakapi.select.multi_select(db, "samples", {"study_id":study_id})
        for sample in samples_ids:
            cell_name = petljakapi.select.multi_select(db, "cells", {"id":sample[6]})[0][1]
            #files = gateway("MUTECT", petljakapi.translate.idtostring(sample[0], "MPS"), scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = db)
            #files = [file for file in files if not file.contains("line_of_origin.txt")]
            ## Get analysis row
            analysis = petljakapi.select.multi_select(db = db, table = "analyses", filters = {"samples_id":sample[0], "pipeline_name":"INDEL", "reference_genome":wildcards.reference}, bench = config["bench"])
            if len(analysis) == 0:
                continue
            analysis = analysis[0]
            aid = petljakapi.translate.idtostring(analysis[0], "MPA")
            ap = analysis[10]
            ## Test if this sample is used as a reference for any other samples
            daughters = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "sample_parent_id", filter_value = sample[0], bench = config["bench"])
            parent_id = sample[5]
            if parent_id is not None and parent_id != "NULL": ## if this sample has a parent ie. it's a daughter
                files = [f"{wildcards.reference}/mutect2/indels.txt"]
                files = [os.path.join(ap, aid, f) for f in files]
                parent_name = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = parent_id)[0][1]
                print(files, sample)
                with open(output.daughter, "a") as f:
                    f.write("%s\n" % '\t'.join([files[0], sample[1], parent_name, cell_name, sample[4], str(sample[7])]))
            #if len(daughters) > 0 or parent_id == "NULL" or parent_id is None: ## ie. if this is a parent
            #    files = [os.path.join(ap, aid, wildcards.reference, "mutect2", "indels.txt")]
            #    with open(output.parental, "a") as f:
            #        f.write("%s\n" % '\t'.join([files[0], sample[1], cell_name, sample[4]]))

rule INDEL_CELLLINE_IMPORT_DB_TABLE:
    input:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/{tool}/indels.vcf",
    output:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/INDEL/{analysis}/{reference}/{tool}/indels.vcf.dbtmp",
    resources:
        runtime = 10,
        threads = 1,
        mem_mb = 5000,
        slurm_partition = config["clusters"][config["parts"]]["short"]
    priority: 1009
    params:
        script_path = INDEL_PATH + "/scripts/load_mutations_table.py",
        db = db,
    shell:
        """
        python {params.script_path} -d {params.db} \
                                    -g {wildcards.reference} \
                                    -a {wildcards.analysis} \
                                    -s {wildcards.sample} \
                                    -u {wildcards.study} \
                                    -t ID \
                                    -v {input.vcf}
        """