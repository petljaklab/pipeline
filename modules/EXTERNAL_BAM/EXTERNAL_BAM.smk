import logging
import re
import shutil
import os
EXTERNAL_BAM_VERSION = "1.0.0"
GATK_ALIGNER_VER = "1.0.1"

def get_external_bam_path(wildcards):
    if wildcards.study == "MPP000003":
        basepath = "/gpfs/data/petljaklab/data/external_data/raw_sequencing_archive/MPP000003/broad_data/bams/Bams/"
        iden = petljakapi.translate.stringtoid(wildcards.sample)
        res = petljakapi.select.multi_select(db = db, table = "samples", filters = {"id":iden}, bench = config["bench"])[0]
        rname = res[1]
        path = os.path.join(basepath, rname + ".bam")
    return(path)

def collate_output_fastqs(wildcards):
    out_path = checkpoints.SPLIT_BAM.get(**wildcards).output[0]
    newpath = re.sub("split_bams", "splitfq", out_path)
    return expand(newpath + "/{i}_1.fq",
                  i=glob_wildcards(os.path.join(out_path, "{i}.bam")).i)

def run_rname_to_rid(wildcards):
    sample_id = petljakapi.translate.stringtoid(wildcards.sample)
    selection = petljakapi.select.multi_select(db = db, table = "runs", filters = {"sample_id":sample_id, "rname":wildcards.run}, bench = config["bench"])
    return(selection[0][0])



rule SPLIT_BAM:
    input:
        get_external_bam_path
    output:
        temp(directory(SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/EXTERNAL_BAM/split_bams/"))
    singularity:
        "/gpfs/data/petljaklab/containers/samtools/samtools_1.18.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/EXTERNAL_BAM/split_bams/split.log"
    threads: 8
    resources:
        threads = 8,
        cpus = 8,
        mem_mb = 8000,
        runtime = 240,
        slurm_partition = "petljaklab,cpu_short",
        iotasks = 5,
    shell:
        "samtools split {input} -f {output}/{wildcards.sample}-%#.%. -@ {resources.threads}"

def get_split_bam(wildcards):
    outs = checkpoints.SPLIT_BAM.get(**wildcards).output[0]
    ## Get sample ID
    rid = petljakapi.translate.stringtoid(wildcards.run)
    rname = petljakapi.select.multi_select(db = db, table = "runs", filters = {"id":rid}, bench = config["bench"])[0][1]
    ## Use sample ID to get run ID & rname
    bamfile = os.path.join(outs, f"{rname}.bam")
    #print(bamfile)
    return(bamfile)

rule EXTERNAL_BAM_TO_FASTQ:
    input:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/EXTERNAL_BAM/split_bams/"
        #get_split_bam
        #SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/EXTERNAL_BAM/{analysis}/split_bams/{rname}.bam"
    output:
        reads1 = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/reads1.fq"),
        reads2 = temp(SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/reads2.fq"),
        readsU = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/readsU.fq",
        readsS = SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/readsS.fq",
    singularity:
        f"/gpfs/data/petljaklab/containers/samtools/samtools_1.18.sif"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/runs/{run}/analyses/EXTERNAL_BAM/{analysis}/fq/split.log",
    resources:
        iotasks = 2,
        mem_mb = 20000,
        runtime = 600,
        slurm_partition = "petljaklab,cpu_short"
    params:
        rname = lambda w: petljakapi.select.multi_select(db = db, table = "runs", filters = {"id":petljakapi.translate.stringtoid(w.run)}, bench = config["bench"])[0][1]
    priority: 1
    shell:
        "samtools collate -u -O {input}/{params.rname}.bam | samtools fastq -1 {output.reads1} -2 {output.reads2} -0 {output.readsU} -s {output.readsS} - &> {log};"

#rule ADD_EXTERNAL_RUNS:
#    input:
#        collate_output_fastqs
#    output:
#        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/EXTERNAL_BAM/{analysis}/splitfq/split.done"
#    resources:
#        slurm_partition = "cpu_short"
#    run:
#        for i in input:
#            i1 = i
#            i2 = i[:-4] + "_2.fq"
#            iu = i[:-4] + "_U.fq"
#            ## Get the filename
#            rname = re.sub(".*/([^.]).bam", "\\1", i)
#            ## Numeric IDs
#            studyid = petljakapi.translate.stringtoid(wildcards.study)
#            sampleid = petljakapi.translate.stringtoid(wildcards.sample)
#            ## Add run to the db, dont set path yet since we need the run_id/analysis_id for that
#            runid = petljakapi.inserts.generic_insert({"study_id":studyid, "sample_id":sampleid, "rname":rname, "source":"external_bam"}, "runs", "petljakdb_devel")[0][0]
#            run_name = petljakapi.translate.idtostring(runid, "MPR")
#            analysis_entry = {
#                "pipeline_name":"EXTERNAL_BAM",
#                "pipeline_version":EXTERNAL_BAM_VERSION,
#                "analysis_type":"FASTQ",
#                "input_table":"runs",
#                "runs_id":run_name
#            }
#            analysis_id = petljakapi.inserts.analysis_insert(analysis_entry, "analyses")[0][0]
#            analysis_name = petljakapi.translate.idtostring(analysis_id, "MPA")
#            ## Get path
#            copypath = PROD_DIR + f"studies/{wildcards.study}/samples/{wildcards.sample}/runs/{run_name}/analyses/EXTERNAL_BAM/{analysis_name}/fq/"
#            os.makedirs(copypath, exist_ok = True)
#            shutil.copy(i1, copypath + "reads1.fq")
#            shutil.copy(i2, copypath + "reads2.fq")
#            shutil.copy(iu, copypath + "readsU.fq")
#            petljakapi.update.update(db = "petljakdb_devel", table = "analyses", filters = {"id":analysis_id}, update_col = "analysis_dir", update_val = copypath)

