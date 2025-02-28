rule GRIDSS:
    input:
        fa = lambda wildcards: ALN_REFERENCES[wildcards.reference],
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)[0],
        parent_merge = lambda wildcards: parent_cell(wildcards),
    output:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/gridss.vcf",
        bam = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/assembly.bam",
    params:
        workdir = lambda w: SCRATCH_DIR + f"studies/{w.study}/samples/{w.sample}/analyses/SV/{w.analysis}/{w.reference}/gridss",
    singularity:
        f"/gpfs/data/petljaklab/containers/purple/new_container.sif"
    resources:
        threads = 2,
        cpus = 2,
        cpu = 2,
        mem_mb = 16000,
        slurm_partition = "petljaklab",
        runtime = 24*60*5,
    threads: 2
    benchmark:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/gridss.resources",
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/gridss.log",
    shell:
        """
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        ln -s {input.cell_merge} {params.workdir}/tumor.cram
        ln -s {input.cell_merge}.crai {params.workdir}/tumor.cram.crai
        ln -s {input.parent_merge} {params.workdir}/normal.cram
        ln -s {input.parent_merge}.crai {params.workdir}/normal.cram.crai
        gridss --labels parent,daughter \
               --output {output.vcf}.gz \
               --reference {input.fa} \
               --threads {threads} \
               --assembly {output.bam} \
               --picardoptions VALIDATION_STRINGENCY=LENIENT \
               --jvmheap 14g \
               --otherjvmheap 14g \
               -w {params.workdir} \
               {params.workdir}/normal.cram \
               {params.workdir}/tumor.cram
        gunzip {output.vcf}.gz
        """

rule GRIDSS_SINGLESAMP:
    input:
        fa = lambda wildcards: ALN_REFERENCES[wildcards.reference],
        cell_merge = lambda wildcards: gateway("WGS_MERGE_BAM", wildcards.sample, scratch_dir = SCRATCH_DIR, prod_dir = PROD_DIR, db = "petljakdb", ref = wildcards.reference)[0],
    output:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss_singlesamp/gridss.vcf",
        bam = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss_singlesamp/assembly.bam",
    params:
        workdir = lambda w: SCRATCH_DIR + f"studies/{w.study}/samples/{w.sample}/analyses/SV/{w.analysis}/{w.reference}/gridss_singlesamp",
    singularity:
        f"/gpfs/data/petljaklab/containers/purple/new_container.sif"
    resources:
        threads = 2,
        cpus = 2,
        cpu = 2,
        mem_mb = 16000,
        slurm_partition = "petljaklab",
        runtime = 24*60*5,
    threads: 2
    benchmark:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss_singlesamp/gridss.resources",
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss_singlesamp/gridss.log",
    shell:
        """
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        ln -s {input.cell_merge} {params.workdir}/tumor.cram
        ln -s {input.cell_merge}.crai {params.workdir}/tumor.cram.crai
        gridss --labels parent \
               --output {output.vcf}.gz \
               --reference {input.fa} \
               --threads {threads} \
               --assembly {output.bam} \
               --picardoptions VALIDATION_STRINGENCY=LENIENT \
               --jvmheap 14g \
               --otherjvmheap 14g \
               -w {params.workdir} \
               {params.workdir}/tumor.cram
        gunzip {output.vcf}.gz
        """

rule GRIDSS_SINGLESAMP_FILTER:
    input:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss_singlesamp/gridss.vcf",
    output:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss_singlesamp/gridss_filtered.vcf",
    resources:
        threads = 1,
        cpus = 1,
        cpu = 1,
        mem_mb = 4000,
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        runtime = 60,
    log: SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss_singlesamp/gridss_filtered.log",
    singularity: 
        "/gpfs/data/petljaklab/containers/gatk/gatk_4.4.0.0.sif"
    shell:
        "gatk SelectVariants --exclude-filtered -V {input.vcf} -O {output.vcf} &> {log}"

rule GRIDSS_JUNCTION_FILTER:
    input:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss_singlesamp/gridss_filtered.vcf",
    output:
        rds = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss_singlesamp/somatic.filtered.gnomAD.sv.rds"
    resources:
        threads = 1,
        cpus = 1,
        cpu = 1,
        mem_mb = 64000,
        slurm_partition = config["clusters"][config["parts"]]["fat"]
    singularity:
        "/gpfs/data/petljaklab/containers/ggnome/ggnome_latest.sif"
    params:
        outdir = lambda w: SCRATCH_DIR + f"studies/{w.study}/samples/{w.sample}/analyses/SV/{w.analysis}/{w.reference}/gridss_singlesamp/"
    shell:
        """
        Rscript modules/SV/scripts/junction_filter.R \
        --gnomAD /gpfs/data/petljaklab/resources/hg19/pipeline_resources/sv/sv/GRCh19.gnomad_v2.1_sv.merged.rds \
        --sv {input.vcf} \
        --padding 1000 \
        --pon /gpfs/data/petljaklab/resources/hg19/pipeline_resources/sv/sv/ggnome_gridss_pon_hg19_new.rds \
        -o {params.outdir}
        """

rule GRIDSS_SOMATIC_FILTER:
    input:
        vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/gridss.vcf",
    output:
        vcf_h = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/high_confidence_somatic.vcf.bgz",
        vcf_lh = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/high_low_confidence_somatic.vcf.bgz",
    resources:
        threads = 1,
        cpus = 1,
        cpu = 1,
        mem_mb = 10000,
        slurm_partition = "petljaklab",
        runtime = 60*4,
    singularity: "/gpfs/data/petljaklab/containers/purple/new_container.sif"
    log: vcf_h = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/gridss_somatic.log",
    params:
        vcf_output = lambda w: f"{SCRATCH_DIR}studies/{w.study}/samples/{w.sample}/analyses/SV/{w.analysis}/{w.reference}/gridss/high_confidence_somatic.vcf",
        vcf_hl_output = lambda w: f"{SCRATCH_DIR}studies/{w.study}/samples/{w.sample}/analyses/SV/{w.analysis}/{w.reference}/gridss/high_low_confidence_somatic.vcf",
        outdir = lambda w: f"{SCRATCH_DIR}studies/{w.study}/samples/{w.sample}/analyses/SV/{w.analysis}/{w.reference}/gridss/"
    shell:
        """
        gridss_somatic_filter -p /gpfs/data/petljaklab/resources/hg19/pipeline_resources/sv/sv/ \
                              -i {input.vcf} \
                              --plotdir {params.outdir} \
                              -o {params.vcf_output} \
                              --fulloutput {params.vcf_hl_output} \
                              --scriptdir /opt/gridss/ 2> {log}
        
        """

rule ANNOTATE_GRIDSS_OUTPUT:
    input:
        fa = lambda wildcards: ALN_REFERENCES[wildcards.reference],
        vcf_lh = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/high_low_confidence_somatic.vcf.bgz",
    output:
        f_vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/high_low_confidence_annotated.vcf",
    resources:
        threads = 1,
        cpus = 1,
        cpu = 1,
        mem_mb = 3000,
        slurm_partition = "petljaklab",
        runtime = 10,
    shell:
        """
        module load r/4.1.2
        Rscript /gpfs/data/petljaklab/lculibrk_prj/pipeljak/modules/SV/scripts/simple_sv_annotation.R -v {input.vcf_lh} -r {input.fa} -o {output.f_vcf}
        """

rule ANNOTATED_TO_TBL:
    input:
        f_vcf = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/high_low_confidence_annotated.vcf",
    output:
        f_tbl = SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/high_low_confidence_annotated.txt",
    resources:
        mem_mb = 3000,
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        runtime = 60 * 4,
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/SV/{analysis}/{reference}/gridss/high_low_confidence_totable.log"
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    shell:
        """
        gatk VariantsToTable -V {input.f_vcf} \
            -O {output.f_tbl} -raw &> {log}
        """

rule SV_META_TABLES:
    output:
        daughter = SCRATCH_DIR + "studies/{study}/analyses/SV/{reference}/daughters_table.txt",
        parents = SCRATCH_DIR + "studies/{study}/analyses/SV/{reference}/parents_table.txt",
    threads: 1
    resources:
        slurm_partition = config["clusters"][config["parts"]]["dev"],
        cpus = 1,
        threads = 1,
        mem_mb = 2000,
        runtime = 15,
    run:
        if not os.path.exists(f"{SCRATCH_DIR}studies/{wildcards.study}/analyses/SV/{wildcards.reference}"):
            os.makedirs(f"{SCRATCH_DIR}studies/{wildcards.study}/analyses/SV/{wildcards.reference}")
        study_id = petljakapi.translate.stringtoid(wildcards.study)
        samples_ids = petljakapi.select.multi_select(db, "samples", {"study_id":study_id})
        for sample in samples_ids:
            cell_name = petljakapi.select.multi_select(db, "cells", {"id":sample[6]})[0][1]
            ## Get analysis row
            analysis = petljakapi.select.multi_select(db = db, table = "analyses", filters = {"samples_id":sample[0], "pipeline_name":"SV", "reference_genome":wildcards.reference}, bench = config["bench"])
            if len(analysis) == 0:
                continue
            analysis = analysis[0]
            aid = petljakapi.translate.idtostring(analysis[0], "MPA")
            ap = analysis[10]
            ## Test if this sample is used as a reference for any other samples
            daughters = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "sample_parent_id", filter_value = sample[0], bench = config["bench"])
            parent_id = sample[5]
            if parent_id is not None and parent_id != "NULL": ## if this sample has a parent ie. it's a daughter
                files = f"{wildcards.reference}/gridss/high_low_confidence_annotated.vcf"
                files = os.path.join(ap, aid, files)
                parent_name = petljakapi.select.simple_select(db = db, table = "samples", filter_column = "id", filter_value = parent_id)[0][1]
                with open(output.daughter, "a") as f:
                    f.write("%s\n" % '\t'.join([files, sample[1], parent_name, cell_name, sample[4], str(sample[7])]))
            if len(daughters) > 0 or parent_id == "NULL" or parent_id is None: ## ie. if this is a parent
                files = f"{wildcards.reference}/gridss_singlesamp/"
                files = os.path.join(ap, aid, files)
                with open(output.parents, "a") as f:
                    f.write("%s\n" % '\t'.join([files, sample[1], cell_name, sample[4]]))