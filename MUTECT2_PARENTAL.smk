rule COPY_RENAME_DAUGHTER_CALLS:
    input:
        SCRATCH_DIR + "studies/{study}/samples/{daughter_sample}/analyses/MUTECT_CELLLINE/{daughter_analysis}/mutect2/{reference}/germ/raw_{chrom}.vcf"
    output:
        SCRATCH_DIR + "studies/{study}/samples/{daughter_sample}/analyses/MUTECT_CELLLINE/{daughter_analysis}/mutect2/{reference}/germ/renamed_{chrom}.vcf"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{daughter_sample}/analyses/MUTECT_CELLLINE/{daughter_analysis}/mutect2/{reference}/germ/renamed_{chrom}.log"
    resources:
        threads = 1,
        mem_mb = 5000,
        runtime = 15,
        slurm_partition = "cpu_dev"
    shell:
    ## First copy over the header
    ## Next, grep out the header, cut out the parental column, then use awk to rename the daughter column to {parental-sample-name}_d, then remove variants with ref/alt length longer than the reads
    ## Otherwise, it does not play nice with mutect2
        """
        grep '##' {input} > {output} 2> {log};
        TUMOR_NAME=$(grep '##tumor_sample' {input} | sed 's/##tumor_sample=//g');
        column_number=$(awk -v search="$TUMOR_NAME" -F'\\t' '{{ for (i=1; i<=NF; i++) {{ if ($i == search) {{ print i; exit }} }} }}' {input})
        grep -v '##' {input} | cut -f 1,2,3,4,5,6,7,8,9,$column_number | awk -v OFS="\\t" 'NR==1{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, "daughter"}}NR>1{{print $0}}' | awk -v OFS="\\t" 'length($4) <= 50 && length($5) <= 50{{print $0}}' >> {output} 2>> {log}
        """    

def all_daughter_calls_function(wildcards):
    ## Function to take in wildcards (ie. the sample) and output a list of all the daughter calls as renamed above
    ##
    daughters = petljakapi.select.multi_select(db = db, table = "samples", filters = {"sample_parent_id":petljakapi.translate.stringtoid(wildcards.sample)})
    paths = []
    for daughter in daughters:
        ## id
        daught_id = daughter[0]
        ## Run gateway to ensure the analysis dir is created
        gateway("MUTECT", petljakapi.translate.idtostring(daught_id, "MPS"), SCRATCH_DIR, config["PROD_DIR"], db = "petljakdb")
        ## Now we need the relevant analysis ID for the daughter sample
        analysis = petljakapi.select.multi_select(db = db, table = "analyses", filters = {"pipeline_name":"MUTECT_CELLLINE", "samples_id":daught_id})
        p = analysis[0][10]
        chromstring = "renamed_" + str(wildcards.chrom) + ".vcf"
        vcf_path = os.path.join(p, petljakapi.translate.idtostring(analysis[0][0], "MPA"), "mutect2", wildcards.reference, "germ", chromstring)
        paths.append(vcf_path)
    return(paths)


rule COMBINE_MUTECT2_CELLLINE_DAUGHTER_CALLS:
    input:
        all_daughter_calls_function
    output:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/daughters_merged_vcf/{chrom}.vcf"
    singularity:
        f"/gpfs/data/petljaklab/containers/gatk/gatk_{GATK_VERSION}.sif"
    resources:
        threads = 1,
        mem_mb = 5000,
        runtime = 15,
        slurm_partition = "cpu_dev"
    log:
        SCRATCH_DIR + "studies/{study}/samples/{sample}/analyses/MUTECT_CELLLINE/{analysis}/mutect2/{reference}/daughters_merged_vcf/{chrom}.log"
    params:
        inputlist = lambda wildcards, input: " -I ".join([input]) if isinstance(input, str) else " -I ".join(input)
    shell:
        "gatk MergeVcfs -I {params.inputlist} -O {output} &> {log}"