### This script is there to fix analyses that were not properly marked complete in the past.
import petljakapi.select
import petljakapi.translate
import petljakapi.update
import sys
import os
import datetime
from lib.input_functions import gateway

#analysis_name, given_id, scratch_dir, prod_dir, db = "petljakdb", bench = False, ref = "hg19"
scratch_dir = "/gpfs/data/petljaklab/data/"
prod_dir = scratch_dir

analyses = petljakapi.select.multi_select("petljakdb", "analyses", filters = {}, headers= True)
anames = analyses[0]
analyses = analyses[1:]
print(anames)

def mark_analysis_complete(analysis_id, path, newval, dryrun = False):
    petljakapi.update.update("petljakdb", "analyses", {"id":analysis_id}, "analysis_complete", newval, dryrun = dryrun)
    petljakapi.update.update("petljakdb", "analyses", {"id":analysis_id}, "analysis_dir", os.path.dirname(path), dryrun = dryrun)


dryrun = False

for ana in analyses:
    aid = petljakapi.translate.idtostring(ana[0], "MPA")
    #print(ana)
    if ana[anames.index("runs_id")] is not None:
        term_id = petljakapi.translate.idtostring(ana[anames.index("runs_id")], "MPR")
    elif ana[anames.index("samples_id")] is not None:
        term_id = petljakapi.translate.idtostring(ana[anames.index("samples_id")], "MPS")
    ana[anames.index("reference_genome")]
    ## Can't use gateway - need to manually pattern-assign based on existing pipelines
    ## This is what happens when you're lazy and don't implement this in the first place
    #critical_files = gateway(ana[anames.index("analysis_type")], term_id, scratch_dir, prod_dir, "petljakdb", bench = False, ref = ana[anames.index("reference_genome")], dryrun = True)
    #print(critical_files)
    pipe = ana[anames.index("pipeline_name")]
    atype = ana[anames.index("analysis_type")]
    path = ana[anames.index("analysis_dir")]
    hgver = ana[anames.index("reference_genome")]
    validate = False
    if pipe == "GATK_BAM":
        ## Check if the terminal path already exists, and if so, it's easy
        if os.path.exists(os.path.join(path, "merged.cram")):
            validate = True
            mark_analysis_complete(ana[0], os.path.join(path, "merged.cram"), "True", dryrun = dryrun)
        elif os.path.exists(os.path.join(path, aid, "merge", hgver, "merged.cram")):
            validate = True
            mark_analysis_complete(ana[0], os.path.join(path, aid, "merge", hgver, "merged.cram"), "True", dryrun = dryrun)
        else:
            validate = False
            mark_analysis_complete(ana[0], os.path.join(path, "merged.cram"), "False", dryrun = dryrun)
    elif atype == "FASTQ":
        validate = True
        mark_analysis_complete(ana[0], os.path.join(path, "filler"), "False", dryrun = dryrun)
    elif pipe == "MUTECT_CELLLINE":
        #MPA007162/mutect2/hg19/proc/
        if os.path.exists(os.path.join(path, aid, "mutect2", hgver, "proc/variants_final.vcf")):
            mark_analysis_complete(ana[0], os.path.join(path, aid, "mutect2", hgver, "proc/variants_final.vcf"), "True", dryrun = dryrun)
            validate = True
        else:
            samp_id = petljakapi.select.multi_select("petljakdb", "samples", filters = {"id":ana[anames.index("samples_id")]})[0]
            if samp_id[5] is None:
                if os.path.exists(os.path.join(path, aid, "mutect2", hgver, "parental/filtered.vcf")):
                    mark_analysis_complete(ana[0], os.path.join(path, aid, "mutect2", hgver, "parental/filtered.vcf"), "True", dryrun = dryrun)
                    validate = True
    elif pipe == "EXTERNAL_BAM":
        validate = True
    elif pipe == "LOAD_EXT_BAM":
        validate = True
    elif pipe == "INDEL":
        #MPA007183/hg19/mutect2/indels.vcf
        if os.path.exists(os.path.join(path, aid, hgver, "mutect2/indels.vcf")):
            validate = True
            mark_analysis_complete(ana[0], os.path.join(path, aid, hgver, "mutect2/indels.vcf"), "True", dryrun = dryrun)
    elif pipe == "MUTECT_BIOP":
        if os.path.exists(os.path.join(path, aid, "mutect2", hgver, "biop/filtered.vcf")):
            validate = True
            mark_analysis_complete(ana[0], os.path.join(path, aid, "mutect2", hgver, "biop/filtered.vcf"), "True", dryrun = dryrun)
    elif pipe == "UDSEQ_BAM":
        #MPA008416/merge/hg38/calls/calls_snv_private.txt
        if os.path.exists(os.path.join(path, aid, "merge", hgver, "calls/calls_snv_private.txt")):
            mark_analysis_complete(ana[0], os.path.join(path, aid, "merge", hgver, "calls/calls_snv_private.txt"), "True", dryrun = dryrun)
            validate = True
    elif pipe == "SV":
        validate = False
    else:
        print(ana)
    #validate = [os.path.exists(f) for f in critical_files]
#    if not validate:
#        if atype not in ["WGS_MERGE_BAM", "FASTQ", "MUTECT_CELLLINE"]:
#            print(ana)
    #        sys.exit()