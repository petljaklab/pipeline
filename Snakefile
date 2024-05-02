import os
import petljakapi
from petljakapi import connection
import petljakapi.select
import petljakapi.inserts
from lib.input_functions import gateway

configfile: "config.yaml"
locals().update(config)

## Get basedir
basedir = os.path.abspath(workflow.basedir)

include: os.path.join(basedir, "modules/GATK_BAM/GATK_BAM.smk")
include: os.path.join(basedir, "modules/SRA/SRA.smk")
include: os.path.join(basedir, "modules/EGA/EGA.smk")
include: os.path.join(basedir, "modules/TESTS/TESTS.smk")
include: os.path.join(basedir, "modules/MUTECT_CELLLINE/MUTECT_CELLLINE.smk")
include: os.path.join(basedir, "modules/EXTERNAL_BAM/EXTERNAL_BAM.smk")
include: os.path.join(basedir, "modules/EXTERNAL_BAM/PREP_EXTERNAL_BAM.smk")

#{prod_dir} OR {scratch_dir}/
#   └──studies/
#       └──flatfile.txt
#       └──{studies.id}/
#           └──analyses/
#               └──{pipeline_name}
#                   └──{analysis_id}
#           └──samples/
#               └──flatfile.txt
#               └──{samples.id}/
#                   └──analyses/
#                       └──{pipeline_name}
#                           └──{analysis.id}
#print(snakemake.get_argument_parser().parse_args().target)
#print(DAG.requested_files)

#def pipeline_init(terminal_id, id_type = None, analysis_name):
    ## Takes a given ID and analysis name and does a few things
    ## First thing we need to do is get an analysis ID, either from existing matching analysis or make a new db entry with pending status
    ## Return the relevant path with AIDs filled in
    ## e.g. file path: "{arbitrary path}/studies/MPP000001/samples/MPS000001/runs/MPS000001/analyses/TESTS/MPA000001/flatfile.txt"
    

wildcard_constraints:
    sample="MPS[0-9]+"

rule all:
    input:
        __file__