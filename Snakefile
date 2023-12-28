import os
import petljakapi
from petljakapi import connection
import petljakapi.select
import petljakapi.inserts


configfile: "config.yaml"
locals().update(config)

## Get basedir
basedir = os.path.abspath(workflow.basedir)

include: os.path.join(basedir, "modules/GATK_BAM/GATK_BAM.smk")
include: os.path.join(basedir, "modules/SRA/SRA.smk")
include: os.path.join(basedir, "modules/TESTS/TESTS.smk")

rule all:
    input:
        __file__