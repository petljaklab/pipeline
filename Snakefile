import os
import petljakdb

configfile: "config.yaml"
locals().update(configfile)

## Get basedir
basedir = os.path.abspath(__file__)

include: os.path.join(basedir, "modules/GATK_BAM/GATK_BAM.smk")
include: os.path.join(basedir, "modules/SRA/SRA.smk")

rule all:
    input:
        "__file__"