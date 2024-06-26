# petljaklab bioinformatics pipelines

Bioinformatics workflows for the [Petljak Lab](https://petljaklab.com/) at the Perlmutter Cancer Center, written by [Luka Culibrk](https://github.com/lculibrk)

## Quick Execute
```python ./executor.py [--idfile ID_FILE | --id SAMPLE_ID] --pipeline PIPELINE_NAME [snakemake options]```

## Overview

The petljaklab pipelines are a collection of [Snakemake](https://snakemake.readthedocs.io/en/stable/) modules, each of which does a specific task in genomics, e.g. WGS alignment, or variant calling in a cell line parent-daughter context. 

Granted, many different workflows make, for example, a VCF. Which workflow is run, is determined from metadata that is stored on the [petljakDB](https://github.com/petljaklab/petljakdb). 

The executor function, `executor.py` is given a sample (or list of samples) alongside an endpoint. For example, `WGS_MERGE_BAM` is the endpoint for making a WGS cram. 
The function determines based on the sample metadata on the petljakdb which module should be used to generate the cram file. 
Similarly, for a parent-daughter cell line experiment, the pipeline determines the sample's matched parental sample as well as its sibling clones based on the database metadata. 
Refer to the database/API repo for information on how the data are stored there.

Current endpoints:
- FASTQ generation
- Alignment
- SNV calling with Mutect2
- Indel calling with Mutect2, Strelka2, and Varscan2
- Somatic (convenience endpoint for SNV + Indel + other endpoints to be added)

Each endpoint has at least one module that can generate that endpoint. 

FASTQ:
- SRA-hosted data (`SRA`)
- EGA-hosted data (`EGA`)
- Local BAM files to be remapped (`EXTERNAL_BAM`)

Alignment:
- GATK WGS best practices CRAM (`GATK_BAM`)


SNV calling:
- Parent-daughter cell line (`MUTECT_CELLLINE`)

Indel calling:
- Tumor-normal (`INDEL`)

The `SOMATIC` endpoint forces both SNV and Indel endpoints to be executed.


