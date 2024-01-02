## Probably this should all go in the database eventually but for now this is easier for prototyping

db_deps = {
    "FASTQ": [
        "studies",
        "samples",
        "runs"
    ],
    "WGS_MERGE_BAM": [
        "studies",
        "samples"
    ]
}

## Describes the endpoint of a pipeline
## Controlled vocabulary: [FASTQ, MERGE_BAM, ]
module_outputs = {
    "TESTS":"FASTQ",
    "SRA":"FASTQ",
    "GATK_BAM":"WGS_MERGE_BAM"
}

filetypes = {
    "FASTQ": [
        "TESTS",
        "SRA"
    ],
    "WGS_MERGE_BAM": [
        "GATK_BAM"
    ]
}

## Empty list if the module requires no input from a module
## Controlled vocabulary: [FASTQ, WGS_MERGE_BAM]
module_inputs = {
    "TESTS":[],
    "SRA":[],
    "GATK_BAM":["FASTQ"]
}

## e.g. we can see that GATK_BAM requires a FASTQ module as input and outputs a WGS_MERGE_BAM