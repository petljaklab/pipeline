#!/bin/bash

module load r/4.1.2 pandoc/2.2.3.2 latex/2019

help()
{
    echo "knit_sbs.bash -d DAUGHTER_TABLE -o OUTPUT_PATH"
    echo ""
    echo "Usage:"
    echo "  -d  daughter_table.txt  Table of daughters used in filter_mutations.R"
    echo "  -o  output_path path where outputs will be written"
}


if [ "$#" -ne 4 ]; then
    help
    exit
fi

wd=$(pwd)


while getopts d:o: flag
do
    case "${flag}" in
        d) daughter=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

# Check if the input path is absolute or not
if [[ "$daughter" = /* ]]; then
    :
else
    # Append the input path to another path
    daughter="$wd/$daughter"
fi

if [[ "$output" = /* ]]; then
    :
else
    # Append the input path to another path
    output="$wd/$output"
fi

if [[ -f "$daughter" ]]; then
    :
else
    echo "Daughter table not found. Tried: $daughter"
fi

echo $daughter
echo $output

odir=$(dirname $output)

echo "Rscript -e \"rmarkdown::render('/gpfs/data/petljaklab/lculibrk_prj/pipeljak/modules/MUTECT_CELLLINE/scripts/sbs_report.Rmd', output_file=\"$output\", knit_root_dir = \"$odir\", param=list(args=c(\"daughter\" = \"$daughter\")))\""
Rscript -e "rmarkdown::render('/gpfs/data/petljaklab/lculibrk_prj/pipeljak/modules/MUTECT_CELLLINE/scripts/sbs_report.Rmd', output_file=\"$output\", knit_root_dir = \"$odir\", param=list(args=c(\"daughter\" = \"$daughter\")))"
