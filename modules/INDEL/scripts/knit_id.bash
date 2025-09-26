#!/bin/bash

help()
{
    echo "knit_id.bash -d DAUGHTER_TABLE -p PARENTS_TABLE -o OUTPUT_PATH"
    echo ""
    echo "Usage:"
    echo "  -d  daughter_table.txt  Table of daughters used in filter_mutations.R"
    echo "  -p  parent_table.txt    Table of parents used in filter_mutations.R"
    echo "  -o  output_path path where outputs will be written"
}

echo "1"

if [ "$#" -ne 6 ]; then
    help
    exit
fi

wd=$(pwd)


while getopts d:p:o: flag
do
    case "${flag}" in
        d) daughter=${OPTARG};;
        p) parent=${OPTARG};;
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

# Same for parent
if [[ "$parent" = /* ]]; then
    :
else
    # Append the input path to another path
    parent="$wd/$parent"
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
echo $parent
echo $output

odir=$(dirname $output)

echo ".libPaths(c(.libPaths(), \"/gpfs/share/apps/R/4.2.2/lib64/R/library\")); Rscript -e \"rmarkdown::render('/gpfs/data/petljaklab/lculibrk_prj/pipeljak/modules/INDEL/scripts/id_report.Rmd', output_file=\"$output\", knit_root_dir = \"$odir\", param=list(args=c(\"daughter\" = \"$daughter\", \"parent\" = \"$parent\")))\""
Rscript -e ".libPaths(c(.libPaths(), \"/gpfs/share/apps/R/4.2.2/lib64/R/library\")); rmarkdown::render('/gpfs/data/petljaklab/lculibrk_prj/pipeljak/modules/INDEL/scripts/id_report.Rmd', output_file=\"$output\", knit_root_dir = \"$odir\", param=list(args=c(\"daughter\" = \"$daughter\", \"parent\" = \"$parent\")))"
