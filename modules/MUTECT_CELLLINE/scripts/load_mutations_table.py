import argparse
import petljakapi
from petljakapi.connection import connection
import petljakapi.select
import petljakapi.inserts
import petljakapi.translate
import petljakapi.update
import os

parser = argparse.ArgumentParser(
    prog="load_mutations_table.py",
    description="Loads mutations from the given table into the petljakdb"
)

parser.add_argument("-a", "--analysis_id", help = "analysis ID")
parser.add_argument("-s", "--sample_id", help = "sample ID")
parser.add_argument("-u", "--study_id", help = "study ID")
parser.add_argument("-v", "--vcf", help = "path to VCF to be loaded")
parser.add_argument("-g", "--genome", help = "name of genome")
parser.add_argument("-t", "--mutation_type", help = "type of mutation (SBS, ID, SV)")
parser.add_argument("-d", "--db", help = "database name")

"""
+------------------+------------------+------+-----+---------+----------------+
| Field            | Type             | Null | Key | Default | Extra          |
+------------------+------------------+------+-----+---------+----------------+
| id               | int(10) unsigned | NO   | PRI | NULL    | auto_increment |
| genome           | varchar(255)     | YES  |     | NULL    |                |
| chromosome       | varchar(255)     | YES  |     | NULL    |                |
| start_position   | int(10) unsigned | YES  |     | NULL    |                |
| end_position     | int(10) unsigned | YES  |     | NULL    |                |
| ref              | varchar(255)     | YES  |     | NULL    |                |
| alt              | varchar(255)     | YES  |     | NULL    |                |
| ref_depth        | int(10) unsigned | YES  |     | NULL    |                |
| alt_depth        | int(10) unsigned | YES  |     | NULL    |                |
| mutation_type    | varchar(255)     | YES  |     | NULL    |                |
| mutation_filter  | varchar(255)     | YES  |     | NULL    |                |
| mutation_comment | varchar(255)     | YES  |     | NULL    |                |
| analysis_id      | int(10) unsigned | YES  | MUL | NULL    |                |
| sample_id        | int(10) unsigned | YES  | MUL | NULL    |                |
| study_id         | int(10) unsigned | YES  | MUL | NULL    |                |
+------------------+------------------+------+-----+---------+----------------+
"""

args = parser.parse_args()
## Write the mutation file to disk as a temp file, then use mySQL to load it in

## TODO: Uncomment later
genome = args.genome
db = args.db
analysis_id, sample_id, study_id = (petljakapi.translate.stringtoid(args.analysis_id),
                                    petljakapi.translate.stringtoid(args.sample_id),
                                    petljakapi.translate.stringtoid(args.study_id))

with open(args.vcf) as f:
    vcf = f.readlines()

vcf = [line.strip().split("\t") for line in vcf]

## Get and parse the header
it = 0
while vcf[it][0][:2] == "##":
    it = it+1

chrom_col = vcf[it].index("#CHROM")
pos_col = vcf[it].index("POS")
ref_col = vcf[it].index("REF")
alt_col = vcf[it].index("ALT")
filter_col = vcf[it].index("FILTER")
format_col = vcf[it].index("FORMAT")
geno_col = vcf[it].index("tumor")

#print(vcf[it:it+5])
## For now we will only import PASS variants
table_to_import = []
## Multiallelic are all filtered by Mutect2, but we will store this for future use
multiallelic_list = []
for line in vcf[(it+1):]:
    if line[filter_col] != "PASS":
        continue
    fmt = line[format_col].split(":")
    ad_ind = fmt.index("AD")
    lab_ind = fmt.index("LAB")
    geno = line[geno_col].split(":")
    try:
        ref_depth, alt_depth = geno[ad_ind].split(",")
    except ValueError:
        ## Multiallelic
        multiallelic_list.append(line)
        continue
    table_to_import.append([genome, line[chrom_col], line[pos_col], line[pos_col], line[ref_col], line[alt_col], ref_depth, alt_depth, args.mutation_type, line[filter_col], geno[lab_ind], str(analysis_id), str(sample_id), str(study_id)])    

table_to_export = [",".join(line) + "\n" for line in table_to_import]

temp_path = args.vcf + ".dbtmp"

with open(temp_path, "w") as f:
    f.writelines(table_to_export)

cursor = connection.cursor()
petljakapi.dbs.chdb(db, cursor)
cursor.execute(f"""
    LOAD DATA LOCAL INFILE '{temp_path}'
    INTO TABLE mutations
    FIELDS TERMINATED BY ','
    LINES TERMINATED BY '\n'
    (genome, chromosome, start_position, end_position, ref, alt, ref_depth, alt_depth, mutation_type, mutation_filter, mutation_comment, analysis_id, sample_id, study_id)
""")
connection.commit()
connection.close()