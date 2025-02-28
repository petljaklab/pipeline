library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(argparse)
library(data.table)

parser <- ArgumentParser()
parser$add_argument("-v", "--vcf", help="gridss vcf")
parser$add_argument("-o", "--outvcf", help="Output VCF")
parser$add_argument("-r", "--ref", help = "reference genome")
args <- parser$parse_args()


f = VariantAnnotation::readVcf(args$vcf, args$ref)

new_info_field <- DataFrame(
	Number = "1",
	Type = "String",
	Description = "Simple SV variant annotation",
	row.names = "SIMPLE_TYPE"
)

hdr <- header(f)
info(hdr) <- rbind(info(hdr), new_info_field)
header(f) = hdr

gr <- breakpointRanges(f)
svtype <- simpleEventType(gr)
info(f)$SIMPLE_TYPE = NA_character_ 
info(f[gr$sourceId])$SIMPLE_TYPE <- svtype





#f <- VCF(rowRanges = gr, colData = colData(f), fixed = fixed(f), info = info(f), header = hdr)

VariantAnnotation::writeVcf(f, args$outvcf)
