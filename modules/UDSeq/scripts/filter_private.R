#!/usr/bin/env Rscript
library(getopt)
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help'     , 'h', 0, "logical",
	'variants' , 'v', 1, "character",
	'genotypes', 'g', 1, "character",
	'output'   , 'o', 1, "character",
	'threads'  , 't', 2, "integer"
), byrow=TRUE, ncol=4)


args = getopt(spec)
if(!all(c("variants", "genotypes", "output") %in% names(args))){
	print(getopt(spec, usage = T))
	stop()
}

if(!'threads' %in% names(args)){
	threads = 1
}else{
	threads = args$threads
}

library(data.table)
message(paste0("Using ", threads, " threads for data.table"))
setDTthreads(threads)
library(stringr)

vars = fread(args$variants)
geno = fread(args$genotypes)

geno[,n:=.N, by = c("V1", "V2")]
bad = unique(geno[n > 1])
bad = bad[,1:2]
setnames(bad, c("CHROM", "POS"))
setkey(bad, CHROM, POS)
setkey(vars, CHROM, POS)
bad_muts = vars[bad][!is.na(TUMOR.DP)]
bad_muts$FILTER = "shared"

good = unique(geno[n == 1])
good = good[,1:2]
setnames(good, c("CHROM", "POS"))
setkey(good, CHROM, POS)

good_muts = vars[good][!is.na(TUMOR.DP)]

muts = rbind(good_muts, bad_muts)
muts = muts[order(CHROM, POS)]
fwrite(muts, args$output, sep = "\t")