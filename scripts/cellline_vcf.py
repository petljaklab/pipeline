import argparse
import re
import os

parser = argparse.ArgumentParser()

parser.add_argument("-v", "--vcf", help = "Input VCF")
parser.add_argument("-c", "--cmdline", help = "commandline used to generate the vcf we're reheadering")
#parser.add_argument("-o", "--output", required=False)

args = parser.parse_args()

## We're editing the header
pattern = r'##'

# Open the text file for reading
## This is a hack but it's few enough files to be ok
lines = []
with open(args.vcf, 'r') as file:
    for line in file:
        # Strip newline characters from the line
        line = line.strip()
        # Check if the line matches the regular expression pattern
        if re.match(pattern, line):
            lines.append(line)  # You can process or store the matched line here
        else:
            break  # Stop reading further lines if the pattern doesn't match

## Now first handle the new filters
filt_ind = max([index for index, value in enumerate(lines) if re.match("##FILTER", value)])+1
lines.insert(filt_ind, r'##FILTER=<ID=shared_external,Description="variant found in unrelated daughters">')
lines.insert(filt_ind, r'##FILTER=<ID=shared_parental,Description="variant found in >50% of study parents">')
lines.insert(filt_ind, r'##FILTER=<ID=depth,Description="variant has <15 read depth in parent">')


format_ind = max([index for index, value in enumerate(lines) if re.match("##FORMAT", value)]) + 1
lines.insert(format_ind, r'##FORMAT=<ID=LAB,Number=1,Type=String,Description="Label describing sharedness of mutation among cell line parents and other daughters">')

## Add the command line to this
with open(args.cmdline, "r") as file:
    cmdline = file.readlines()
mtime = os.path.getmtime(args.cmdline)

if len(cmdline) != 1:
    raise ValueError("cmdline file does not have exactly one line!")

cmdline = cmdline[0].strip()

info_ind = min([index for index, value in enumerate(lines) if re.match("##INFO", value)]) + 1

cmdline_entry = f"FilteringCommand=<ID=cellline_filtering,CommandLine={cmdline},Version=\"1.0\",Date=\"{mtime}\""
cmdline = r'##' + cmdline_entry
lines.insert(info_ind, cmdline)

[print(line) for line in lines]

##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
#Version="4.4.0.0",Date="May 20, 2024 at 10:31:06 AM GMT">
###GATKCommandLine=<ID=Mutect2,CommandLine="Mutect2 --f1r2-tar-gz /gpfs/data/petljaklab/data/studies/MPP000001/samples/MPS000229/analyses/MUTECT_CELLLINE/MPA003497/mutect2/hg19/std/f1r2_X.tar.gz --normal-sample MPS000228 --panel-of-normals /gpfs/data/petljaklab/resources/hg19/pipeline_resources/somatic_celline/reference_vcf/pon.vcf --germline-resource /gpfs/data/petljaklab/resources/hg19/pipeline_resources/somatic_celline/reference_vcf/gnomad.vcf --native-pair-hmm-threads 1 --output /gpfs/data/petljaklab/data/studies/MPP000001/samples/MPS000229/analyses/MUTECT_CELLLINE/MPA003497/mutect2/hg19/std/raw_X.vcf --intervals X --input /gpfs/data/petljaklab/data/studies/MPP000001/samples/MPS000229/analyses/GATK_BAM/MPA001217/merge/hg19/merged.cram --input /gpfs/data/petljaklab/data/studies/MPP000001/samples/MPS000228/analyses/GATK_BAM/MPA001216/merge/hg19/merged.cram --reference /gpfs/data/petljaklab/resources/hg19/pipeline_resources/genome/fasta/hg19.fa --f1r2-median-mq 50 --f1r2-min-bq 20 --f1r2-max-depth 200 --flow-likelihood-parallel-threads 0 --flow-likelihood-optimized-comp false --flow-use-t0-tag false --flow-probability-threshold 0.003 --flow-remove-non-single-base-pair-indels false --flow-remove-one-zero-probs false --flow-quantization-bins 121 --flow-fill-empty-bins-value 0.001 --flow-symmetric-indel-probs false --flow-report-insertion-or-deletion false --flow-disallow-probs-larger-than-call false --flow-lump-probs false --flow-retain-max-n-probs-base-format false --flow-probability-scaling-factor 10 --flow-order-cycle-length 4 --flow-number-of-uncertain-flows-to-clip 0 --flow-nucleotide-of-first-uncertain-flow T --keep-boundary-flows false --genotype-pon-sites false --genotype-germline-sites false --af-of-alleles-not-in-resource -1.0 --mitochondria-mode false --mutect3-training-mode false --mutect3-ref-downsample 10 --mutect3-alt-downsample 20 --mutect3-non-artifact-ratio 20 --tumor-lod-to-emit 3.0 --initial-tumor-lod 2.0 --pcr-snv-qual 40 --pcr-indel-qual 40 --max-population-af 0.01 --downsampling-stride 1 --callable-depth 10 --max-suspicious-reads-per-alignment-start 0 --normal-lod 2.2 --ignore-itr-artifacts false --gvcf-lod-band -2.5 --gvcf-lod-band -2.0 --gvcf-lod-band -1.5 --gvcf-lod-band -1.0 --gvcf-lod-band -0.5 --gvcf-lod-band 0.0 --gvcf-lod-band 0.5 --gvcf-lod-band 1.0 --minimum-allele-fraction 0.0 --independent-mates false --flow-mode NONE --disable-adaptive-pruning false --kmer-size 10 --kmer-size 25 --dont-increase-kmer-sizes-for-cycles false --allow-non-unique-kmers-in-ref false --num-pruning-samples 1 --min-dangling-branch-length 4 --recover-all-dangling-branches false --max-num-haplotypes-in-population 128 --min-pruning 2 --adaptive-pruning-initial-error-rate 0.001 --pruning-lod-threshold 2.302585092994046 --pruning-seeding-lod-threshold 9.210340371976184 --max-unpruned-variants 100 --linked-de-bruijn-graph false --disable-artificial-haplotype-recovery false --enable-legacy-graph-cycle-detection false --debug-assembly false --debug-graph-transformations false --capture-assembly-failure-bam false --num-matching-bases-in-dangling-end-to-recover -1 --error-correction-log-odds -Infinity --error-correct-reads false --kmer-length-for-read-error-correction 25 --min-observations-for-kmer-to-be-solid 20 --likelihood-calculation-engine PairHMM --base-quality-score-threshold 18 --dragstr-het-hom-ratio 2 --dont-use-dragstr-pair-hmm-scores false --pair-hmm-gap-continuation-penalty 10 --expected-mismatch-rate-for-read-disqualification 0.02 --pair-hmm-implementation FASTEST_AVAILABLE --pcr-indel-model CONSERVATIVE --phred-scaled-global-read-mismapping-rate 45 --disable-symmetric-hmm-normalizing false --disable-cap-base-qualities-to-map-quality false --enable-dynamic-read-disqualification-for-genotyping false --dynamic-read-disqualification-threshold 1.0 --native-pair-hmm-use-double-precision false --flow-hmm-engine-min-indel-adjust 6 --flow-hmm-engine-flat-insertion-penatly 45 --flow-hmm-engine-flat-deletion-penatly 45 --pileup-detection false --pileup-detection-enable-indel-pileup-calling false --num-artificial-haplotypes-to-add-per-allele 5 --artifical-haplotype-filtering-kmer-size 10 --pileup-detection-snp-alt-threshold 0.1 --pileup-detection-indel-alt-threshold 0.5 --pileup-detection-absolute-alt-depth 0.0 --pileup-detection-snp-adjacent-to-assembled-indel-range 5 --pileup-detection-bad-read-tolerance 0.0 --pileup-detection-proper-pair-read-badness true --pileup-detection-edit-distance-read-badness-threshold 0.08 --pileup-detection-chimeric-read-badness true --pileup-detection-template-mean-badness-threshold 0.0 --pileup-detection-template-std-badness-threshold 0.0 --bam-writer-type CALLED_HAPLOTYPES --dont-use-soft-clipped-bases false --override-fragment-softclip-check false --min-base-quality-score 10 --smith-waterman JAVA --emit-ref-confidence NONE --max-mnp-distance 1 --force-call-filtered-alleles false --reference-model-deletion-quality 30 --soft-clip-low-quality-ends false --allele-informative-reads-overlap-margin 2 --smith-waterman-dangling-end-match-value 25 --smith-waterman-dangling-end-mismatch-penalty -50 --smith-waterman-dangling-end-gap-open-penalty -110 --smith-waterman-dangling-end-gap-extend-penalty -6 --smith-waterman-haplotype-to-reference-match-value 200 --smith-waterman-haplotype-to-reference-mismatch-penalty -150 --smith-waterman-haplotype-to-reference-gap-open-penalty -260 --smith-waterman-haplotype-to-reference-gap-extend-penalty -11 --smith-waterman-read-to-haplotype-match-value 10 --smith-waterman-read-to-haplotype-mismatch-penalty -15 --smith-waterman-read-to-haplotype-gap-open-penalty -30 --smith-waterman-read-to-haplotype-gap-extend-penalty -5 --flow-assembly-collapse-hmer-size 0 --flow-assembly-collapse-partial-mode false --flow-filter-alleles false --flow-filter-alleles-qual-threshold 30.0 --flow-filter-alleles-sor-threshold 3.0 --flow-filter-lone-alleles false --flow-filter-alleles-debug-graphs false --min-assembly-region-size 50 --max-assembly-region-size 300 --active-probability-threshold 0.002 --max-prob-propagation-distance 50 --force-active false --assembly-region-padding 100 --padding-around-indels 75 --padding-around-snps 20 --padding-around-strs 75 --max-extension-into-assembly-region-padding-legacy 25 --max-reads-per-alignment-start 50 --enable-legacy-assembly-region-trimming false --interval-set-rule UNION --interval-padding 0 --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --max-variants-per-shard 0 --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --sites-only-vcf-output false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --gcs-project-for-requester-pays  --disable-tool-default-read-filters false --minimum-mapping-quality 20 --max-read-length 2147483647 --min-read-length 30 --disable-tool-default-annotations false --enable-all-annotations false",Version="4.4.0.0",Date="May 20, 2024 at 10:31:06 AM GMT">