library("ExomeDepth")
library("GenomicRanges")
library("optparse")
 
option_list <- list(
    make_option(c('--bam'), type = "character", default = NULL),
    make_option(c('--sample-name'), type = "character", default = NULL),  
    make_option(c('--min-mapq'), type = "numeric", default = NULL),
    make_option(c('--cohort-rdata'), type = "character", default = NULL),
    make_option(c('--out-base'), type = "character", default = NULL)

)

opt <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)

# Load: normal_bams, bed_df, ref_fasta, cohort_count_df, cohort_count_matrix
load(opt$cohort_rdata) 

# Counts for all bams
bam_counts <- getBamCounts(
    bed.frame = bed_df, bam.files = opt$bam, min.mapq = opt$min_mapq, referenceFasta = ref_fasta
)

# create count matrix
case_count_df <- as(bam_counts[, colnames(bam_counts)], "data.frame")
#case_count_df$chromosome <- gsub(as.character(case_count_df$space), pattern = "chr", replacement = "")
case_count_df$chromosome <- as.character(case_count_df$space)

# sample count is 7th column
case_count_matrix <- as.matrix(case_count_df[, 7])
# choose best reference for sample
closest_reference <- select.reference.set(
    test.counts = case_count_matrix,
    reference.counts = cohort_count_matrix,
    bin.length = (cohort_count_df$end - cohort_count_df$start)/1000,
    n.bins.reduced = 10000
)

reference_matrix <- as.matrix(cohort_count_df[, closest_reference$reference.choice, drop = FALSE])
selected_reference <- apply(X = reference_matrix, MAR = 1, FUN = sum)

# call CNVs
cnv_calls <- new(
    'ExomeDepth',
    test = case_count_matrix[, 1],
    reference = selected_reference,
    formula = 'cbind(test, reference) ~ 1'
)

cnv_calls <- CallCNVs(
    x = cnv_calls,
    transition.probability = 0.01,
    chromosome = cohort_count_df$chromosome,
    start = cohort_count_df$start,
    end = cohort_count_df$end,
    name = cohort_count_df$names
)


# save data
rdata_path <- paste(opt$out_base, ".Rdata", sep = "")
print(rdata_path)
save(closest_reference, cnv_calls, file = rdata_path)

cnv_table <- cnv_calls@CNV.calls
cnv_table$chrom <- cnv_table$chromosome

cnv_call_path <- paste(opt$out_base, ".txt", sep = "")
write.table(cnv_calls@CNV.calls, file = cnv_call_path, sep= "\t", row.names = FALSE, quote = FALSE)