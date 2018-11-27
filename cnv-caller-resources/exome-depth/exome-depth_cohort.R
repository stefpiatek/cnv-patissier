library("ExomeDepth")
library("GenomicRanges")
library("optparse")
 
option_list <- list(
    make_option(c('--bam-table'), type = "character", default = NULL),
    make_option(c('--capture-bed'), type = "character", default = NULL),
    make_option(c('--ref-fasta'), type = "character", default = NULL),
    make_option(c('--min-mapq'), type = "numeric", default = NULL),
    make_option(c('--out-path'), type = "character", default = NULL)

)

opt <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)

bam_df <- read.table(opt$bam_table, header = TRUE, stringsAsFactors = FALSE)
normal_bams <- bam_df$path
ref_fasta <- opt$ref_fasta

bed_df <- read.table(
    opt$capture_bed, header = FALSE, stringsAsFactors = FALSE, col.names = c("chromosome", "start", "end", "name")
)
bed_df$blank <- 0
bed_df$strand <- "+"
bed_df$exome_depth_roi <- "exome_depth_roi"
bed_df$gene_name <- paste(bed_df$name, "test", sep = "_")


# Counts for all bams
bam_counts <- getBamCounts(
    bed.frame = bed_df, bam.files = normal_bams, min.mapq = opt$min_mapq, referenceFasta = ref_fasta
)

# create count matrix
cohort_count_df <- as(bam_counts[, colnames(bam_counts)], "data.frame")
cohort_count_df$chromosome <- gsub(as.character(cohort_count_df$space), pattern = "chr", replacement = "")
cohort_count_matrix <- as.matrix(cohort_count_df[, grep(names(cohort_count_df), pattern = "*.bam")])

# save data
save(normal_bams, bed_df, ref_fasta, cohort_count_df, cohort_count_matrix, file = opt$out_path)
