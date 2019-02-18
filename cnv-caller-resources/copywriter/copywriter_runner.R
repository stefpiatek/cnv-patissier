library(CopywriteR)
library(optparse)
library(dplyr)
library(tidyr)
library(stringr)

option_list <- list(
    make_option(c('--max-cpu'), type = "numeric", default = NULL),
    make_option(c('--output-path'), type = "character", default = NULL),
    make_option(c('--capture-regions'), type = "character", default = NULL),
    make_option(c('--chromosome-prefix'), type = "character", default = NULL)

)

opt <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)

# Set up parameters
all_samples <- read.table(file.path(opt$output_path, "all_samples.txt"), header = TRUE, stringsAsFactors = FALSE)
bp.param <- SnowParam(workers = opt$max_cpu, type = "SOCK")

# Create proprocessed reference
preCopywriteR(
    output.folder = file.path(opt$output_path),
    bin.size = 50000,
    ref.genome = "hg19",
    prefix = opt$chromosome_prefix
)


if(opt$chromosome_prefix == "chr"){
    reference_name <- "hg19_50kb_chr"
} else {
    reference_name <- "hg19_50kb"
}


# Calculate read depth
CopywriteR(
    sample.control = all_samples,
    destination.folder = file.path(opt$output_path),
    reference.folder = file.path(opt$output_path, reference_name),
    capture.regions.file = opt$capture_regions,
    bp.param = bp.param,
    keep.intermediary.files = FALSE
)

# segment
plotCNA(destination.folder = file.path(opt$output_path), sample.plot = all_samples)

# Output segment results
load(file.path(opt$output_path, "CNAprofiles", "segment.Rdata"))
result_file <- file.path(opt$output_path, "results.txt")

segment.CNA.object$output %>%
    # remove self-comparisons and split out columns to have unknown and control used
    filter(seg.mean != 0) %>%  
    mutate(
        ID = str_replace_all(ID, "log2.", ""),
        start = round(loc.start),
        end = round(loc.end)
    ) %>%
    select(-c(loc.start, loc.end)) %>%
    separate(ID, c("unknown", "control"), ".vs.") %>%
    write.table(file = result_file, sep = "\t", row.names = FALSE, quote = FALSE)
