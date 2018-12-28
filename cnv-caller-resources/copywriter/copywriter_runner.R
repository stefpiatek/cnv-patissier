library(CopywriteR)
library(optparse)

option_list <- list(
    make_option(c('--max-cpu'), type = "numeric", default = NULL),
    make_option(c('--output-path'), type = "character", default = NULL),
    make_option(c('--capture-regions'), type = "character", default = NULL)
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
    prefix = "chr"
)

# Calculate read depth
CopywriteR(
    sample.control = all_samples,
    destination.folder = file.path(opt$output_path),
    reference.folder = file.path(opt$output_path, "hg19_50kb_chr"),
    capture.regions.file = opt$capture_regions,
    bp.param = bp.param,
    keep.intermediary.files = FALSE
)

# segment
plotCNA(destination.folder = file.path(opt$output_path), sample.plot = all_samples)