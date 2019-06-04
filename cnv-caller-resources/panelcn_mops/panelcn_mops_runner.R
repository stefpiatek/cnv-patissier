library(dplyr)
library(optparse)
library(panelcn.mops)
library(purrr)


detect_cnvs <- function(sample_number) {
  "
    Function to detect CNV for each column of unknown Granges object
    Takes sample_number 
    Returns table of only positive CNVs with sample name
  "
  # Set up test data
  test_and_control <- unknowns[, sample_number]
  elementMetadata(test_and_control) <- cbind(elementMetadata(test_and_control), elementMetadata(normals))
  # Detect CNVs
  result_list <- runPanelcnMops(test_and_control, countWindows = count_windows, selectedGenes = selected_gene, maxControls = 30)
  # Output results from S4 object to table
  file_names <- colnames(elementMetadata(unknowns))
  raw_results <- createResultTable(
    resultlist = result_list, XandCB = test_and_control,
    countWindows = count_windows,
    selectedGenes = selected_gene,
    sampleNames = file_names
  )
  
  # Process table for output
  cnv_results <- raw_results[[1]] %>%
    filter(CN != "CN2") %>%
    mutate(Sample = unknown_samples$sample_name[sample_number]) %>%
    group_by(Sample, Chr, Gene, CN) %>%
    summarise(Exon = paste(min(as.character(Exon)), max(as.character(Exon))),
              Start = min(Start, End),
              End = max(Start, End),
              RC = mean(RC),
              medRC = mean(medRC),
              `RC.norm` = mean(`RC.norm`),
              `medRC.norm`  = mean(`medRC.norm`),
              lowQual = min(as.character(lowQual))
              ) %>%
    return()
}


# Parse options
option_list <- list(
  make_option(c("--output-path"), type = "character"),
  make_option(c("--gene"), type = "character"),
  make_option(c("--chrom-prefix"), type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)

# Set up parameters
bed <- file.path(opt$output_path, "capture.bed")
split_bed <- file.path(opt$output_path, "capture_split.bed")
splitROIs(bed, split_bed)

all_samples <- read.table(file.path(opt$output_path, "samples.tsv"), header = TRUE, stringsAsFactors = FALSE)
unknown_samples <- all_samples[all_samples$sample_type == "unknown", ]
normal_samples <- all_samples[all_samples$sample_type == "normal_panel", ]
selected_gene <- opt$gene

# Get counts
count_windows <- getWindows(split_bed)
normals <- countBamListInGRanges(countWindows = count_windows, bam.files = normal_samples$bam_path, read.width = 150)
unknowns <- countBamListInGRanges(countWindows = count_windows, bam.files = unknown_samples$bam_path, read.width = 150)

# Run CNV calling
called_cnvs <- map_dfr(.x = seq(1, nrow(unknown_samples)), .f = detect_cnvs) %>%
  filter(!is.na(Sample))

# Save results
write.table(called_cnvs, file = file.path(opt$output_path, "calls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
