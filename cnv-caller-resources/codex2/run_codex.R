library(optparse)
library(CODEX2)

option_list <- list(
    make_option(c("--output-path"), type = "character", default = NULL),
    make_option(c("--capture-bed"), type = "character", default = NULL),
    make_option(c("--chrom-prefix"), type = "character", default = NULL)

)

opt <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)

# Set up parameters
bed_df <- read.table(
    opt$capture_bed, header = FALSE, stringsAsFactors = FALSE, col.names = c("chromosome", "start", "end", "name")
)
all_samples <- read.table(file.path(opt$output_path, "samples.tsv"), header = TRUE, stringsAsFactors = FALSE)

normal_panel_index <- as.numeric(rownames(all_samples)[all_samples$sample_type == "normal_panel"])

bams_object <- getbambed(
    bamdir = all_samples$bam_path, 
    bedFile = opt$capture_bed, 
    sampname = all_samples$sample_name,
    projectname = "cnv-pat"
)
ref <- bams_object$ref

# GC and mappability
if(opt$chrom_prefix == "chr"){
    genome <- BSgenome.Hsapiens.UCSC.hg19
    mapp <- getmapp(ref, genome)
} else {
    # Use hg19 mappability scores if grch37
    hg19_genome <- BSgenome.Hsapiens.UCSC.hg19
    mapp <- getmapp(ref, hg19_genome)
    library(BSgenome.Hsapiens.1000genomes.hs37d5)
    genome <- BSgenome.Hsapiens.1000genomes.hs37d5
}

gc <- getgc(ref, genome)
gene <- bed_df$name
values(ref) <- cbind(values(ref), DataFrame(gc, mapp, gene))  

# coverage
coverage_object <- getcoverage(bams_object, mapqthres = 20)
Y <- coverage_object$Y
write.csv(Y, file = file.path(opt$output_path, "coverage.csv"), quote = FALSE)

# qc
qc_object <- qc(
    Y, all_samples$sample_name, ref, cov_thresh = c(20, Inf), length_thresh = c(20, Inf), 
    mapp_thresh = 0.9, gc_thresh = c(20, 80)
)

Y_qc <- qc_object$Y_qc 
sampname_qc <- qc_object$sampname_qc
ref_qc <- qc_object$ref_qc 
qcmat <- qc_object$qcmat
gc_qc <- ref_qc$gc

write.table(qcmat, file = file.path(opt$output_path, "qcmat.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


# Library size estimate for each sample
Y_non_zero <- Y_qc[apply(Y_qc, 1, function(x) {!any(x == 0)}), ]
pseudo_sample <- apply(Y_non_zero, 1, function(x) {prod(x ^ (1 / length(x)))})  # changed to avoid product being maximum value r can take (> 1.7e308)
N <- apply(apply(Y_non_zero, 2, function(x) {x / pseudo_sample}), 2, median)

# Reference sample normalisation
norm_object <- normalize_codex2_ns(Y_qc = Y_qc, gc_qc = gc_qc, K = 1:10, norm_index = normal_panel_index, N = N)


Yhat <- norm_object$Yhat
AIC <- norm_object$AIC
BIC <- norm_object$BIC
RSS <- norm_object$RSS


# CBS segmentation per gene for targeted seq
source("/mnt/cnv-caller-resources/codex2/segment_targeted.R")
opt_K <- which.max(BIC)
final_call <- matrix(ncol = 14, nrow = 0)
colnames(final_call) = c(
    "sample_name","chr","gene","cnv", "st_bp","ed_bp","length_kb", "st_exon","ed_exon","raw_cov",
    "norm_cov","copy_no","lratio","mBIC"
)

for(genei in unique(ref_qc$gene)){
  cat("Segmenting gene", genei,"\n")
  gene_index = which(ref_qc$gene == genei)
  yi <- Y_qc[gene_index, , drop=FALSE]
  yhati <- Yhat[[opt_K]][gene_index, , drop=FALSE]
  refi <- ref_qc[gene_index]
  final_calli <- segment_targeted(yi, yhati, sampname_qc, refi, genei, lmax = length(gene_index), mode = "fraction") 
  final_call <- rbind(final_call, final_calli)
}

cn <- (as.numeric(as.matrix(final_call[, "copy_no"])))

# removing calls with fractional copy numbers close to 2
cn_filter <- (cn <= 1.7) | (cn >=2.3) 

final_call  <- final_call[cn_filter, ]
if(class(final_call) == "character"){
  # if only 1 CNV, becomes a named character vector so convert to df
  final_call <- data.frame(as.list(final_call))  
}
length_exon <- as.numeric(final_call[, "ed_exon"])-as.numeric(final_call[,"st_exon"])+1
final_call <- cbind(final_call[, 1:7], length_exon, final_call[, 10:14])

write.table(final_call, file = file.path(opt$output_path, "calls.txt"), sep="\t", quote = FALSE, row.names = FALSE)
