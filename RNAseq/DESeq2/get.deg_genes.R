#!/usr/bin/env Rscript

# Usage:
#   Rscript merge_deg.R --mode union --padj 0.05 --lfc 1 --direction both \
#       res1.tsv res2.tsv res3.tsv -o merged_genes.txt
#
# Options:
#   --mode       "union" or "intersect"
#   --padj       padj cutoff (default: 0.05)
#   --lfc        log2FoldChange cutoff (absolute value, default: 0)
#   --direction  "both", "up", or "down" [default both]
#   -o           output file (default: merged_genes.txt)

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--mode"), type="character", default="union",
              help="merge mode: union or intersect [default %default]"),
  make_option(c("--padj"), type="double", default=0.05,
              help="padj cutoff [default %default]"),
  make_option(c("--lfc"), type="double", default=0,
              help="absolute log2FoldChange cutoff [default %default]"),
  make_option(c("--direction"), type="character", default="both",
              help="direction: both, up, or down [default %default]"),
  make_option(c("-o", "--out"), type="character", default="merged_genes.txt",
              help="output file name [default %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, positional_arguments = TRUE)

mode <- tolower(opt$options$mode)
padj_cutoff <- opt$options$padj
lfc_cutoff <- opt$options$lfc
direction <- tolower(opt$options$direction)
out_file <- opt$options$out
files <- opt$args

if (length(files) < 2) {
  stop("Please provide at least two DESeq2 result files.")
}

gene_lists <- lapply(files, function(f) {
  df <- read.delim(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  # detect gene column
  if (!"gene" %in% colnames(df)) {
    if ("transcript" %in% colnames(df)){
      df$gene <- df$transcript
    }
    else if ("X" %in% colnames(df)) {
      df$gene <- df$X
    } else if ("row.names" %in% colnames(df)) {
      df$gene <- df$row.names
    } else {
      df$gene <- rownames(df)
    }
  }
  
  if (!"padj" %in% colnames(df) | !"log2FoldChange" %in% colnames(df)) {
    stop(paste("File", f, "must contain columns: padj, log2FoldChange"))
  }
  
  # filtering
  df_filt <- subset(df, !is.na(padj) & padj < padj_cutoff & abs(log2FoldChange) >= lfc_cutoff)
  
  if (direction == "up") {
    df_filt <- subset(df_filt, log2FoldChange >= lfc_cutoff)
  } else if (direction == "down") {
    df_filt <- subset(df_filt, log2FoldChange <= -lfc_cutoff)
  }
  
  unique(df_filt$gene)
})

if (mode == "union") {
  merged_genes <- Reduce(union, gene_lists)
} else if (mode == "intersect") {
  merged_genes <- Reduce(intersect, gene_lists)
} else {
  stop("Invalid mode. Use 'union' or 'intersect'.")
}

write.table(merged_genes, file=out_file,
            quote=FALSE, row.names=FALSE, col.names=FALSE)

cat(sprintf("Done. %d genes written to %s\n", length(merged_genes), out_file))
