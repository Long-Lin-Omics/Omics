#!/usr/bin/env Rscript

# Load necessary library
library(rtracklayer)
library(data.table)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if input and output are provided
if (length(args) < 3) {
    stop("Usage: Rscript convert_peaks_to_bigbed.R input_peaks.bed output.bb fa.fai")
}

# Assign input and output file paths
peak_file <- args[1]
bb_file <- args[2]
chrom_sizes_file <- args[3]
chrom_sizes <- fread(chrom_sizes_file, header = FALSE, col.names = c("chrom", "size","offset","linebases","linewidth"))

# Import peak file (BED/narrowPeak format)
col_names <- c("chrom", "start", "end", "name", "score", "strand",
               "signalValue", "pValue", "qValue", "peak")

peak_data <- fread(peak_file, col.names = col_names, na.strings = ".", fill = TRUE)

gr <- GRanges(seqnames = peak_data$chrom,
              ranges = IRanges(start = peak_data$start, end = peak_data$end),
              strand = "*",
              score = peak_data$score)

chrom_sizes_filtered <- chrom_sizes[chrom_sizes$chrom %in% seqlevels(gr), ]
chrom_sizes_filtered <- chrom_sizes_filtered[match(seqlevels(gr), chrom_sizes_filtered$chrom), ]

seqlengths(gr) <- setNames(chrom_sizes_filtered$size, chrom_sizes_filtered$chrom)


# Export as BigBed
export(gr, bb_file, format = "BigBed")

cat("Conversion complete: ", peak_file, " → ", bb_file, "\n")

